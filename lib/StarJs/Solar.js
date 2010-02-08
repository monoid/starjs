/*global StarJs: false */
/*jslint onevar: false */
/**
 * @depends ../StarJs.js
 * @depends Math.js
 * @depends Vector.js
 * @depends Kepler.js
 */
StarJs.Solar = (typeof StarJs.Solar === 'undefined') ? {} : StarJs.Solar;

StarJs.Solar.EPS = 23.43920111 * StarJs.Math.DEG2RAD;

StarJs.Solar.approxMoon = function (mct) {
    var l0, l, ls, f, d, dl, s, h, n, ml, mb, me;

    l0 = StarJs.Math.frac(0.606433 + 1336.855225 * mct);
    l  = StarJs.Math.PI2 * StarJs.Math.frac(0.374897 + 1325.552410 * mct);
    ls = StarJs.Math.PI2 * StarJs.Math.frac(0.993133 + 99.997361 * mct);
    d  = StarJs.Math.PI2 * StarJs.Math.frac(0.827361 + 1236.853086 * mct);
    f  = StarJs.Math.PI2 * StarJs.Math.frac(0.259086 + 1342.227825 * mct);

    dl = +22640 * Math.sin(l) - 4586 * Math.sin(l - 2 * d) + 2370 * Math.sin(2 * d) +
        769 * Math.sin(2 * l) - 668 * Math.sin(ls) - 412 * Math.sin(2 * f) +
        -212 * Math.sin(2 * l - 2 * d) - 206 * Math.sin(l + ls - 2 * d) +
        192 * Math.sin(l + 2 * d) - 165 * Math.sin(ls - 2 * d) +
        -125 * Math.sin(d) - 110 * Math.sin(l + ls) + 148 * Math.sin(l - ls) +
        -55 * Math.sin(2 * f - 2 * d);
    s = f + (dl + 412 * Math.sin(2 * f) + 541 * Math.sin(ls)) * StarJs.Math.ARCS2RAD;
    h = f - 2 * d;
    n = -526 * Math.sin(h) + 44 * Math.sin(l + h) - 31 * Math.sin(h - l) +
        -23 * Math.sin(ls + h) + 11 * Math.sin(h - ls) - 25 * Math.sin(f - 2 * l) +
        21 * Math.sin(f - l);
    ml = StarJs.Math.PI2 * StarJs.Math.frac(l0 + dl / 1296.0e3);
    mb = (18520.0 * Math.sin(s) + n) * StarJs.Math.ARCS2RAD;
    me = StarJs.Vector.Matrix3.r_x(-StarJs.Solar.EPS).apply((new StarJs.Vector.Polar3(ml, mb)).toVector3());
    me = new StarJs.Vector.Polar3(me);
    return {ra: me.phi, dec: me.theta};

};

StarJs.Solar.approxSun = function (mct) {
    var l, m, m2, me, se;
    m2 = StarJs.Math.frac(0.993133 + 99.997361 * mct);
    m = StarJs.Math.PI2 * m2;
    l = StarJs.Math.PI2 * StarJs.Math.frac(
        0.7859453  + m2 + (6893.0 * Math.sin(m) +
                           72.0 * Math.sin(2 * m) +
                           6191.2 * mct) / 1296.0e3);
    me = StarJs.Vector.Matrix3.r_x(-StarJs.Solar.EPS).apply((new StarJs.Vector.Polar3(l, 0)).toVector3());
    me = new StarJs.Vector.Polar3(me);
    return {ra: me.phi, dec: me.theta};
};


StarJs.Solar.EVENTS = [{
    body: 'moon',
    name: 'day',
    title: "Moon",
    sinh0: Math.sin((+8.0 / 60.0) * StarJs.Math.DEG2RAD),
    posFunc: StarJs.Solar.approxMoon
}, {
    body: 'sun',
    name: 'day',
    title: "Sun",
    sinh0: Math.sin((-50.0 / 60.0) * StarJs.Math.DEG2RAD),
    posFunc: StarJs.Solar.approxSun
}, {
    body: 'sun',
    name: 'twilightC',
    title: "Civil twilight",
    sinh0: Math.sin(-6.0 * StarJs.Math.DEG2RAD),
    posFunc: StarJs.Solar.approxSun
}, {
    body: 'sun',
    name: 'twilightN',
    title: "Nautical twilight",
    sinh0: Math.sin(-12.0 * StarJs.Math.DEG2RAD),
    posFunc: StarJs.Solar.approxSun
}, {
    body: 'sun',
    name: 'twilightA',
    title: "Astronomical twilight",
    sinh0: Math.sin(-18.0 * StarJs.Math.DEG2RAD),
    posFunc: StarJs.Solar.approxSun
}];

/** Calculate Sun and Moon events: rise, set and twilight.
 *
 * start: start time (MJD)
 * end: end day (MJD)
 * lon, lat: geographical coordinates
 * tz: timezone offset function (converts UTC to local).
 *
 * Returns array of objects, each object describes particular day in form:
 *
 *  {
 *    moon: { 
 *      midnight: -0.4575,
 *      day: {
 *          rise: 10.2689,
 *          set: 21.0516
 *      }
 *    },
 *    sun:  { 
 *      midnight: -0.3877,
 *      day:       { rise: 6.3015, set: 20.56 },
 *      twilightA: { rise: 3.9416, set: 22.9039 },
 *      twilightN: { ... },
 *      twilightC: { ... }
 *    }
 *  }
 *
 * 'moon' and 'sun' objects have property 'midnight' that gives sinAlt
 * value of Moon and Sun at beginning of the day.  If rise and set
 * values absent, boolean 'alwaysAbove' field is set.
 */
StarJs.Solar.sunAndMoonEvents = function (start, end, lambda, phi, tz) {
    function sinAlt(approxBody, mjd, lambda, cphi, sphi) {
        var t, pos, tau;
        t = (mjd - 51544.5) / 36525.0;
        pos = approxBody(t);
        tau = StarJs.Time.gmst(mjd) + lambda - pos.ra;
        return sphi * Math.sin(pos.dec) +
            cphi * Math.cos(pos.dec) * Math.cos(tau);
    }

    var EVENTS = StarJs.Solar.EVENTS;

    var result = [], today, h, j, ptime, ctime, ntime, H = 1.0 / 24.0, posp = {}, pos0 = {}, posn = {}, cphi = Math.cos(phi), sphi = Math.sin(phi), name, evt, yp = {sun: {}, moon: {}}, y0 = {sun: {}, moon: {}}, yn = {sun: {}, moon: {}}, interp = {};
    if (typeof tz === 'undefined') {
        tz = function (a) {
            return a;
        };
    }

    ptime = start;

    posp.moon = sinAlt(StarJs.Solar.approxMoon, ptime, lambda, cphi, sphi);
    posp.sun  = sinAlt(StarJs.Solar.approxSun,  ptime, lambda, cphi, sphi);

    for (j = 0; j < EVENTS.length; j += 1) {
        evt = EVENTS[j];
        name = evt.name;

        yp[evt.body][name] = posp[evt.body] - evt.sinh0;
    }

    while (ptime < end) {
        today = {
            moon: {
                midnight: posp.moon
            },
            sun: {
                midnight: posp.sun
            }
        };
        for (h = 1; h < 24; h += 2) {
            ctime = ptime + H;
            ntime = ctime + H; // ntime = ctime + 2 * H;

            // Calc new positions...
            pos0.moon = sinAlt(StarJs.Solar.approxMoon, ctime, lambda,
                               cphi, sphi);
            pos0.sun = sinAlt(StarJs.Solar.approxSun, ctime, lambda,
                              cphi, sphi);
            posn.moon = sinAlt(StarJs.Solar.approxMoon, ntime, lambda,
                               cphi, sphi);
            posn.sun = sinAlt(StarJs.Solar.approxSun, ntime, lambda,
                              cphi, sphi);

            for (j = 0; j < EVENTS.length; j += 1) {
                evt = EVENTS[j];
                name = evt.name;
                today[evt.body][name] = today[evt.body][name] || {};

                y0[evt.body][name] = pos0[evt.body] - evt.sinh0;
                yn[evt.body][name] = posn[evt.body] - evt.sinh0;

                interp = StarJs.Math.quadInterpolation(yp[evt.body][name], y0[evt.body][name], yn[evt.body][name]);

                switch (interp.roots.length) {
                case 0:
                    // No roots
                    break;
                case 1:
                    // Single root
                    if (yp[evt.body][name] < 0.0) {
                        today[evt.body][name].rise = h + interp.roots[0];
                    } else {
                        today[evt.body][name].set  = h + interp.roots[0];
                    }
                    break;
                case 2:
                    // Two roots
                    if (interp.ye < 0.0) {
                        today[evt.body][name].rise = h + interp.roots[1];
                        today[evt.body][name].set  = h + interp.roots[0];
                    } else {
                        today[evt.body][name].rise = h + interp.roots[0];
                        today[evt.body][name].set  = h + interp.roots[1];
                    }
                    break;
                }
            }

            // Next interval
            yp = yn;
            yn = {moon: {}, sun: {}};

            ptime = ntime;
        }

        for (j = 0; j < EVENTS.length; j += 1) {
            evt = EVENTS[j];
            name = evt.name;

            if (!today[evt.body][name].set && !today[evt.body][name].rise) {
                today[evt.body][name].alwaysAbove = (pos0[evt.body] > evt.sinh0);
            }
        }
        // Next day
        ptime = (start += 1.0);

        result.push(today);
    }

    return result;
};


/**********************************************************************
 *
 * Solar sytem plantes.
 *
 */

StarJs.Solar.Sun = function () {
    this.name = 'Sun';
};

StarJs.Solar.Sun.$POS = new StarJs.Vector.Vector3(0, 0, 0);

StarJs.Solar.GM = StarJs.Math.sqr(0.01720209895);

StarJs.Solar.Sun.prototype.keplerCoord = function (t) {
    return StarJs.Solar.Sun.$POS;
};

StarJs.Solar.Body = function (name, params) {
    this.name = name;
    this.a  = params.a;
    this.ec = params.e;
    this.m0 = params.M0;
    this.n  = params.n;
    this.omega = params.O;
    this.i  = params.i;
    this.w  = params.w;
    this.t0 = params.T0;
    // ...
};

StarJs.Solar.Body.prototype.keplerCoord = function (t) {
    var P = 1.3970;
    var m = StarJs.Math.DEG2RAD * (this.m0 + this.n * (t - this.t0));
    var r = StarJs.Kepler.elliptic(StarJs.Solar.GM, m, this.a, this.ec);
    var pqr = StarJs.Kepler.gaussVec(
        StarJs.Math.DEG2RAD * (this.omega + P * t),
        StarJs.Math.DEG2RAD * this.i,
        StarJs.Math.DEG2RAD * (this.w - this.omega));
    return pqr.apply(r);
};

StarJs.Solar.BODIES = {
    'Sun':     new StarJs.Solar.Sun(),
    'Mercury': new StarJs.Solar.Body('Mercury', {
      a :  0.387099, e : 0.205634, M0 : 174.7947, n : 149472.6738,
      O :  48.331,   i : 7.0048,   w  :  77.4552, T0 : 0.0
    }),
    'Venus':   new StarJs.Solar.Body('Venus', {
      a :  0.723332, e : 0.006773, M0 :  50.4071, n  : 58517.8149,
      O :  76.680,   i : 3.3946,   w  : 131.5718, T0 : 0.0
    }),
    'Earth':   new StarJs.Solar.Body('Earth', {
      a :  1.000000, e : 0.016709, M0 : 357.5256, n  : 35999.3720,
      O : 174.876,   i : 0.0000,   w  : 102.9400, T0 : 0.0
    }),
    'Mars':   new StarJs.Solar.Body('Mars', {
      a :  1.523692, e : 0.093405, M0 :  19.3879, n  : 19140.3023,
      O :  49.557,   i : 1.8496,   w  : 336.0590, T0 : 0.0
    }),
    'Jupiter':   new StarJs.Solar.Body('Jupiter', {
      a :  5.204267, e : 0.048775, M0 :  18.8185, n  :  3033.6272,
      O : 100.4908,  i : 1.3046,   w  :  15.5576, T0 : 0.0
    }),
    'Saturn':   new StarJs.Solar.Body('Saturn', {
      a :  9.582018, e : 0.055723, M0 : 320.3477, n  :  1213.8664,
      O : 113.6427,  i : 2.4852,   w  :  89.6567, T0 : 0.0
    }),
    'Uranus':   new StarJs.Solar.Body('Uranus', {
      a : 19.229412, e : 0.044406, M0 : 142.9559, n  :   426.9282,
      O :  73.9893,  i : 0.7726,   w  : 170.5310, T0 : 0.0
    }),
    'Neptune':   new StarJs.Solar.Body('Neptune', {
      a : 30.103658, e : 0.011214, M0 : 267.7649, n  :   217.9599,
      O : 131.7942,  i : 1.7680,   w  :  37.4435, T0 : 0.0
    })
};
