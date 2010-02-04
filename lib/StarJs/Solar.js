/*global StarJs: false */
/*jslint onevar: false */
/**
 * @depends ../StarJs.js
 * @depends Math.js
 * @depends Vector.js
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
    var PI0_5 = Math.PI/2, PI1_5 = 1.5*Math.PI;

    function sinAlt(approxBody, mjd, lambda, cphi, sphi) {
        var t, pos, tau;
        t = (mjd - 51544.5) / 36525.0;
        pos = approxBody(t);
        tau = StarJs.Time.gmst(mjd) + lambda - pos.ra;
        return {
            sa: sphi * Math.sin(pos.dec) +
                cphi * Math.cos(pos.dec) * Math.cos(tau),
            pos: pos
        };
    }

    function phase(sun, moon) {
        var TAU_SUN = 8.32 / (1440.0 * 36535.0);
        var diff = moon.phi - sunpos.phi; // We ignore speed of light here
        return StarJs.Math.mod(diff + Math.PI, StarJs.Math.PI2) - Math.PI;
    }

    var EVENTS = StarJs.Solar.EVENTS;

    var result = [], today, h, j, ptime, ctime, ntime, H = 1.0 / 24.0, posp = {}, pos0 = {}, posn = {}, cphi = Math.cos(phi), sphi = Math.sin(phi), name, evt, yp = {sun: {}, moon: {}}, y0 = {sun: {}, moon: {}}, yn = {sun: {}, moon: {}}, interp = {}, php, ph0, phn;
    if (typeof tz === 'undefined') {
        tz = function (a) {
            return a;
        };
    }

    ptime = start;

    posp.moon = sinAlt(StarJs.Solar.approxMoon, ptime, lambda, cphi, sphi);
    posp.sun  = sinAlt(StarJs.Solar.approxSun,  ptime, lambda, cphi, sphi);
    php = phase(posp.sun.pos, posp.moon.pos);

    for (j = 0; j < EVENTS.length; j += 1) {
        evt = EVENTS[j];
        name = evt.name;

        yp[evt.body][name] = posp[evt.body].sa - evt.sinh0;
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

            // Calculate phases
            ph0 = phase(pos0.sun.pos, pos0.moon.pos);
            phn = phase(posn.sun.pos, posn.moon.pos);

            for (j = 0; j < EVENTS.length; j += 1) {
                evt = EVENTS[j];
                name = evt.name;
                today[evt.body][name] = today[evt.body][name] || {};

                y0[evt.body][name] = pos0[evt.body].sa - evt.sinh0;
                yn[evt.body][name] = posn[evt.body].sa - evt.sinh0;

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


            // TODO: convert to [-PI,PI) range!!!  Or something more
            // difficult is required?
            
            // New moon
            interp = StarJs.Math.quadInterpolation(php, ph0, phn);
            // First quarter
            interp = StarJs.Math.quadInterpolation(php-PI0_5, ph0-PI0_5, phn-PI0_5);
            // Full moon
            interp = StarJs.Math.quadInterpolation(php-Math.PI, ph0-Math.PI, phn-Math.PI);

            // Third quarter
            interp = StarJs.Math.quadInterpolation(php-PI1_5, ph0-PI1_5, phn-PI1_5);


            // Prepare for next interval
            yp = yn;
            yn = {moon: {}, sun: {}};

            ptime = ntime;
            php = phn;
        }

        for (j = 0; j < EVENTS.length; j += 1) {
            evt = EVENTS[j];
            name = evt.name;

            if (!today[evt.body][name].set && !today[evt.body][name].rise) {
                today[evt.body][name].alwaysAbove = (pos0[evt.body].sa > evt.sinh0);
            }
        }
        // Next day
        ptime = (start += 1.0);

        result.push(today);
    }

    return result;
};
