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

    l0 = StarJs.Math.frac(0.606433 + 1336.85525 * mct);
    l  = StarJs.Math.PI2 * StarJs.Math.frac(0.374897 + 1325.552410 * mct);
    ls = StarJs.Math.PI2 * StarJs.Math.frac(0.993133 + 99.997361 * mct);
    d  = StarJs.Math.PI2 * StarJs.Math.frac(0.827361 + 1236.853086 * mct);
    f  = StarJs.Math.PI2 * StarJs.Math.frac(0.259086 + 1342.227825 * mct);

    dl = +22640 * Math.sin(l) - 4686 * Math.sin(l - 2 * d) + 2370 * Math.sin(2 * d) +
        769 * Math.sin(2 * l) - 668 * Math.sin(ls) - 412 * Math.sin(2 * f) +
        -212 * Math.sin(2 * l - 2 * d) - 206 * Math.sin(l + ls - 2 * d) +
        192 * Math.sin(l + 2 * d) + 192 * Math.sin(l + 2 * d) - 165 * Math.sin(ls - 2 * d) +
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
    name: 'moon',
    title: "Moon",
    sinh0: Math.sin((+8.0 / 60.0) * StarJs.Math.DEG2RAD),
    body: 'moon',
    posFunc: StarJs.Solar.approxMoon
}, {
    name: 'sun',
    title: "Sun",
    sinh0: Math.sin((-50.0 / 60.0) * StarJs.Math.DEG2RAD),
    body: 'sun',
    posFunc: StarJs.Solar.approxSun
}, {
    name: 'twilightC',
    title: "Civil twilight",
    sinh0: Math.sin(-6.0 * StarJs.Math.DEG2RAD),
    body: 'sun',
    posFunc: StarJs.Solar.approxSun
}, {
    name: 'twilightN',
    title: "Nautical twilight",
    sinh0: Math.sin(-12.0 * StarJs.Math.DEG2RAD),
    body: 'sun',
    posFunc: StarJs.Solar.approxSun
}, {
    name: 'twilightA',
    title: "Astronomical twilight",
    sinh0: Math.sin(-18.0 * StarJs.Math.DEG2RAD),
    body: 'sun',
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
 * {
 *   moon: { rise: 22.3 },
 *   sun:  { rise: 4.5, set: 21.89 },
 *   twilightA: { rise: 3.10, set: 23.6 }
 * }
 *
 * (This is a purely artificial example, I haven't checked its correctnes).
 *
 *  'moon' and 'sun' objects also have property 'midnight' that gives
 *  sinAlt value of Moon and Sun at beginning of the day (so you can
 *  tell if body never rises or never sets this day if no event was
 *  registered).
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

    var result = [], today, h, j, ptime, ctime, ntime, H = 1.0 / 24.0, posp = {}, pos0 = {}, posn = {}, cphi = Math.cos(phi), sphi = Math.sin(phi), name, evt, yp = {}, y0 = {}, yn = {}, interp = {}, midmoon, midsun;
    if (typeof tz === 'undefined') {
        tz = function (a) {
            return a;
        };
    }

    // This is a sketchy implementation that sure doesn't work.
    ptime = start;

    posp.moon = sinAlt(StarJs.Solar.approxMoon, ptime, lambda, cphi, sphi);
    posp.sun  = sinAlt(StarJs.Solar.approxSun,  ptime, lambda, cphi, sphi);

    for (j = 0; j < EVENTS.length; j += 1) {
        evt = EVENTS[j];
        name = evt.name;

        yp[name] = posp[evt.body] - evt.sinh0;
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
                today[name] = today[name] || {};

                y0[name] = pos0[evt.body] - evt.sinh0;
                yn[name] = posn[evt.body] - evt.sinh0;

                interp = StarJs.Math.quadInterpolation(yp[name], y0[name], yn[name]);

                switch (interp.roots.length) {
                case 0:
                    // No roots
                    break;
                case 1:
                    // Single root
                    if (yp[name] < 0.0) {
                        today[name].rise = h + interp.roots[0];
                    } else {
                        today[name].set  = h + interp.roots[0];
                    }
                    break;
                case 2:
                    // Two roots
                    if (interp.ye < 0.0) {
                        today[name].rise = h + interp.roots[1];
                        today[name].set  = h + interp.roots[0];
                    } else {
                        today[name].rise = h + interp.roots[0];
                        today[name].set  = h + interp.roots[1];
                    }
                    break;
                }
            }

            // Next interval
            yp = yn;
            yn = {};

            ptime = ntime;
        }
        // Next day
        ptime = (start += 1.0);

        result.push(today);
    }

    return result;
};
