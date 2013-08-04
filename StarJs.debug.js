StarJs = {};
StarJs.Math = {};

/**
 * Scale factor for converting degrees to radians.
 * @const
 * @type {number}
 */
StarJs.Math.DEG2RAD = Math.PI / 180.0;
/**
 * Scale factor for converting radians to degrees.
 * @const
 * @type {number}
 */
StarJs.Math.RAD2DEG = 180.0 / Math.PI;
/**
 * Scale factor for converting radians to arcseconds.
 * @const
 * @type {number}
 */
StarJs.Math.RAD2ARCS = 648000.0 / Math.PI;
/**
 * Scale factor for converting arcseconds to radians.
 * @const
 * @type {number}
 */
StarJs.Math.ARCS2RAD = Math.PI / 648000.0;
/**
 * 2*PI
 * @const
 * @type {number}
 */
StarJs.Math.PI2 = 2.0 * Math.PI;

/**
 * Scale factor for converting radians to hour measure.
 * @const
 * @type {number}
 */
StarJs.Math.ANGLE2HMS = 12.0 / Math.PI;

/**
 * Return square of x (x*x)
 * @param {number} x argument
 * @return {number} x*x
 */
StarJs.Math.sqr = function (x) {
    return x * x;
};

/**
 * Return fractional part of argument
 * @param {number} x argument
 * @return {number} fractioal part
 */
StarJs.Math.frac = function (x) {
    return x - Math.floor(x);
};

/**
 * Modulo
 * @param {number} dividend
 * @param {number} r divisor
 * @return {number} remainder
 */
StarJs.Math.mod = function (x, r) {
    return r * StarJs.Math.frac(x / r);
};

/**
 * Degree-minute-second object.
 * Either 1 or 4 arguments are accepted.
 * @constructor
 * @param {number} sign_or_angle sign if four values are passed, angle
 *                 in degrees otherwise
 * @param {number} degree Integer degree part (optional)
 * @param {number} minute Integer minute part (optional)
 * @param {number} second Seconds (optional)
 */
StarJs.Math.Dms = function (sign_or_angle, degree, minute, second) {
    if (arguments.length === 4) {
        this['sign'] = sign_or_angle;
        this['degree'] = degree;
        this['minute'] = minute;
        this['second'] = second;
    } else {
        var angle = sign_or_angle;
        this['sign'] = 1;
        if (angle < 0) {
            angle = -angle;
            this['sign'] = -1;
        }
        this['degree'] = Math.floor(angle);
        angle = 60 * (angle - this['degree']);
        this['minute'] = Math.floor(angle);
        angle = 60 * (angle - this['minute']);
        this['second'] = angle;
    }
};

/**
 * Convert angle DMS (degree-minute-second) to float degree value.
 */
StarJs.Math.dms2deg = function (sign, deg, minute, second) {
    return sign * (deg + minute / 60.0 + second / 3600.0);
};

/**
 * Convert angle (degree-minute-second) to float degree value.
 */
StarJs.Math.Dms.prototype.dms2deg = function () {
    return StarJs.Math.dms2deg(this['sign'], this['degree'], this['minute'], this['second']);
};

/**
 * Convert float degree value to DMS (degree-minute-second).
 * @param {number} deg Angle
 */
StarJs.Math.Dms.deg2dms = function (deg) {
    return new StarJs.Math.Dms(deg);
};

/**
 * Convert radians to hour measure.
 * @param {number} angle Angle in radians.
 */
StarJs.Math.angle2hms = function (angle) {
    angle = StarJs.Math.mod(angle, StarJs.Math.PI2);
    var a = StarJs.Math.ANGLE2HMS * angle, res = StarJs.Math.deg2dms(a);
    // Change field names and remove sign field as it is always 1.
    res['hour'] = res['deg'];
    delete res['deg'];
    delete res['sign'];
    
    return res;
};

StarJs.Math.quadInterpolation = function (ym, y0, yp) {
    var a = 0.5 * (yp + ym) - y0, b = 0.5 * (yp - ym), c = y0, xe = -b / (2 * a), ye = (a * xe + b) * xe + c;
    var dis = b * b - 4 * a * c, roots = [], dx, r1, r2;
    if (dis >= 0) {
        dx = 0.5 * Math.sqrt(dis) / Math.abs(a);

        r1 = xe - dx;
        r2 = xe + dx;
        
        if (Math.abs(r1) <= 1.0) {
            roots.push(r1);
        }
        if (Math.abs(r2) <= 1.0) {
            roots.push(r2);
        }
    }
    return {'xe': xe, 'ye': ye, 'roots': roots};
};

/**
 * Hyperbolic sinus.
 */
StarJs.Math.sinh = function (x) {
    return (Math.exp(x) - Math.exp(-x))/2;
};

/**
 * Hyperbolic cosinus.
 */
StarJs.Math.cosh = function (x) {
    return (Math.exp(x) + Math.exp(-x))/2;
};

StarJs.Vector = {};

/** @constructor */
StarJs.Vector.Vector3 = function (x, y, z) {
    this['x'] = x;
    this['y'] = y;
    this['z'] = z;
};

(function (p) {
    p.len = function () {
        return Math.sqrt(this['x'] * this['x'] + this['y'] * this['y'] + this['z'] * this['z']);
    };

    p.add = function (v3) {
        return new StarJs.Vector.Vector3(this['x'] + v3['x'],
                                         this['y'] + v3['y'],
                                         this['z'] + v3['z']);
    };

    p.sub = function (v3) {
        return new StarJs.Vector.Vector3(this['x'] - v3['x'],
                                         this['y'] - v3['y'],
                                         this['z'] - v3['z']);
    };

    p.neg = function () {
        return new StarJs.Vector.Vector3(-this['x'], -this['y'], -this['z']);
    };

    p.scale = function (a) {
        return new StarJs.Vector.Vector3(a * this['x'], a * this['y'], a * this['z']);
    };

    p.clone = function () {
        return new StarJs.Vector.Vector3(this['x'], this['y'], this['z']);
    };
}(StarJs.Vector.Vector3.prototype));

/** @constructor */
StarJs.Vector.Polar3 = function (az_or_v3, elev, rad) {
    var alen = arguments.length;
    if (alen === 2) {
        rad = 1.0;
    }
    if (alen === 2 || alen === 3) {
        this['phi'] = az_or_v3;
        this['theta'] = elev;
        this['rad'] = rad;
    } else {
        var v3 = az_or_v3;
        var rho2 = v3['x'] * v3['x'] + v3['y'] * v3['y'];
        this['rad'] = Math.sqrt(rho2 + v3['z'] * v3['z']);
        this['phi'] = (v3['x'] === 0.0 && v3['y'] === 0.0) ?  0.0 : Math.atan2(v3['y'], v3['x']);
        if (this['phi'] < 0.0) {
            this['phi'] += StarJs.Math.PI2;
        }
        var rho = Math.sqrt(rho2);
        this['theta'] = (v3['z'] === 0.0 && rho === 0.0) ? 0.0 : Math.atan2(v3['z'], rho);
    }
};

(function (p) {
    p.toVector3 = function () {
        var ct = Math.cos(this['theta']);
        var rad = this['rad']
        var x = rad * ct * Math.cos(this['phi']);
        var y = rad * ct * Math.sin(this['phi']);
        var z = rad * Math.sin(this['theta']);
        return new StarJs.Vector.Vector3(x, y, z);
    };

    p.clone = function () {
        return new StarJs.Vector.Polar3(this['rad'], this['phi'], this['theta']);
    };
}(StarJs.Vector.Polar3.prototype));

/** @constructor */
StarJs.Vector.Matrix3 = function (v1, v2, v3) {
    if (arguments.length === 3) {
        this['mat'] = [[v1['x'], v1['y'], v1['z']], [v2['x'], v2['y'], v2['z']], [v3['x'], v3['y'], v3['z']]];
    } else {
        if (v2) {
            this['mat'] = v1;
        } else {
            this['mat'] = [[v1[0][0], v1[0][1], v1[0][2]],
                        [v1[1][0], v1[1][1], v1[1][2]],
                        [v1[2][0], v1[2][1], v1[2][2]]];
        }
    }
};

(function (p) {
    p.apply = function (v) {
        var l = this['mat'][0];
        var x = l[0] * v['x'] + l[1] * v['y'] + l[2] * v['z'];

        l = this['mat'][1];
        var y = l[0] * v['x'] + l[1] * v['y'] + l[2] * v['z'];

        l = this['mat'][2];
        var z = l[0] * v['x'] + l[1] * v['y'] + l[2] * v['z'];

        return new StarJs.Vector.Vector3(x, y, z);
    };

    p.mult = function (matrix) {
        var res = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        var mat = matrix['mat'];
        for (var i = 0; i < 3; i += 1) {
            var tline = this['mat'][i];
            var rline = res[i];
            for (var j = 0; j < 3; j += 1) {
                for (var k = 0; k < 3; k += 1) {
                    rline[j] += tline[k] * mat[k][j];
                }
            }
        }
        return new StarJs.Vector.Matrix3(res, true);
    };
}(StarJs.Vector.Matrix3.prototype));

StarJs.Vector.Matrix3.r_x = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[1.0, 0.0, 0.0],
                                      [0.0,  cp,  sp],
                                      [0.0, -sp,  cp]]);
};

StarJs.Vector.Matrix3.r_y = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[ cp, 0.0, -sp],
                                      [0.0, 1.0, 0.0],
                                      [ sp, 0.0,  cp]]);
};

StarJs.Vector.Matrix3.r_z = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[ cp,  sp, 0.0],
                                      [-sp,  cp, 0.0],
                                      [0.0, 0.0, 1.0]]);
    
};
StarJs.Time = {};

StarJs.Time['DEFAULT_JULIAN_DATE'] = {'year': 1582, 'month': 10, 'day': 4};
StarJs.Time['DEFAULT_JULIAN_JD'] = 2299161;
StarJs.Time['JD_MJD'] = 2400000.5;

StarJs.Time.time2mjd = function (t) {
    if (typeof t !== 'number') {
        t = t.getTime();
    }
    return t / 86400000.0 + 40587.0;
};

StarJs.Time.mjd2dt = function (mjd, jul) {
    if (typeof jul === 'undefined') {
        jul = StarJs.Time['DEFAULT_JULIAN_JD'];
    }
    var a = Math.floor(mjd) + 2400001, b, c, d, e, f, day, mon, year, hour, t;
    if (a < jul) {
        // Julian calendar
        b = 0;
        c = a + 1524;
    } else {
        // Grigorian calendar
        b = Math.floor((a - 1867216.25) / 36524.25);
        c = a + b - Math.floor(b / 4) + 1525;
    }
    d = Math.floor((c - 122.1) / 365.25);
    e = 365 * d + Math.floor(d / 4);
    f = Math.floor((c - e) / 30.6001);
    day = c - e - Math.floor(30.6001 * f);
    mon = f - 1 - 12 * Math.floor(f / 14);
    year = d - 4715 - Math.floor((7 + mon) / 10);
    hour = 24.0 * (mjd - Math.floor(mjd));
    t = StarJs.Time.hour2hms(hour);
    t['year'] = year;
    t['month'] = mon;
    t['day'] = day;
    return t;
};

StarJs.Time.hms2hour = function (h, m, s) {
    return h + (m / 60.0) + (s / 3600.0);
};

StarJs.Time.hour2hms = function (hour) {
    var dms = new StarJs.Math.Dms(hour);
    return {'hour': dms['degree'], 'minute': dms['minute'], 'second': dms['second']};
};

StarJs.Time.dt2mjd = function (dt, jul) {
    if (typeof jul === 'undefined') {
        jul = StarJs.Time['DEFAULT_JULIAN_DATE'];
    }

    var year = dt['year'], mon = dt['month'], b;
    if (mon <= 2) {
        mon += 12;
        year -= 1;
    }
    if (year <= jul['year'] && mon <= jul['month'] && dt['day'] <= jul['day']) {
        // Julian
        b = -2 + Math.floor((year + 4716) / 4) - 1179;
    } else {
        // Gregorian
        b = Math.floor(year / 400) - Math.floor(year / 100) + Math.floor(year / 4);
    }
    var mjdMidnight = 365 * year - 679004 + b + Math.floor(30.6001 * (mon + 1)) + dt['day'];
    var frac = StarJs.Time.hms2hour(dt['hour'], dt['minute'], dt['second']) / 24.0;
    return mjdMidnight + frac;
};

StarJs.Time.mjd2jct = function (mjd) {
    return (mjd - 51544.5) / 36525.0;
};

StarJs.Time.gmst = function (mjd) {
    /* TODO: move to global */
    var SECS = 86400; // 24*60*60 -- number of seconds in day;
    var mjd0 = Math.floor(mjd), ut = SECS * (mjd - mjd0), t0, t, gmst;
    t0 = StarJs.Time.mjd2jct(mjd0);
    t = StarJs.Time.mjd2jct(mjd);
    gmst = 24110.54841 + 8640184.812866 * t0 + 1.0027379093 * ut +
	(0.093104 - (6.2e-6) * t) * t * t;
    return StarJs.Math.PI2 / SECS * StarJs.Math.mod(gmst, SECS);
};
StarJs.Coord = {};

/** Precession matrix (in ecliptical coordinates) from epoch t1 to
 *  epoch t2.
 */
StarJs.Coord.precessionEclMatrix = function (t1, t2) {
    var dt = t2 - t1, p1, p2, pa;
    p1 = 174.876383889 * StarJs.Math.DEG2RAD +
	(((3289.4789 + 0.60622 * t1) * t1) +
	 ((-869.8089 - 0.50491 * t1) + 0.03536 * dt) * dt) * StarJs.Math.ARCS2RAD;
    p2 = ((47.0029 - (0.06603 - 0.000598 * t1) * t1) +
	  ((-0.03302 + 0.00598 * t1) + 0.000060 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    pa = ((5029.0966 + (2.22226 - 0.000042 * t1) * t1) +
	  ((1.11113 - 0.000042 * t1) - 0.000006 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    return StarJs.Vector.Matrix3.r_z(-(p1 + pa)).mult(StarJs.Vector.Matrix3.r_x(p2)).mult(StarJs.Vector.Matrix3.r_z(p1));
};

/** Precession matrix (in equatorial coordinates) from epoch t1 to
 *  epoch t2.
 */
StarJs.Coord.precessionEquMatrix = function (t1, t2) {
    var dt = t2 - t1, zeta, z, theta;
    zeta = ((2306.2181 + (1.39656 - 0.000139 * t1) * t1) +
	    ((0.30188 - 0.000344 * t1) + 0.017998 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    z = zeta + ((0.79280 + 0.000411 * t1) + 0.000205 * dt) * dt * dt * StarJs.Math.ARCS2RAD;
    theta  = ((2004.3109 - (0.85330 + 0.00217 * t1) * t1) -
	      ((0.42665 + 0.000217 * t1) + 0.041833 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    return StarJs.Vector.Matrix3.r_z(-z).mult(StarJs.Vector.Matrix3.r_y(theta)).mult(StarJs.Vector.Matrix3.r_z(-zeta));
    
};

/** Oliquity of the ecliptic.
    You may use StarJs.Solar.EPS as constant approximation.
 */
StarJs.Coord.eclipticObliquity = function (jct) {
    return (23.43929111 - (46.8150 + (0.00059 - 0.001813 * jct) * jct) * jct / 3600.0) * StarJs.Math.DEG2RAD;
};

/** Matrix for conversion from equatorial to ecliptic coordinate system.
 */
StarJs.Coord.equ2eclMatrix = function (jct) {
    var eps = StarJs.Coord.eclipticObliquity(jct);
    return StarJs.Vector.Matrix3.r_x(eps);
};

/** Matrix for conversion from ecliptic to equatorial coordinate system.
 */
StarJs.Coord.ecl2equMatrix = function (jct) {
    var eps = StarJs.Coord.eclipticObliquity(jct);
    return StarJs.Vector.Matrix3.r_x(-eps);
};

/** Convert from equatorial to horizontal coordinate system.
 */
StarJs.Coord.equ2hor = function (dec, tau, lat) {
    var hor_vec = StarJs.Vector.Matrix3.r_y(Math.PI / 2 - lat).apply(new StarJs.Vector.Polar3(tau, dec).toVector3());
    var hor_pol = new StarJs.Vector.Polar3(hor_vec);
    return {'h': hor_pol['theta'], 'az': hor_pol['phi']};
};

/** Convert from horizontal to equatorial coordinate system.
 */
StarJs.Coord.hor2equ = function (h, az, lat) {
    var equ_vec = StarJs.Vector.Matrix3.r_y(-Math.PI / 2 + lat).apply(new StarJs.Vector.Polar3(az, lat).toVector3());
    var equ_pol = new StarJs.Vector.Polar3(equ_vec);
    return {'dec': equ_pol['theta'], 'tau': equ_pol['phi']};
};
StarJs.Kepler = {};

/** Default number of iterations for iterative method. */
StarJs.Kepler['DEFAULT_ITERATIONS'] = 100;
/** Default iteration precision for iterative method. */
StarJs.Kepler['DEFAULT_PRECISION'] = 1e-9;

function _acceleratedModule(global, forein, heap) {
    "use asm";

    var PI = 3.1415926535897932385;
    var sin = global.Math.sin;
    var cos = global.Math.cos;
    var abs = global.Math.abs;
    var pow = global.Math.pow;
    var sqrt = global.Math.sqrt;
    var log = global.Math.log;
    var exp = global.Math.exp;

    var HEAPD = new global.Float64Array(heap);

    function mod(a, b) {
        a = +a;
        b = +b;
        return +(a % b); // TODO: check sign issues
    }
    
    function cosh (x) {
        x = +x;
        return (exp(x) + exp(-x))/2.0;
    };

    function eccAnomaly(m, ec, maxiter, prec) {
        m = +m;
        ec = +ec;
        maxiter = maxiter|0;
        prec = +prec;

        var f = 0.0;
        var sinm = 0.0;
        var e0 = 0.0;

        m = mod(m, 2.0*PI);

        /* Iterative solution of Kepler equation with Newthon's method.
         */
        //var e0 = 0; /* TODO: initial value depends on eccentricity */
        /* Gary R. Smith.  A simple, efficient starting value for the
         * iterative solution of Kepler's equation.  // Celestial
         * Mechanics and Dynamical Astronomy.  Volume 19, Number 2 (1979)
         * DOI: 10.1007/BF01796088.
         */
        sinm = sin(m);
        e0 = +(m + ec * sinm / (1.0 - sin(m+ec) + sinm));

        do {
            f = +(e0 - ec * sin(e0) - m);
            e0 = +(e0 - f / (1.0 - ec * cos(e0)));
            maxiter = (maxiter - 1)|0;
        } while ((~(~maxiter) > ~(~0)) & (f > prec));

        if (~(~maxiter) < ~(~0)) {
            e0 = 1.0/0.0;
        }
        return +e0;
    }

    function stumpff(e2, eps) {
        e2 = +e2;
        eps = +eps;

        var c1 = 0.0;
        var c2 = 0.0;
        var c3 = 0.0;
        var add = 1.0;
        var n = 1.0;

        do {
            c1 = c1 + add;
            add = add / (2.0*n);
            c2 = c2 + add;
            add = add / (2.0*n + 1.0);
            c3 = c3 + add;
            add = add * -e2;
            n = n + 1.0;
        } while (abs(add) >= eps);

        HEAPD[3] = c1;
        HEAPD[4] = c2;
        HEAPD[5] = c3;
    }

    function parabolic(gm, t0, t, q, ec, maxiter, prec) {
        gm = +gm;
        t0 = +t0;
        t = +t;
        q = +q;
        ec = +ec;
        maxiter = maxiter|0;
        prec = +prec;

        var e2 = 0.0;
        var e20 = 0.0;
        var fac = 0.0;
        var k = 0.0;
        var tau = 0.0;
        var a = 0.0;
        var b = 0.0;
        var u = 0.0;
        var u2 = 0.0;
        var r = 0.0;

        fac = 0.5 * ec;
        k = sqrt(gm/(q*(1.0+ec)));
        tau = sqrt(gm)*(t-t0);

        do {
            e20 = e2;
            a = 1.5 * sqrt(fac/(q*q*q))*tau;
            b = pow(sqrt(a*a + 1.0)+a, 1.0/3.0);
            u = b - 1.0/b;
            u2 = u*u;
            e2 = u2*(1.0-ec)/fac;
            stumpff(e2, prec);
            fac = 3.0 * ec * HEAPD[5];
        } while (abs(e2-e20) >= prec)

        HEAPD[0] = q*(1.0 - u2*HEAPD[4]/fac);
        HEAPD[1] = q*sqrt((1.0+ec)/fac)*u*HEAPD[3];
        HEAPD[2] = 0.0;
    }

    function hypAnom(mh, e, maxiter, prec) {
        mh = +mh;
        e = +e;
        maxiter = maxiter|0;
        prec = +prec;

        var i = 0;
        var h = 0.0;
        var f = 0.0;

        h = log(2.0*Math.abs(mh)/e+1.8);
        if (mh < 0.0) h = -h;

        do {
            f = e*StarJs.Math.sinh(h) - h - mh;
            h -= f/(e*cosh(h) - 1.0);
            i = i + 1;
            if (i === maxiter) {
                return 1.0/0.0;
            }
        } while (abs(f) > prec*(1.0 + abs(h+mh)));

        return h;
    };

    return {
        eccAnomaly: eccAnomaly,
        stumpff: stumpff,
        parabolic: parabolic,
        hypAnom: hypAnom
    };
}

var MEM = new ArrayBuffer(4096);
var HEAPD = new Float64Array(MEM);

_accelerated = _acceleratedModule(window, {}, MEM);

/** Compute eccentric anomaly for elliptical orbit.
 * @param m Mean anomaly.
 * @param ec Eccentricity (ec < 1)
 * @param maxiter Optional number of iterations.
 */
StarJs.Kepler.eccAnomaly = function (m, ec, maxiter) {
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    e0 = _accelerated.eccAnomaly(m, ec, maxiter, prec);

    return (isNaN(e0)) ? null : e0;
};

/** Compute position of a body on elliptical orbit.
 * @param gm Gravity constant.
 * @param m Mean anomaly.
 * @param a Semi-major axis.
 * @param ec Eccentricity.
 */
StarJs.Kepler.elliptic = function (gm, m, a, ec) {
    var k = Math.sqrt(gm / a), e = StarJs.Kepler.eccAnomaly(m, ec);
    var cosE = Math.cos(e), sinE = Math.sin(e);
    var fac = Math.sqrt((1.0 - ec) * (1.0 + ec));

    return new StarJs.Vector.Vector3(a * (cosE - ec), a * fac * sinE, 0.0);
    // var rho = 1.0 - ec * cosE;
    // var vel = StarJs.Vector.Vector3(-k * sinE / rho, k * fac * cosE / rho, 0.0);
};

/** Compute position of a body on parabolic orbit.
 @param gm Gravity constant.
 @param t0   time of pericenter
 @param t    time to calculate position for
 @param q    pericenter distance
 @param ec Eccentricity.
 @param maxiter Optional maximal number of iterations.
 */
StarJs.Kepler.parabolic = function (gm, t0, t, q, ec, maxiter) {
    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    _accelerated.parabolic(gm, t0, t, q, ec, maxiter, prec);    

    return new StarJs.Vector.Vector3(HEAPD[0], HEAPD[1], HEAPD[2]);
//     // res is position
//     r = q * (1+u2*c.c2*ec/fac);
//     var vel = new StarJs.Vector.Vector3(-k*res.y/r,
//                                         k*(res.x/r+ec),
//                                         0.0);
};

/** Compute hyperbolic anomaly. */
StarJs.Kepler.hypAnom = function (mh, e, maxiter) {
    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    h = _accelerated.hypAnom(mh, e, maxiter, prec);
    return h;
};

/** Compute position of a body on hyperbolic orbit.
 * @param gm Gravity constant
 * @param t0 Time of pericenter
 * @param t Time
 * @param a Semi-major axis
 * @param e Eccentricity.
 */
StarJs.Kepler.hyperbolic = function (gm, t0, t, a, e) {
    var k, mh, h, ch, sh, rho, fac;

    a = Math.abs(a);
    k = Math.sqrt(gm/a);
    mh = k*(t-t0)/a;
    h = _accelerated.hypAnom(mh, e,
                             StarJs.Kepler['DEFAULT_ITERATIONS'],
                             StarJs.Kepler['DEFAULT_PRECISION']);
    ch = StarJs.Math.cosh(h);
    sh = StarJs.Math.sinh(h);
    fac = Math.sqrt((e+1)*(e-1)); // (e*e-1) ?
    rho = e*ch - 1;
    return new StarJs.Vector.Vector3(a*(e-ch), a*fac*sh, 0);
//     // ret is position
//     vel = new StarJs.Vector.Vector3(-k*sh/rho, k*fac*ch/rho, 0);
};

/** Calculate Keplerian orbital position.
    @param gm   GM
    @param t0   time of pericenter
    @param t    time to calculate position for
    @param q    pericenter distance
    @param e    eccentricity
    @param pqr  Gauss' vector matrix
 */
StarJs.Kepler.keplerPos = function (gm, t0, t, q, e, pqr) {
    var M0 = 0.1, EPS = 0.1, m, delta, tau, invax, r;

    delta = Math.abs(1-e);
    invax = delta / q;
    tau = Math.sqrt(gm)*(t-t0);
    m = tau * Math.sqrt(invax*invax*invax);

    if ((m < M0) && (delta < EPS)) {
        r = StarJs.Kepler.parabolic(gm, t0, t, q, e);
    } else if (e < 1.0) {
        r = StarJs.Kepler.elliptic(gm, m, 1.0/invax, e);
    } else {
        r = StarJs.Kepler.hyperbolic(gm, t0, t, 1.0/invax, e);
    }

    return pqr.apply(r);
};

/** Return Gauss vectors matrix for given elements.
    @param omega longitude of the ascending node (radians)
    @param i     inclination (radians)
    @param w     argument of pericenter (radians)
 */
StarJs.Kepler.gaussVec = function (omega, i, w) {
    return StarJs.Vector.Matrix3.r_z(-omega).mult(StarJs.Vector.Matrix3.r_x(-i))
        .mult(StarJs.Vector.Matrix3.r_z(-w));
};
StarJs.Solar = {};

/**
 * Earth orbit's axial tilt.
 * @const
 * @type {number}
 */
StarJs.Solar.EPS = 23.43920111 * StarJs.Math.DEG2RAD;

/**
 * Return approximate position of Moon, quickly.
 * @param {number} mct Julian centuries from J2000. TODO
 * @return Object with fields ra and dec.
 */
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
    return {'ra': me['phi'], 'dec': me['theta']};

};

/**
 * Return approximate position of Moon, quickly.
 * @param {number} mct Julian centuries from J2000.
 * @return Object with fields ra and dec.
 */
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
    return {'ra': me['phi'], 'dec': me['theta']};
};

/** @const */
StarJs.Solar.EVENTS = [{
    'body': 'moon',
    'name': 'day',
    'title': "Moon",
    'sinh0': Math.sin((+8.0 / 60.0) * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxMoon
}, {
    'body': 'sun',
    'name': 'day',
    'title': "Sun",
    'sinh0': Math.sin((-50.0 / 60.0) * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
}, {
    'body': 'sun',
    'name': 'twilightC',
    'title': "Civil twilight",
    'sinh0': Math.sin(-6.0 * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
}, {
    'body': 'sun',
    'name': 'twilightN',
    'title': "Nautical twilight",
    'sinh0': Math.sin(-12.0 * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
}, {
    'body': 'sun',
    'name': 'twilightA',
    'title': "Astronomical twilight",
    'sinh0': Math.sin(-18.0 * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
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
        tau = StarJs.Time.gmst(mjd) + lambda - pos['ra'];
        return sphi * Math.sin(pos['dec']) +
            cphi * Math.cos(pos['dec']) * Math.cos(tau);
    }

    /** @const */
    var EVENTS = StarJs.Solar.EVENTS;

    var result = [], today, h, j, ptime, ctime, ntime, H = 1.0 / 24.0, posp = {}, pos0 = {}, posn = {}, cphi = Math.cos(phi), sphi = Math.sin(phi), name, evt, yp = {'sun': {}, 'moon': {}}, y0 = {'sun': {}, 'moon': {}}, yn = {'sun': {}, 'moon': {}}, interp = {};
    if (typeof tz === 'undefined') {
        tz = function (a) {
            return a;
        };
    }

    ptime = start;

    posp['moon'] = sinAlt(StarJs.Solar.approxMoon, ptime, lambda, cphi, sphi);
    posp['sun']  = sinAlt(StarJs.Solar.approxSun,  ptime, lambda, cphi, sphi);

    for (j = 0; j < EVENTS.length; j += 1) {
        evt = EVENTS[j];
        name = evt['name'];

        yp[evt['body']][name] = posp[evt['body']] - evt['sinh0'];
    }

    while (ptime < end) {
        today = {
            'moon': {
                'midnight': posp['moon']
            },
            'sun': {
                'midnight': posp['sun']
            }
        };
        for (h = 1; h < 24; h += 2) {
            ctime = ptime + H;
            ntime = ctime + H; // ntime = ctime + 2 * H;

            // Calc new positions...
            pos0['moon'] = sinAlt(StarJs.Solar.approxMoon, ctime, lambda,
                               cphi, sphi);
            pos0['sun'] = sinAlt(StarJs.Solar.approxSun, ctime, lambda,
                              cphi, sphi);
            posn['moon'] = sinAlt(StarJs.Solar.approxMoon, ntime, lambda,
                               cphi, sphi);
            posn['sun'] = sinAlt(StarJs.Solar.approxSun, ntime, lambda,
                              cphi, sphi);

            for (j = 0; j < EVENTS.length; j += 1) {
                evt = EVENTS[j];
                name = evt['name'];
                today[evt['body']][name] = today[evt['body']][name] || {};

                y0[evt['body']][name] = pos0[evt['body']] - evt['sinh0'];
                yn[evt['body']][name] = posn[evt['body']] - evt['sinh0'];

                interp = StarJs.Math.quadInterpolation(yp[evt['body']][name], y0[evt['body']][name], yn[evt['body']][name]);

                switch (interp['roots'].length) {
                case 0:
                    // No roots
                    break;
                case 1:
                    // Single root
                    if (yp[evt['body']][name] < 0.0) {
                        today[evt['body']][name]['rise'] = h + interp['roots'][0];
                    } else {
                        today[evt['body']][name]['set']  = h + interp['roots'][0];
                    }
                    break;
                case 2:
                    // Two roots
                    if (interp['ye'] < 0.0) {
                        today[evt['body']][name]['rise'] = h + interp['roots'][1];
                        today[evt['body']][name]['set']  = h + interp['roots'][0];
                    } else {
                        today[evt['body']][name]['rise'] = h + interp['roots'][0];
                        today[evt['body']][name]['set']  = h + interp['roots'][1];
                    }
                    break;
                }
            }

            // Next interval
            yp = yn;
            yn = {'moon': {}, 'sun': {}};

            ptime = ntime;
        }

        for (j = 0; j < EVENTS.length; j += 1) {
            evt = EVENTS[j];
            name = evt['name'];

            if (!today[evt['body']][name]['set'] && !today[evt['body']][name]['rise']) {
                today[evt['body']][name]['alwaysAbove'] = (pos0[evt['body']] > evt['sinh0']);
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

/** @constructor */
StarJs.Solar.Sun = function () {
    this['name'] = 'Sun';
};

/** @const */
StarJs.Solar.Sun.$POS = new StarJs.Vector.Vector3(0, 0, 0);

/** @const */
StarJs.Solar.GM = StarJs.Math.sqr(0.01720209895);

/**
 * Cooordinates of Sun are constant.
 */
StarJs.Solar.Sun.prototype['keplerCoord'] = function (t) {
    return StarJs.Solar.Sun.$POS;
};

/** @constructor */
StarJs.Solar.Body = function (name, params) {
    this['name'] = name;
    this['a']  = params['a'];
    this['ec'] = params['e'];
    this['m0'] = params['M0'];
    this['n']  = params['n'];
    this['omega'] = params['O'];
    this['i']  = params['i'];
    this['w']  = params['w'];
    this['t0'] = params['T0'];
    // ...
};

StarJs.Solar.Body.prototype['keplerCoord'] = function (t) {
    var P = 1.3970;
    var m = StarJs.Math.DEG2RAD * (this['m0'] + this['n'] * (t - this['t0']));
    var r = StarJs.Kepler.elliptic(StarJs.Solar.GM, m, this['a'], this['ec']);
    var pqr = StarJs.Kepler.gaussVec(
        StarJs.Math.DEG2RAD * (this['omega'] + P * t),
        StarJs.Math.DEG2RAD * this['i'],
        StarJs.Math.DEG2RAD * (this['w'] - this['omega']));
    return pqr.apply(r);
};

/** @const */
StarJs.Solar.BODIES = {
    'Sun':     new StarJs.Solar.Sun(),
    'Mercury': new StarJs.Solar.Body('Mercury', {
      'a' :  0.387099, 'e' : 0.205634, 'M0' : 174.7947, 'n'  : 149472.6738,
      'O' :  48.331,   'i' : 7.0048,   'w'  :  77.4552, 'T0' : 0.0
    }),
    'Venus':   new StarJs.Solar.Body('Venus', {
      'a' :  0.723332, 'e' : 0.006773, 'M0' :  50.4071, 'n'  : 58517.8149,
      'O' :  76.680,   'i' : 3.3946,   'w'  : 131.5718, 'T0' : 0.0
    }),
    'Earth':   new StarJs.Solar.Body('Earth', {
      'a' :  1.000000, 'e' : 0.016709, 'M0' : 357.5256, 'n'  : 35999.3720,
      'O' : 174.876,   'i' : 0.0000,   'w'  : 102.9400, 'T0' : 0.0
    }),
    'Mars':   new StarJs.Solar.Body('Mars', {
      'a' :  1.523692, 'e' : 0.093405, 'M0' :  19.3879, 'n'  : 19140.3023,
      'O' :  49.557,   'i' : 1.8496,   'w'  : 336.0590, 'T0' : 0.0
    }),
    'Jupiter':   new StarJs.Solar.Body('Jupiter', {
      'a' :  5.204267, 'e' : 0.048775, 'M0' :  18.8185, 'n'  : 3033.6272,
      'O' : 100.4908,  'i' : 1.3046,   'w'  :  15.5576, 'T0' : 0.0
    }),
    'Saturn':   new StarJs.Solar.Body('Saturn', {
      'a' :  9.582018, 'e' : 0.055723, 'M0' : 320.3477, 'n'  : 1213.8664,
      'O' : 113.6427,  'i' : 2.4852,   'w'  :  89.6567, 'T0' : 0.0
    }),
    'Uranus':   new StarJs.Solar.Body('Uranus', {
      'a' : 19.229412, 'e' : 0.044406, 'M0' : 142.9559, 'n'  : 426.9282,
      'O' :  73.9893,  'i' : 0.7726,   'w'  : 170.5310, 'T0' : 0.0
    }),
    'Neptune':   new StarJs.Solar.Body('Neptune', {
      'a' : 30.103658, 'e' : 0.011214, 'M0' : 267.7649, 'n'  : 217.9599,
      'O' : 131.7942,  'i' : 1.7680,   'w'  :  37.4435, 'T0' : 0.0
    })
};
