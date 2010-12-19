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

