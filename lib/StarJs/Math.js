StarJs.Math = {};
StarJs['Math'] = StarJs.Math;

/** @const */
StarJs.Math.DEG2RAD = Math.PI / 180.0;
StarJs.Math['DEG2RAD'] = Math.PI / 180.0;
/** @const */
StarJs.Math.RAD2DEG = 180.0 / Math.PI;
StarJs.Math['RAD2DEG'] = 180.0 / Math.PI;
/** @const */
StarJs.Math.RAD2ARCS = 648000.0 / Math.PI;
StarJs.Math['RAD2ARCS'] = 648000.0 / Math.PI;
/** @const */
StarJs.Math.ARCS2RAD = Math.PI / 648000.0;
StarJs.Math['ARCS2RAD'] = Math.PI / 648000.0;
/** @const */
StarJs.Math.PI2 = 2.0 * Math.PI;
StarJs.Math['PI2'] = 2.0 * Math.PI;

/** @const */
StarJs.Math.ANGLE2HMS = 12.0 / Math.PI;

StarJs.Math.sqr = function (x) {
    return x * x;
};
StarJs.Math['sqr'] = StarJs.Math.sqr;

StarJs.Math.frac = function (x) {
    return x - Math.floor(x);
};
StarJs.Math['frac'] = StarJs.Math.frac;

StarJs.Math.mod = function (x, r) {
    return r * StarJs.Math.frac(x / r);
};
StarJs.Math['mod'] = StarJs.Math.mod;

/** @constructor */
StarJs.Math.Dms = function (sign_or_angle, degree, minute, second) {
    if (arguments.length === 4) {
        this.sign = sign_or_angle;
        this.degree = degree;
        this.minute = minute;
        this.second = second;
    } else {
        var angle = sign_or_angle;
        this.sign = 1;
        if (angle < 0) {
            angle = -angle;
            this.sign = -1;
        }
        this.degree = Math.floor(angle);
        angle = 60 * (angle - this.degree);
        this.minute = Math.floor(angle);
        angle = 60 * (angle - this.minute);
        this.second = angle;
    }
};
StarJs.Math['Dms'] = StarJs.Math.Dms;

StarJs.Math.dms2deg = function (sign, deg, minute, second) {
    return sign * (deg + minute / 60.0 + second / 3600.0);
};
StarJs.Math['dms2deg'] = StarJs.Math.dms2deg;

StarJs.Math.Dms.prototype.dms2deg = function () {
    return StarJs.Math.dms2deg(this.sign, this.degree, this.minute, this.second);
};
StarJs.Math.Dms.prototype['dms2deg'] = StarJs.Math.Dms.prototype.dms2deg;

StarJs.Math.Dms.deg2dms = function (deg) {
    return new StarJs.Math.Dms(deg);
};
StarJs.Math.Dms['deg2dms'] = StarJs.Math.Dms.deg2dms;

StarJs.Math.angle2hms = function (angle) {
    angle = StarJs.Math.mod(angle, StarJs.Math.PI2);
    var a = StarJs.Math.ANGLE2HMS * angle, res = StarJs.Math.deg2dms(a);
    // Change field names and remove sign field as it is always 1.
    res.hour = res.deg;
    delete res.deg;
    delete res.sign;
    
    return res;
};
StarJs.Math['angle2hms'] = StarJs.Math.angle2hms;

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
    return {xe: xe, ye: ye, roots: roots};
};
StarJs.Math['quadInterpolation'] = StarJs.Math.quadInterpolation;
