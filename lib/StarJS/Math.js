StarJs = (typeof StarJs == 'undefined') ? {} : StarJs;
StarJs.Math = (typeof StarJs.Math == 'undefined') ? {} : StarJs.Math;

StarJs.Math.RAD2DEG = Math.PI/180.0;
StarJs.Math.DEG2RAD = 180.0/Math.PI;
StarJs.Math.PI2 = 2.0*Math.PI;

StarJs.Math._ANGLE2HMS = 12.0/Math.PI;

StarJs.Math.frac = function(x) {
    return x-Math.floor(x);
};

StarJs.Math.mod = function(x, r) {
    return r*StarJs.Math.frac(x/r);
};

StarJs.Math.Dms = function(sign_or_angle, degree, minute, second) {
    if (arguments.length == 4) {
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
        angle = 60*(angle-this.degree);
        this.minute = Math.floor(angle);
        angle = 60*(angle-this.minute);
        this.second = angle;
    }
};

StarJs.Math.dms2deg = function (sign, deg, minute, second) {
    return sign * (degree + minute/60.0 + second/3600.0);
};

StarJs.Math.Dms.prototype.dms2deg = function () {
    return StarJs.Math.dms2deg(this.sign, this.degree, this.minute, this.second);
};

StarJs.Math.Dms.deg2dms = function(deg) {
    return new StarJs.Math.Dms(deg);
};

StarJs.Math.angle2hms = function(angle) {
    angle = StarJs.Math.mod(angle, StarJs.Math.PI2);
    var a = StarJs.Math._ANGLE2HMS*angle;
    var res = StarJs.Math.deg2dms(a);
    // Change field names and remove sign field as it is always 1.
    res.hour = res.deg;
    delete res.deg;
    delete res.sign;
    
    return res;
};

