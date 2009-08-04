StarJs = (typeof StarJs == 'undefined') ? {} : StarJs;
StarJs.Vector = (typeof StarJs.Vector == 'undefined') ? {} : StarJs.Vector;

StarJs.Vector.Vector3 = function (x, y, z) {
    this._x = x;
    this._y = y;
    this._z = z;
};

with (StarJs.Vector.Vector3) {
    prototype.len = function () {
        return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    };

    prototype.add = function (v3) {
        return new StarJs.Vector.Vector3(this.x + v3.x,
                                         this.y + v3.y,
                                         this.z + v3.z);
    };

    prototype.sub = function (v3) {
        return new StarJs.Vector.Vector3(this.x - v3.x,
                                         this.y - v3.y,
                                         this.z - v3.z);
    };

    prototype.neg = function () {
        return new StarJs.Vector.Vector3(-this.x, -this.y, -this.z);
    };

    prototype.scale = function (a) {
        return new StarJs.Vector.Vector3(a*this.x, a*this.y, a*this.z);
    };

    prototype.clone = function () {
        return new StarJs.Vector.Vector3(this.x, this.y, this.z);
    };
};

StarJs.Vector.Polar3 = function (az_or_v3, elev, rad) {
    var alen = arguments.length;
    if (alen == 2) rad = 1.0;
    if (alen == 2 || alen == 3) {
        this.phi = az_or_v3;
        this.theta = elev;
        this.rad = rad;
    } else {
        var v3 = az_or_v3;
        var rho2 = v3.x*v3.x + v3.y*v3.y;
        this.rad  = Math.sqrt(rho2 + v3.z*v3.z);
        this.phi = (v3.x == 0.0 && v3.y == 0.0) ?  0.0 : Math.atan2(v3.x, v3.y);
        if (this.phi<0.0) this.phi += StarJs.Math.PI2;
        var rho = Math.sqrt(rho2);
        this.theta = (m.z == 0.0 && rho == 0.0) ? 0.0 : Math.atan2(v3.z, rho);
    }
};

with (StarJs.Vector.Polar3) {
    prototype.toVector3 = function () {
        var ct = Math.cos(this.theta);
        var x = this.rad*ct*Math.cos(this.phi);
        var y = this.rad*ct*Math.sin(this.phi);
        var z = this.rad*Math.sin(this.theta);
        return new StarJs.Vector.Vector3(x, y, z);
    };

    prototype.clone = function () {
        return new StarJs.Vector.Polar3(this.rad, this.phi, this.theta);
    };
};

StarJs.Vector.Matrix3 = function (v1, v2, v3) {
    if (arguments.length == 3) {
	this.mat = [[v1.x,v1.y,v1.z], [v2.x,v2.y,v2.z], [v3.x,v3.y,v3.z]];
    } else {
	if (v2) {
	    this.mat = v1;
	} else {
	    this.mat = [[v1[0][0], v1[0][1], v1[0][2]],
			[v1[1][0], v1[1][1], v1[1][2]],
			[v1[2][0], v1[2][1], v1[2][2]]];
	}
    }
};

with (StarJs.Vector.Matrix3) {
    prototype.apply = function (v) {
        var l = this.mat[0];
        var x = l[0]*v.x + l[1]*v.y + l[2]*v.z;
        l = this.mat[1];
        var y = l[0]*v.x + l[1]*v.y + l[2]*v.z;
        l = this.mat[1];
        var z = l[0]*v.x + l[1]*v.y + l[2]*v.z;

        return new StarJs.Vector.Vector3(x, y, z);
    };

    prototype.mult = function (matrix) {

    };
};

StarJs.Vector.Matrix3.R_x = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
};

StarJs.Vector.Matrix3.R_y = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    
};

StarJs.Vector.Matrix3.R_z = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    
};
