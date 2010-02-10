StarJs.Vector = {};
StarJs['Vector'] = StarJs.Vector;

/** @constructor */
StarJs.Vector.Vector3 = function (x, y, z) {
    this['x'] = x;
    this['y'] = y;
    this['z'] = z;
};
StarJs.Vector['Vector3'] = StarJs.Vector.Vector3;

(function (p) {
    p.len = function () {
        return Math.sqrt(this['x'] * this['x'] + this['y'] * this['y'] + this['z'] * this['z']);
    };
    p['len'] = p.len;

    p.add = function (v3) {
        return new StarJs.Vector.Vector3(this['x'] + v3['x'],
                                         this['y'] + v3['y'],
                                         this['z'] + v3['z']);
    };
    p['add'] = p.add;

    p.sub = function (v3) {
        return new StarJs.Vector.Vector3(this['x'] - v3['x'],
                                         this['y'] - v3['y'],
                                         this['z'] - v3['z']);
    };
    p['sub'] = p.sub;

    p.neg = function () {
        return new StarJs.Vector.Vector3(-this['x'], -this['y'], -this['z']);
    };
    p['neg'] = p.neg;

    p.scale = function (a) {
        return new StarJs.Vector.Vector3(a * this['x'], a * this['y'], a * this['z']);
    };
    p['scale'] = p.scale;

    p.clone = function () {
        return new StarJs.Vector.Vector3(this['x'], this['y'], this['z']);
    };
    p['clone'] = p.clone;
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
StarJs.Vector['Polar3'] = StarJs.Vector.Polar3;

(function (p) {
    p.toVector3 = function () {
        var ct = Math.cos(this['theta']);
        var rad = this['rad']
        var x = rad * ct * Math.cos(this['phi']);
        var y = rad * ct * Math.sin(this['phi']);
        var z = rad * Math.sin(this['theta']);
        return new StarJs.Vector.Vector3(x, y, z);
    };
    p['toVector'] = p.toVector;

    p.clone = function () {
        return new StarJs.Vector.Polar3(this['rad'], this['phi'], this['theta']);
    };
    p['clone'] = p.clone;
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
StarJs.Vector['Matrix3'] = StarJs.Vector.Matrix3;

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
    p['apply'] = p.apply;

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
    p['mult'] = p.mult;
}(StarJs.Vector.Matrix3.prototype));

StarJs.Vector.Matrix3.r_x = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[1.0, 0.0, 0.0],
                                      [0.0,  cp,  sp],
                                      [0.0, -sp,  cp]]);
};
StarJs.Vector.Matrix3['r_x'] = StarJs.Vector.Matrix3.r_x;

StarJs.Vector.Matrix3.r_y = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[ cp, 0.0, -sp],
                                      [0.0, 1.0, 0.0],
                                      [ sp, 0.0,  cp]]);
};
StarJs.Vector.Matrix3['r_y'] = StarJs.Vector.Matrix3.r_y;

StarJs.Vector.Matrix3.r_z = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[ cp,  sp, 0.0],
                                      [-sp,  cp, 0.0],
                                      [0.0, 0.0, 1.0]]);
    
};
StarJs.Vector.Matrix3['r_y'] = StarJs.Vector.Matrix3.r_y;
