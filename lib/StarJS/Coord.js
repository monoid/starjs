StarJs = (typeof StarJs == 'undefined') ? {} : StarJs;
StarJs.Coord = (typeof StarJs.Coord == 'undefined') ? {} : StarJs.Coord;

StarJs.Coord.precessionEclMatrix = function (t1, t2) {
    var dt = t2-t1, p1, p2, pa;
    p1 = 174.876383889*StarJs.Math.DEG2RAD +
	(((3289.4789+0.60622*t1)*t1) +
	 ((-869.8089-0.50491*t1) + 0.03536*dt)*dt)*StarJs.Math.ARCS2RAD;
    p2 = ((47.0029-(0.06603-0.000598*t1)*t1)+
	  ((-0.03302+0.00598*t1)+0.000060*dt)*dt)*dt*StarJs.Math.ARCS2RAD;
    pa = ((5029.0966+(2.22226-0.000042*t1)*t1)+
	  ((1.11113-0.000042*t1)-0.000006*dt)*dt)*dt*StarJs.Math.ARCS2RAD;
    return StarJs.Vector.Matrix3.r_z(-(p1+pa)).mult(StarJs.Vector.Matrix3.r_x(p2)).mult(StarJs.Vector.Matrix3.r_z(p1));
};

StarJs.Coord.precessionEquMatrix = function (t1, t2) {
    var dt = t2-t1, zeta, z, theta;
    zeta = ((2306.2181+(1.39656-0.000139*t1)*t1)+
	    ((0.30188-0.000344*t1)+0.017998*dt)*dt)*dt*StarJs.Math.ARCS2RAD;
    z = zeta + ((0.79280+0.000411*t1)+0.000205*dt)*dt*dt*StarJs.Math.ARCS2RAD;
    theta  = ((2004.3109-(0.85330+0.00217*t1)*t1)-
	      ((0.42665+0.000217*t1)+0.041833*dt)*dt)*dt*StarJs.Math.ARCS2RAD;
    return StarJs.Vector.Matrix3.r_z(-z).mult(StarJs.Vector.Matrix3.r_y(theta)).mult(StarJs.Vector.Matrix3.r_z(-zeta));
    
};

StarJs.Coord.equ2eclMatrix = function (jct) {
    var eps = (23.43929111 - (46.8150 + (0.00059-0.001813*jct)*jct)*jct/3600.0)*StarJs.DEG2RAD;
    return StarJs.Vector.Matrix3.r_x(eps);
};

StarJs.Coord.ecl2equMatrix = function (jtc) {
    var eps = (23.43929111 - (46.8150 + (0.00059-0.001813*jct))*jct)*jct*StarJs.ARCS2RAD;
    return StarJs.Vector.Matrix3.r_x(-eps);
};

StarJs.Coord.equ2hor = function (dec, tau, lat) {
    var hor_vec = StarJs.Vector.Matrix3.r_y(Math.PI/2-lat).apply(new StarJs.Vector.Polar3(tau, dec).toVector3());
    var hor_pol = new StarJs.Vector.Polar3(hor_vec);
    return {h: hor_pol.theta, az: hor_pol.phi};
};

StarJs.Coord.hor2equ = function (h, az, lat) {
    var equ_vec = StarJs.Vector.Matrix3.r_y(-Math.PI/2+lat).apply(new StarJs.Vector.Polar3(tau, dec).toVector3());
    var equ_pol = new StarJs.Vector.Polar3(hor_vec);
    return {dec: equ_pol.theta, tau: equ_pol.phi};
};

