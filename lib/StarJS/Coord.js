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
    return StarJs.Vector.r_z(-(p1+pa)).mult(StarJs.Vector.r_x(p2)).mult(StarJs.Vector.r_z(p1));
};

StarJs.Coord.precessionEquMatrix = function (t1, t2) {
    var dt = t2-t1, zeta, z, theta;
    zeta = ((2306.2181+(1.39656-0.000139*t1)*t1)+
	    ((0.30188-0.000344*t1)+0.017998*dt)*dt)*dt*StarJs.Math.ARCS2RAD;
    z = zeta + ((0.79280+0.000411*t1)+0.000205*dt)*dt*dt*StarJs.Math.ARCS2RAD;
    theta  = ((2004.3109-(0.85330+0.00217*t1)*t1)-
	      ((0.42665+0.000217*t1)+0.041833*dt)*dt)*dt*StarJs.Math.ARCS2RAD;
    return StarJs.Vector.r_z(-z).mult(StarJs.Vector.r_y(theta)).mult(StarJs.Vector.r_z(-zeta));
    
};
