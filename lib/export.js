window['StarJs'] = StarJs;

/* Math */
StarJs['Math'] = StarJs.Math;
StarJs.Math['DEG2RAD'] = StarJs.Math.DEG2RAD;
StarJs.Math['RAD2DEG'] = StarJs.Math.RAD2DEG;
StarJs.Math['RAD2ARCS'] = StarJs.Math.RAD2ARCS;
StarJs.Math['ARCS2RAD'] = StarJs.Math.ARCS2RAD;
StarJs.Math['PI2'] = StarJs.Math.PI2;
StarJs.Math['sqr'] = StarJs.Math.sqr;
StarJs.Math['frac'] = StarJs.Math.frac;
StarJs.Math['mod'] = StarJs.Math.mod;
StarJs.Math['Dms'] = StarJs.Math.Dms;
StarJs.Math['dms2deg'] = StarJs.Math.dms2deg;
StarJs.Math.Dms.prototype['dms2deg'] = StarJs.Math.Dms.prototype.dms2deg;
StarJs.Math['deg2dms'] = StarJs.Math.deg2dms;
StarJs.Math['angle2hms'] = StarJs.Math.angle2hms;
StarJs.Math['quadInterpolation'] = StarJs.Math.quadInterpolation;
StarJs.Math['sinh'] = StarJs.Math.sinh;
StarJs.Math['cosh'] = StarJs.Math.cosh;

/* Coord */
StarJs['Coord'] = StarJs.Coord;
StarJs.Coord['precessionEclMatrix'] = StarJs.Coord.precessionEclMatrix;
StarJs.Coord['precessionEquMatrix'] = StarJs.Coord.precessionEquMatrix;
StarJs.Coord['eclipticObliquity'] = StarJs.Coord.eclipticObliquity;
StarJs.Coord['equ2eclMatrix'] = StarJs.Coord.equ2eclMatrix;
StarJs.Coord['ecl2equMatrix'] = StarJs.Coord.ecl2equMatrix;
StarJs.Coord['equ2hor'] = StarJs.Coord.equ2hor;
StarJs.Coord['hor2equ'] = StarJs.Coord.hor2equ;

/* Kepler */
StarJs['Kepler'] = StarJs.Kepler;
StarJs.Kepler['eccAnomaly'] = StarJs.Kepler.eccAnomaly;
StarJs.Kepler['elliptic'] = StarJs.Kepler.elliptic;
StarJs.Kepler['parabolic'] = StarJs.Kepler.parabolic;
StarJs.Kepler['hypAnom'] = StarJs.Kepler.hypAnom;
StarJs.Kepler['hyperbolic'] = StarJs.Kepler.hyperbolic;
StarJs.Kepler['keplerPos'] = StarJs.Kepler.keplerPos;
StarJs.Kepler['gaussVec'] = StarJs.Kepler.gaussVec;

/* Solar */
StarJs['Solar'] = StarJs.Solar;
StarJs.Solar['EPS'] = StarJs.Solar.EPS;
StarJs.Solar['approxMoon'] = StarJs.Solar.approxMoon;
StarJs.Solar['approxSun'] = StarJs.Solar.approxSun;
StarJs.Solar['sunAndMoonEvents'] = StarJs.Solar.sunAndMoonEvents;
StarJs.Solar['Sun'] = StarJs.Solar.Sun;
StarJs.Solar['GM'] = StarJs.Solar.GM;
StarJs.Solar['Body'] = StarJs.Solar.Body;
StarJs.Solar['BODIES'] = StarJs.Solar.BODIES;

/* Time */
StarJs['Time'] = StarJs.Time;
StarJs.Time['time2mjd'] = StarJs.Time.time2mjd;
StarJs.Time['mjd2dt'] = StarJs.Time.mjd2dt;
StarJs.Time['hms2hour'] = StarJs.Time.hms2hour;
StarJs.Time['hour2hms'] = StarJs.Time.hour2hms;
StarJs.Time['dt2mjd'] = StarJs.Time.dt2mjd;
StarJs.Time['mjd2jct'] = StarJs.Time.mjd2jct;
StarJs.Time['gmst'] = StarJs.Time.gmst; 

/* Vector */
StarJs['Vector'] = StarJs.Vector;
StarJs.Vector['Vector3'] = StarJs.Vector.Vector3;

StarJs.Vector.Vector3.prototype['len'] = StarJs.Vector.Vector3.prototype.len;
StarJs.Vector.Vector3.prototype['add'] = StarJs.Vector.Vector3.prototype.add;
StarJs.Vector.Vector3.prototype['sub'] = StarJs.Vector.Vector3.prototype.sub;
StarJs.Vector.Vector3.prototype['neg'] = StarJs.Vector.Vector3.prototype.neg;
StarJs.Vector.Vector3.prototype['scale'] = StarJs.Vector.Vector3.prototype.scale;
StarJs.Vector.Vector3.prototype['clone'] = StarJs.Vector.Vector3.prototype.clone;

StarJs.Vector['Polar3'] = StarJs.Vector.Polar3;

StarJs.Vector.Polar3.prototype['toVector3'] = StarJs.Vector.Polar3.prototype.toVector3;
StarJs.Vector.Polar3.prototype['clone'] = StarJs.Vector.Polar3.prototype.clone;

StarJs.Vector['Matrix3'] = StarJs.Vector.Matrix3;

StarJs.Vector.Matrix3.prototype['apply'] = StarJs.Vector.Matrix3.prototype.apply;
StarJs.Vector.Matrix3.prototype['mult'] = StarJs.Vector.Matrix3.prototype.mult;

StarJs.Vector.Matrix3['r_x'] = StarJs.Vector.Matrix3.r_x;
StarJs.Vector.Matrix3['r_y'] = StarJs.Vector.Matrix3.r_y;
StarJs.Vector.Matrix3['r_y'] = StarJs.Vector.Matrix3.r_y;
