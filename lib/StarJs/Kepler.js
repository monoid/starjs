StarJs.Kepler = {};
StarJs['Kepler'] = StarJs.Kepler;

StarJs.Kepler['DEFAULT_ITERATIONS'] = 100;
StarJs.Kepler['DEFAULT_PRECISION'] = 1e-9;

StarJs.Kepler.eccAnomaly = function (m, ec, maxiter) {
    var i, f, prec = StarJs.Kepler['DEFAULT_PRECISION'];

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    m = StarJs.Math.mod(m, StarJs.Math.PI2);

    /* Iterative solution of Kepler equation with Newthon's method.
     */
    //var e0 = 0; /* TODO: initial value depends on eccentricity */
    /* Gary R. Smith.  A simple, efficient starting value for the
     * iterative solution of Kepler's equation.  // Celestial
     * Mechanics and Dynamical Astronomy.  Volume 19, Number 2 (1979)
     * DOI: 10.1007/BF01796088.
     */
    var sinm = Math.sin(m);
    var e0 = m + ec * sinm / (1 - Math.sin(m + ec) + sinm);

    do {
        f = e0 - ec * Math.sin(e0) - m;
        e0 -= f / (1.0 - ec * Math.cos(e0));
    } while (maxiter-- > 0 && (f > prec));
    return (maxiter > 0) ? e0 : null;
};
StarJs.Kepler['eccAnomaly'] = StarJs.Kepler.eccAnomaly;

StarJs.Kepler.elliptic = function (gm, m, a, ec) {
    var k = Math.sqrt(gm / a), e = StarJs.Kepler.eccAnomaly(m, ec);
    var cosE = Math.cos(e), sinE = Math.sin(e);
    var fac = Math.sqrt((1.0 - ec) * (1.0 + ec));
    return new StarJs.Vector.Vector3(a * (cosE - ec), a * fac * sinE, 0.0);
    // var rho = 1.0 - ec * cosE;
    // var vel = StarJs.Vector.Vector3(-k * sinE / rho, k * fac * cosE / rho, 0.0);
};
StarJs.Kepler['elliptic'] = StarJs.Kepler.elliptic;

StarJs.Kepler.gaussVec = function (omega, i, w) {
    return StarJs.Vector.Matrix3.r_z(-omega).mult(StarJs.Vector.Matrix3.r_x(-i))
        .mult(StarJs.Vector.Matrix3.r_z(-w));
};
StarJs.Kepler['gaussVec'] = StarJs.Kepler.gaussVec;
