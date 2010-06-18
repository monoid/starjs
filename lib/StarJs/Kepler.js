/*global StarJs: false */
/*jslint onevar: false */
/**
 * @depends ../StarJs.js
 * @depends Math.js
 * @depends Vector.js
 */
StarJs.Kepler = (typeof StarJs.Kepler === 'undefined') ? {} : StarJs.Kepler;

/* @const */
StarJs.Kepler.DEFAULT_ITERATIONS = 100;
StarJs.Kepler.DEFAULT_PRECISION = 1e-9;

StarJs.Kepler.eccAnomaly = function (m, ec, maxiter) {
    var i, f, prec = StarJs.Kepler.DEFAULT_PRECISION;

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler.DEFAULT_ITERATIONS;
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
    } while (maxiter > 0 && (f < prec));
    return (maxiter > 0) ? e0 : null;
};

StarJs.Kepler.elliptic = function (gm, m, a, ec) {
    var k = Math.sqrt(gm / a), e = StarJs.Kepler.eccAnomaly(m, e);
    var cosE = Math.cos(e), sinE = Math.sin(e);
    var fac = Math.sqrt((1.0 - ec) * (1.0 + ec));
    return new StarJs.Vector.Vector3(a * (cosE - ec), a * fac * sinE, 0.0);
    // var rho = 1.0 - ec * cosE;
    // var vel = StarJs.Vector.Vector3(-k * sinE / rho, k * fac * cosE / rho, 0.0);
};

StarJs.Kepler.parabolic = function (gm, t0, t, q, ec, maxiter) {
    function stumpff(e2, ret) {
        var eps = StarJs.Kepler.DEFAULT_PRECISION;
        var n, add, c1, c2, c3;

        c1 = c2 = c3 = 0.0; 
        add = n = 1.0;

        do {
            c1 += add;  add /= (2.0*n);
            c2 += add;  add /= (2.0*n+1.0);
            c3 += add;  add *= -E2; 
            n += 1.0;
        }
        while (fabs(add) >= eps);

        ret.c1 = c1;
        ret.c2 = c2;
        ret.c3 = c3;
    }

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    var i = 0, e2 = 0, e20, fac, c = {}, k, tau, a, b, u, u2, r;
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    fac = 0.5 * ec;
    k   = Math.sqrt(gm/(q*(1+ec)));
    tau = Math.sqrt(gm)*(t-t0);
    do {
        ++i;
        e20 = e2;
        a = 1.5 * Math.sqrt(fac/(q*q*q))*tau;
        b = Math.pow(Math.sqrt(a*a + 1)+a, 1/3.0);
        u = b - 1.0/b;
        u2 = u*u;
        e2 = u2*(1-ec)/fac;
        stumpff(e2, c);
        fac = 3.0 * ec * c.c3;
    } while (Math.abs(e2-e20) >= prec);
    r = q * (1+u2*c.c2*ec/fac);
    return new StarJs.Vector.Vector3(q*(1-u2*c.c2/fac),
                                     q*Math.sqrt((1+e)/fac)*u*c.c1,
                                     0);
//     // res is position
//     var vel = new StarJs.Vector.Vector3(-k*res.y/r,
//                                         k*(res.x/r+ec),
//                                         0.0);
};

StarJs.Kepler.gaussVec = function (omega, i, w) {
    return StarJs.Math.Matrix3.r_z(-omega).mult(StarJs.Math.Matrix3.r_x(-i))
        .mult(StarJs.Math.Matrix3.r_z(-w));
};
