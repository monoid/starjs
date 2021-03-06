StarJs.Kepler = {};

/** Default number of iterations for iterative method. */
StarJs.Kepler['DEFAULT_ITERATIONS'] = 100;
/** Default iteration precision for iterative method. */
StarJs.Kepler['DEFAULT_PRECISION'] = 1e-9;

/** Compute eccentric anomaly for elliptical orbit.
 * @param {number} m Mean anomaly.
 * @param {number} ec Eccentricity (ec < 1)
 * @param {number=} maxiter Optional number of iterations.
 */
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

/** Compute position of a body on elliptical orbit.
 * @param {number} gm Gravity constant.
 * @param {number} m Mean anomaly.
 * @param {number} a Semi-major axis.
 * @param {number} ec Eccentricity.
 */
StarJs.Kepler.elliptic = function (gm, m, a, ec) {
    var k = Math.sqrt(gm / a), e = StarJs.Kepler.eccAnomaly(m, ec);
    var cosE = Math.cos(e), sinE = Math.sin(e);
    var fac = Math.sqrt((1.0 - ec) * (1.0 + ec));
    return new StarJs.Vector.Vector3(a * (cosE - ec), a * fac * sinE, 0.0);
    // var rho = 1.0 - ec * cosE;
    // var vel = StarJs.Vector.Vector3(-k * sinE / rho, k * fac * cosE / rho, 0.0);
};

/** Compute position of a body on parabolic orbit.
 @param {number} gm Gravity constant.
 @param {number} t0   time of pericenter
 @param {number} t    time to calculate position for
 @param {number} q    pericenter distance
 @param {number} ec Eccentricity.
 @param {number=} maxiter Optional maximal number of iterations.
 */
StarJs.Kepler.parabolic = function (gm, t0, t, q, ec, maxiter) {
    function stumpff(e2, ret) {
        var eps = StarJs.Kepler['DEFAULT_PRECISION'];
        var n, add, c1, c2, c3;

        c1 = c2 = c3 = 0.0; 
        add = n = 1.0;

        do {
            c1 += add;  add /= (2.0*n);
            c2 += add;  add /= (2.0*n+1.0);
            c3 += add;  add *= -e2; 
            n += 1.0;
        }
        while (Math.abs(add) >= eps);

        ret.c1 = c1;
        ret.c2 = c2;
        ret.c3 = c3;
    }

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    var e2 = 0, e20, fac, c = {}, k, tau, a, b, u, u2, r;
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    fac = 0.5 * ec;
    k   = Math.sqrt(gm/(q*(1+ec)));
    tau = Math.sqrt(gm)*(t-t0);
    do {
        e20 = e2;
        a = 1.5 * Math.sqrt(fac/(q*q*q))*tau;
        b = Math.pow(Math.sqrt(a*a + 1)+a, 1/3.0);
        u = b - 1.0/b;
        u2 = u*u;
        e2 = u2*(1-ec)/fac;
        stumpff(e2, c);
        fac = 3.0 * ec * c.c3;
    } while (Math.abs(e2-e20) >= prec);
    return new StarJs.Vector.Vector3(q*(1-u2*c.c2/fac),
                                     q*Math.sqrt((1+ec)/fac)*u*c.c1,
                                     0);
//     // res is position
//     r = q * (1+u2*c.c2*ec/fac);
//     var vel = new StarJs.Vector.Vector3(-k*res.y/r,
//                                         k*(res.x/r+ec),
//                                         0.0);
};

/** Compute hyperbolic anomaly.
 * @param {number} mh
 * @param {number} e
 * @param {number=} maxiter
 */
StarJs.Kepler.hypAnom = function (mh, e, maxiter) {
    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    var i = 0, h, f;

    h = Math.log(2.0*Math.abs(mh)/e+1.8);
    if (mh < 0.0) h = -h;

    do {
        f = e*StarJs.Math.sinh(h) - h - mh;
        h -= f/(e*StarJs.Math.cosh(h) - 1.0);
        ++i;
        if (i === maxiter) {
            // TODO: throw exception?
            return null;
        }
    } while (Math.abs(f) > prec*(1.0 + Math.abs(h+mh)));

    return h;
};

/** Compute position of a body on hyperbolic orbit.
 * @param {number} gm Gravity constant
 * @param {number} t0 Time of pericenter
 * @param {number} t Time
 * @param {number} a Semi-major axis
 * @param {number} e Eccentricity.
 */
StarJs.Kepler.hyperbolic = function (gm, t0, t, a, e) {
    var k, mh, h, ch, sh, rho, fac;

    a = Math.abs(a);
    k = Math.sqrt(gm/a);
    mh = k*(t-t0)/a;
    h = StarJs.Kepler.hypAnom(mh, e);
    ch = StarJs.Math.cosh(h);
    sh = StarJs.Math.sinh(h);
    fac = Math.sqrt((e+1)*(e-1)); // (e*e-1) ?
    rho = e*ch - 1;
    return new StarJs.Vector.Vector3(a*(e-ch), a*fac*sh, 0);
//     // ret is position
//     vel = new StarJs.Vector.Vector3(-k*sh/rho, k*fac*ch/rho, 0);
};

/** Calculate Keplerian orbital position.
    @param {number} gm   GM
    @param {number} t0   time of pericenter
    @param {number} t    time to calculate position for
    @param {number} q    pericenter distance
    @param {number} e    eccentricity
    @param {number} pqr  Gauss' vector matrix
 */
StarJs.Kepler.keplerPos = function (gm, t0, t, q, e, pqr) {
    var M0 = 0.1, EPS = 0.1, m, delta, tau, invax, r;

    delta = Math.abs(1-e);
    invax = delta / q;
    tau = Math.sqrt(gm)*(t-t0);
    m = tau * Math.sqrt(invax*invax*invax);

    if ((m < M0) && (delta < EPS)) {
        r = StarJs.Kepler.parabolic(gm, t0, t, q, e);
    } else if (e < 1.0) {
        r = StarJs.Kepler.elliptic(gm, m, 1.0/invax, e);
    } else {
        r = StarJs.Kepler.hyperbolic(gm, t0, t, 1.0/invax, e);
    }

    return pqr.apply(r);
};

/** Return Gauss vectors matrix for given elements.
    @param {number} omega longitude of the ascending node (radians)
    @param {number} i     inclination (radians)
    @param {number} w     argument of pericenter (radians)
 */
StarJs.Kepler.gaussVec = function (omega, i, w) {
    return StarJs.Vector.Matrix3.r_z(-omega).mult(StarJs.Vector.Matrix3.r_x(-i))
        .mult(StarJs.Vector.Matrix3.r_z(-w));
};
