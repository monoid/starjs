StarJs.Kepler = {};

/** Default number of iterations for iterative method. */
StarJs.Kepler['DEFAULT_ITERATIONS'] = 100;
/** Default iteration precision for iterative method. */
StarJs.Kepler['DEFAULT_PRECISION'] = 1e-9;

function _acceleratedModule(global, forein, heap) {
    "use asm";

    var PI = 3.1415926535897932385;
    var sin = global.Math.sin;
    var cos = global.Math.cos;
    var abs = global.Math.abs;
    var pow = global.Math.pow;
    var sqrt = global.Math.sqrt;

    var HEAPD = new global.Float64Array(heap);

    function mod(a, b) {
        a = +a;
        b = +b;
        return +(a % b); // TODO: check sign issues
    }
    
    function eccAnomaly(m, ec, maxiter, prec) {
        m = +m;
        ec = +ec;
        maxiter = maxiter|0;
        prec = +prec;

        var f = 0.0;
        var sinm = 0.0;
        var e0 = 0.0;

        m = mod(m, 2.0*PI);

        /* Iterative solution of Kepler equation with Newthon's method.
         */
        //var e0 = 0; /* TODO: initial value depends on eccentricity */
        /* Gary R. Smith.  A simple, efficient starting value for the
         * iterative solution of Kepler's equation.  // Celestial
         * Mechanics and Dynamical Astronomy.  Volume 19, Number 2 (1979)
         * DOI: 10.1007/BF01796088.
         */
        sinm = sin(m);
        e0 = +(m + ec * sinm / (1.0 - sin(m+ec) + sinm));

        do {
            f = +(e0 - ec * sin(e0) - m);
            e0 = +(e0 - f / (1.0 - ec * cos(e0)));
            maxiter = (maxiter - 1)|0;
        } while ((~(~maxiter) > ~(~0)) & (f > prec));

        if (~(~maxiter) < ~(~0)) {
            e0 = 1.0/0.0;
        }
        return +e0;
    }

    function stumpff(e2, eps) {
        e2 = +e2;
        eps = +eps;

        var c1 = 0.0;
        var c2 = 0.0;
        var c3 = 0.0;
        var add = 1.0;
        var n = 1.0;

        do {
            c1 = c1 + add;
            add = add / (2.0*n);
            c2 = c2 + add;
            add = add / (2.0*n + 1.0);
            c3 = c3 + add;
            add = add * -e2;
            n = n + 1.0;
        } while (abs(add) >= eps);

        HEAPD[3] = c1;
        HEAPD[4] = c2;
        HEAPD[5] = c3;
    }

    function parabolic(gm, t0, t, q, ec, maxiter, prec) {
        gm = +gm;
        t0 = +t0;
        t = +t;
        q = +q;
        ec = +ec;
        maxiter = maxiter|0;
        prec = +prec;

        var e2 = 0.0;
        var e20 = 0.0;
        var fac = 0.0;
        var k = 0.0;
        var tau = 0.0;
        var a = 0.0;
        var b = 0.0;
        var u = 0.0;
        var u2 = 0.0;
        var r = 0.0;

        fac = 0.5 * ec;
        k = sqrt(gm/(q*(1.0+ec)));
        tau = sqrt(gm)*(t-t0);

        do {
            e20 = e2;
            a = 1.5 * sqrt(fac/(q*q*q))*tau;
            b = pow(sqrt(a*a + 1.0)+a, 1.0/3.0);
            u = b - 1.0/b;
            u2 = u*u;
            e2 = u2*(1.0-ec)/fac;
            stumpff(e2, prec);
            fac = 3.0 * ec * HEAPD[5];
        } while (abs(e2-e20) >= prec)

        HEAPD[0] = q*(1.0 - u2*HEAPD[4]/fac);
        HEAPD[1] = q*sqrt((1.0+ec)/fac)*u*HEAPD[3];
        HEAPD[2] = 0.0;
    }

    return {
        eccAnomaly: eccAnomaly,
        stumpff: stumpff,
        parabolic: parabolic
    };
}

var MEM = new ArrayBuffer(4096);
var HEAPD = new Float64Array(MEM);

_accelerated = _acceleratedModule(window, {}, MEM);

/** Compute eccentric anomaly for elliptical orbit.
 * @param m Mean anomaly.
 * @param ec Eccentricity (ec < 1)
 * @param maxiter Optional number of iterations.
 */
StarJs.Kepler.eccAnomaly = function (m, ec, maxiter) {
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    e0 = _accelerated.eccAnomaly(m, ec, maxiter, prec);

    return (isNaN(e0)) ? null : e0;
};

/** Compute position of a body on elliptical orbit.
 * @param gm Gravity constant.
 * @param m Mean anomaly.
 * @param a Semi-major axis.
 * @param ec Eccentricity.
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
 @param gm Gravity constant.
 @param t0   time of pericenter
 @param t    time to calculate position for
 @param q    pericenter distance
 @param ec Eccentricity.
 @param maxiter Optional maximal number of iterations.
 */
StarJs.Kepler.parabolic = function (gm, t0, t, q, ec, maxiter) {
    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    _accelerated.parabolic(gm, t0, t, q, ec, maxiter, prec);    

    return new StarJs.Vector.Vector3(HEAPD[0], HEAPD[1], HEAPD[2]);
//     // res is position
//     r = q * (1+u2*c.c2*ec/fac);
//     var vel = new StarJs.Vector.Vector3(-k*res.y/r,
//                                         k*(res.x/r+ec),
//                                         0.0);
};

/** Compute hyperbolic anomaly. */
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
 * @param gm Gravity constant
 * @param t0 Time of pericenter
 * @param t Time
 * @param a Semi-major axis
 * @param e Eccentricity.
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
    @param gm   GM
    @param t0   time of pericenter
    @param t    time to calculate position for
    @param q    pericenter distance
    @param e    eccentricity
    @param pqr  Gauss' vector matrix
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
    @param omega longitude of the ascending node (radians)
    @param i     inclination (radians)
    @param w     argument of pericenter (radians)
 */
StarJs.Kepler.gaussVec = function (omega, i, w) {
    return StarJs.Vector.Matrix3.r_z(-omega).mult(StarJs.Vector.Matrix3.r_x(-i))
        .mult(StarJs.Vector.Matrix3.r_z(-w));
};
