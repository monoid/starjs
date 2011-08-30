function formatDate(date) {
    return date.year+'.'+date.month+'.'+date.day+' '+date.hour+':'+date.minute+':'+Math.floor(date.second);
}

function init(id) {
    var output = document.getElementById(id);

    var mjNow = StarJs.Time.time2mjd(new Date()); // Current Modified Julian day

    var R2D = StarJs.Math.RAD2DEG, D2R = StarJs.Math.DEG2RAD;

    var text = '';

    /* Elements from
       http://scully.cfa.harvard.edu/cgi-bin/returnprepeph.cgi?d=c&o=CK10X010

    C/2010 X1 (Elenin)
Epoch 2011 Aug. 27.0 TT = JDT 2455800.5
T 2011 Sept. 10.7227 TT                                 MPC
q   0.482465             (2000.0)            P               Q
z  -0.000057       Peri.  343.8056      +0.6023511      +0.7980000
 +/-0.000002       Node   323.2267      -0.7287563      +0.5399437
e   1.000028       Incl.    1.8392      -0.3257108      +0.2676879
From 2058 observations 2010 Dec. 10-2011 Aug. 1, mean residual 0".5.
    */
    /*
    CK10X010  2011 09 10.7226  0.482465  1.000028  343.8056  323.2267    1.8392  20110827  10.0  4.0  C/2010 X1 (Elenin)                                       MPC 75713
    */
    var params = {
        t0: 2455800.5 - 2400000.5,
        q: 0.482465,
        z: -0.000057,
        e: 1.000028,
        peri: 343.8056*D2R,
        node: 323.2267*D2R,
        incl: 1.8392*D2R
    };

    // Gauss vectors for comet's orbit
    var pqr = StarJs.Kepler.gaussVec(params.node, params.incl, params.peri);
    
    // Equatorial prcession matrix.
    // It changes slowly, so we use one matrix for whole year
    var precessionMatrix = StarJs.Coord.precessionEquMatrix(0, StarJs.Time.mjd2jct(mjNow));

    for (i = 0; i < 10; ++i) {
        // time i weeks from now
        var mjdi = mjNow+i*7;
        var time = StarJs.Time.mjd2dt(mjdi);

        // Julian centuries after J2000.0.
        var jct = StarJs.Time.mjd2jct(mjdi);

        // Heliocentric ecliptic coordinates (StarJs.Math.Vector) of
        // comet at that time
        var keplerCoord = StarJs.Kepler.keplerPos(StarJs.Solar.GM,
                                                  params.t0,
                                                  mjdi,
                                                  params.q,
                                                  params.e,
                                                  pqr);
        // Heliocentric ecliptic coordinates (x,y,z) of Earth at that
        // time (epoch 2000.0)
        var earthKeplerCoord = StarJs.Solar.BODIES.Earth.keplerCoord(jct);

        // Earth-centric ecliptic coordinates of the comet
        var cometCoord = keplerCoord.sub(earthKeplerCoord);

        var distance = cometCoord.len();
        
        // Earth-centric equatorial coordinates of the comet
        var cometECoord = StarJs.Coord.ecl2equMatrix(jct).apply(cometCoord);

        var cometPolar = new StarJs.Vector.Polar3(cometECoord);

        // Output text
        text += formatDate(time) + ' | ' + (cometPolar.phi*R2D).toFixed(2) + ' | ' + (cometPolar.theta*R2D).toFixed(2) + ' | ' + distance.toFixed(2) + '\n';
    }

    output.innerHTML = text;
}
