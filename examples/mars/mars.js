var lon = 83, lat = 54; // Approximate coordinates of Novosibirsk

function formatDate(date) {
    return date.year+'.'+date.month+'.'+date.day+' '+date.hour+':'+date.minute+':'+Math.floor(date.second);
}

function init(id) {
    var output = document.getElementById(id);

    var mjNow = StarJs.Time.time2mjd(new Date()); // Current Modified Julian day

    var R2D = StarJs.Math.RAD2DEG, D2R = StarJs.Math.DEG2RAD;

    var text = '';

    // Equatorial prcession matrix.
    // It changes slowly, so we use one matrix for whole day
    var precessionMatrix = StarJs.Coord.precessionEquMatrix(0, StarJs.Time.mjd2jct(mjNow));

    for (i = 0; i < 24; ++i) {
        // time i hours from now
        var mjdi = mjNow+i/24;
        var time = StarJs.Time.mjd2dt(mjdi);

        // Julian centuries after J2000.0.
        var jct = StarJs.Time.mjd2jct(mjdi);

        // Heliocentric ecliptic coordinates (StarJs.Math.Vector) of
        // Mars at that time (epoch J2000.0)
        var marsKeplerCoord = StarJs.Solar.BODIES.Mars.keplerCoord(jct);
        // Heliocentric ecliptic coordinates (x,y,z) of Earth at that
        // time (epoch 2000.0)
        var earthKeplerCoord = StarJs.Solar.BODIES.Earth.keplerCoord(jct);

        // Earth-centric ecliptic coordinates of Mars (epoch J2000)
        var marsCoord = marsKeplerCoord.sub(earthKeplerCoord);

        // Earth-centric equatorial coordinates of Mars (epoch J2000)
        var marsECoord = StarJs.Coord.ecl2equMatrix(jct).apply(marsCoord);

        // Convert to current Epoch.  Actually, the difference between
        // J2000 and 2011 is minimal, so this step can be
        // skipped.  It is included here just for completeness.
        marsECoord = precessionMatrix.apply(marsECoord);

        var marsPolar = new StarJs.Vector.Polar3(marsECoord);

        var gmst = StarJs.Time.gmst(mjdi);

        // Object with fields 'h' (height) and 'az' (azimuth, 0 is
        // South, not North, as it is usual in astronomy) of the body
        var skyCoord = StarJs.Coord.equ2hor(marsPolar.theta, (gmst-marsPolar.phi+lon*D2R), lat*D2R);

        // Output text
        text += formatDate(time) + ' | ' + marsPolar.phi*R2D + ' | ' + marsPolar.theta*R2D + ' | '+skyCoord.h*R2D+' | ' + skyCoord.az*R2D + '\n';
    }

    output.innerHTML = text;
}
