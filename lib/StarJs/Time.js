/*global StarJs: false */
/*jslint onevar: false */
/**
 * @depends ../StarJs.js
 * @depends Math.js
 */
StarJs.Time = (typeof StarJs.Time === 'undefined') ? {} : StarJs.Time;

StarJs.Time.DEFAULT_JULIAN_DATE = {year: 1582, month: 10, day: 4};
StarJs.Time.DEFAULT_JULIAN_JD = 2299161;

StarJs.Time.time2mjd = function (t) {
    if (typeof t !== 'number') {
        t = t.getTime();
    }
    return t/86400000.0 + 40587.0;
}

StarJs.Time.mjd2dt = function (mjd, jul) {
    if (typeof jul === 'undefined') {
        jul = StarJs.Time.DEFAULT_JULIAN_JD;
    }
    var a = Math.floor(mjd) + 2400001, b, c, d, e, f, day, mon, year, hour, t;
    if (a < jul) {
        // Julian calendar
        b = 0;
        c = a + 1524;
    } else {
        // Grigorian calendar
        b = Math.floor((a - 1867216.25) / 36524.25);
        c = a + b - Math.floor(b / 4) + 1525;
    }
    d = Math.floor((c - 122.1) / 365.25);
    e = 365 * d + Math.floor(d / 4);
    f = Math.floor((c - e) / 30.6001);
    day = c - e - Math.floor(30.6001 * f);
    mon = f - 1 - 12 * Math.floor(f / 14);
    year = d - 4715 - Math.floor((7 + mon) / 10);
    hour = 24.0 * (mjd - Math.floor(mjd));
    t = StarJs.Time.hour2hms(hour);
    t.year = year;
    t.month = mon;
    t.day = day;
    return t;
};

StarJs.Time.hms2hour = function (h, m, s) {
    return h + (m / 60.0) + (s / 3600.0);
};

StarJs.Time.hour2hms = function (hour) {
    var dms = new StarJs.Math.Dms(hour);
    return {hour: dms.degree, minute: dms.minute, second: dms.second};
};

StarJs.Time.dt2mjd = function (dt, jul) {
    if (typeof jul === 'undefined') {
        jul = StarJs.Time.DEFAULT_JULIAN_DATE;
    }

    var year = dt.year, mon = dt.month, b;
    if (mon <= 2) {
        mon += 12;
        year -= 1;
    }
    if (year <= jul.year && mon <= jul.month && dt.day <= jul.day) {
        // Julian
        b = -2 + Math.floor((year + 4716) / 4) - 1179;
    } else {
        // Gregorian
        b = Math.floor(year / 400) - Math.floor(year / 100) + Math.floor(year / 4);
    }
    var mjdMidnight = 365 * year - 679004 + b + Math.floor(30.6001 * (mon + 1)) + dt.day;
    var frac = StarJs.Time.hms2hour(dt.hour, dt.minute, dt.second) / 24.0;
    return mjdMidnight + frac;
};

StarJs.Time.mjd2jct = function (mjd) {
    return (mjd - 51544.5) / 36525.0;
};

StarJs.Time.gmst = function (mjd) {
    /* TODO: move to global */
    var SECS = 86400; // 24*60*60 -- number of seconds in day;
    var mjd0 = Math.floor(mjd), ut = SECS * (mjd - mjd0), t0, t, gmst;
    t0 = StarJs.Time.mjd2jct(mjd0);
    t = StarJs.Time.mjd2jct(mjd);
    gmst = 24110.54841 + 8640184.812866 * t0 + 1.0027379093 * ut +
	(0.093104 - (6.2e-6) * t) * t * t;
    return StarJs.Math.PI2 / SECS * StarJs.Math.mod(gmst, SECS);
};


