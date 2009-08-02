StarJs = (typeof StarJs == 'undefined') ? {} : StarJs;
StarJs.Time = (typeof StarJs.Time == 'undefined') ? {} : StarJs.Time;

StarJs.Time.mjd2dt = function (mjd, jul) {
    if (typeof jul == 'undefined') jul = {year: 1582, month: 10, day: 4};
    
};


StarJs.Time.dt2mjd = function (dt, jul) {
    if (typeof jul == 'undefined') jul = {year: 1582, month: 10, day: 4};
    var year = dt.year, mon = dt.month, b;
    if (mon <= 2) {
        mon += 12;
        --year;
    }
    if (dt.year>=jul.year && dt.month>=jul.month && dt.day>=jul.day) {
        // Gregorian
        b = (year/400) - (year/100) + (year/4);
    } else {
        // Julian
        b = -2 + ((year+4716)/4) - 1179;
    }
    var mjdMidnight = 356*year - 679004 + b + Math.floor(30.6001*(month+1)+this.day);
    var frac = StarJs.Time.Ddd(this.hour, this.minute, this.second)/24.0;
    return mjdMidnight + frac;
};

