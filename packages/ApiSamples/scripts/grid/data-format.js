// Formatting data differently. Note that this does not change the actual values and does not affect sorting.

let view = grok.shell.addTableView(grok.data.testData('demog', 5000));

// similar effect can be achieved by editing 'format' column tag via UI (Column | Properties).
// see available format strings below
view.grid.col('age').format = 'compact simple currency';
view.grid.col('weight').format = '#.0000';

// NUMBER FORMATS:
// The usual formats, such as "0.000" or "#.0000"
// ALSO, the following formats are understood:
// two digits after comma
// four digits after comma
// max two digits after comma
// scientific
// money
// compact long
// compact
// compact simple currency
// percent
// thousand separator
// full precision
// money($)


// DATE-TIME FORMATS (see https://pub.dev/documentation/intl/latest/intl/DateFormat-class.html)
// ICU Name                   Skeleton
// --------                   --------
//     DAY                          d
// ABBR_WEEKDAY                 E
// WEEKDAY                      EEEE
// ABBR_STANDALONE_MONTH        LLL
// STANDALONE_MONTH             LLLL
// NUM_MONTH                    M
// NUM_MONTH_DAY                Md
// NUM_MONTH_WEEKDAY_DAY        MEd
// ABBR_MONTH                   MMM
// ABBR_MONTH_DAY               MMMd
// ABBR_MONTH_WEEKDAY_DAY       MMMEd
// MONTH                        MMMM
// MONTH_DAY                    MMMMd
// MONTH_WEEKDAY_DAY            MMMMEEEEd
// ABBR_QUARTER                 QQQ
// QUARTER                      QQQQ
// YEAR                         y
// YEAR_NUM_MONTH               yM
// YEAR_NUM_MONTH_DAY           yMd
// YEAR_NUM_MONTH_WEEKDAY_DAY   yMEd
// YEAR_ABBR_MONTH              yMMM
// YEAR_ABBR_MONTH_DAY          yMMMd
// YEAR_ABBR_MONTH_WEEKDAY_DAY  yMMMEd
// YEAR_MONTH                   yMMMM
// YEAR_MONTH_DAY               yMMMMd
// YEAR_MONTH_WEEKDAY_DAY       yMMMMEEEEd
// YEAR_ABBR_QUARTER            yQQQ
// YEAR_QUARTER                 yQQQQ
// HOUR24                       H
// HOUR24_MINUTE                Hm
// HOUR24_MINUTE_SECOND         Hms
// HOUR                         j
// HOUR_MINUTE                  jm
// HOUR_MINUTE_SECOND           jms
// HOUR_MINUTE_GENERIC_TZ       jmv
// HOUR_MINUTE_TZ               jmz
// HOUR_GENERIC_TZ              jv
// HOUR_TZ                      jz
// MINUTE                       m
// MINUTE_SECOND                ms
// SECOND                       s//