// Formatting data differently. Note that this does not change the actual values and does not affect sorting.
// Similar effect can be achieved by editing 'format' column tag via UI (Column | Properties or press F2).
// See available format strings below

let df = grok.data.demo.demog();
let view = grok.shell.addTableView(df);

df.col('age').meta.format = 'compact simple currency';
df.col('height').meta.format = 'scientific';
df.col('weight').meta.format = '#.0000';
df.col('started').meta.format = 'dd.MM.yyyy';

// NUMBER FORMATS (https://datagrok.ai/help/discover/tags#numbers)
// The usual formats, such as "0.000" or "#.0000", are supported.
// In addition, you can choose one of the following formats in the UI:
// |           Format           | Example |
// |----------------------------|---------|
// | int                        | 71      |
// | two digits after comma     | 70.50   |
// | four digits after comma    | 70.5000 |
// | max two digits after comma | 70.5    |
// | scientific                 | 7E1     |
// | money                      | $70.50  |
// | compact long               | 70.5    |
// | compact                    | 70.5    |
// | compact simple currency    | $70.50  |
// | percent                    | 7,050%  |
// | thousand separator         | 71      |
// | full precision             | 70.5    |


// DATETIME FORMATS (see https://datagrok.ai/help/discover/tags#datetime)
// The standard formats representing date and time include:
// |         Format          |          Example            |
// |-------------------------|-----------------------------|
// | MM/dd/yyyy HH:mm:ss.fff | 07/21/1991 00:00:00.000     |
// | M/d/yyyy                | 7/21/1991                   |
// | dd.MM.yyyy              | 21.07.1991                  |
// | M/d/yyyy h:mm tt        | 7/21/1991 12:00 AM          |
// | M/d/yyyy h:mm:ss tt     | 7/21/1991 12:00:00 AM       |
// | yyyy-MM-dd              | 1991-07-21                  |
// | dddd, MMMM d, yyyy      | Sunday, July 21, 1991       |
// | MMM d, yyyy             | Jul 21, 1991                |
// | h:mm:ss                 | 12:00:00                    |
// | h:mm:ss.fff             | 12:00:00.000                |
// | relative                | 29 years ago                |
// | auto                    | Jul 21, 1991                |

// Supported units of time:
// | Symbol |          Meaning                 |          Example           |
// |--------|----------------------------------|----------------------------|
// | yy     | Year without the century         | 00, 01, ..., 20, ..., 99   |
// | yyyy   | Year with the century            | 0001, ..., 2020, ..., 9999 |
// | M      | Month                            | 1, 2, 3, ..., 12           |
// | MM     | Zero-padded month                | 01, 02, 03, ..., 12        |
// | MMM    | Abbreviated month name           | Jan, Feb, Mar, ..., Dec    |
// | MMMM   | Full month name                  | January, ..., December     |
// | d      | Day of the month                 | 1, 2, 3, ..., 31           |
// | dd     | Zero-padded day                  | 01, 02, 03, ..., 31        |
// | ddd    | Abbreviated weekday name         | Mon, ..., Fri, Sat, Sun    |
// | dddd   | Full weekday name                | Monday, ..., Sunday        |
// | h      | Hour (12-hour clock)             | 1, 2, 3, ..., 12           |
// | hh     | Zero-padded hour (12-hour clock) | 01, 02, 03, ..., 12        |
// | H      | Hour (24-hour clock)             | 0, 1, 2, ..., 23           |
// | HH     | Zero-padded hour (24-hour clock) | 00, 01, 02, ..., 23        |
// | m      | Minute                           | 0, 1, 2, ..., 59           |
// | mm     | Zero-padded minute               | 00, 01, 02, ..., 59        |
// | s      | Second                           | 0, 1, 2, ..., 59           |
// | ss     | Zero-padded second               | 00, 01, 02, ..., 59        |
// | f      | Second fraction (1 digit)        | 0, 1, ..., 9               |
// | ff     | Second fraction (2 digits)       | 00, 01, ..., 99            |
// | fff    | Second fraction (3 digits)       | 000, 001, ..., 999         |
// | tt     | 12-hour periods                  | AM, PM                     |
