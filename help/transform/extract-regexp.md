<!-- TITLE: Extract RegExp -->
<!-- SUBTITLE: -->

# Extract regular expression

Matches the specified regular expression against the content of the specified column.

To split column regular expression groups are used. Groups are added as new columns.

## Example

RegExp: (.\*)-(.\*)

Data:

| Planet  |   |
|---------|---|
| Earth-3 |   |
| Mars-4  |   |

Result:

| Planet  | r1    | r2 |
|---------|-------|----|
| Earth-3 | Earth | 3  |
| Mars-4  | Mars  | 4  |

To learn more about regular expressions, please visit
[http://www.regular-expressions.info](http://www.regular-expressions.info).

Click on the _keyboard_arrow_down_ icon in order to access the following commands:

* Filter Matching
* Filter Not Matching
* Select Matching
* Select None
