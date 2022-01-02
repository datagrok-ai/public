<!-- TITLE: Dataframe -->

# DataFrame

DataFrame is a tabular structure with strongly-typed columns of various types. A dataframe class 
[DG.DataFrame](https://datagrok.ai/js-api/classes/dg.dataframe)
is used in virtually any Datagrok extension or application. It operates via a columnar in-memory data engine
which Datagrok implemented from scratch to support highly efficient operation with data in a modern browser
together with fast in-browser data visualizations. 

## Dataframe JavaScript API

Dataframe stores data as list of columns. Constructing, modifying and efficiently accessing data points are
embodied in both [`DG.Column`][101] and [`DG.DataFrame`][100] classes. Event handling, visual aspects of working
with dataframes, fast column selection, handy construction methods and row-based access are provided in
[`DG.DataFrame`][100]. Instances of [`DG.ColumnList`][099], [`DG.RowList`][098],
[`DG.Row`](/js-api/classes/dg.row) and [`DG.Cell`](/js-api/classes/dg.cell) are used as related properties
or functions return values of `DG.Column` and `DG.DataFrame`.

## Dataframe design

A Datagrok dataframe reminds of functionally similar structures in
[Python](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html) and
[R](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/data.frame).

In comparison, Datagrok implementation is considerably optimized. While it isn't possible to control the browser
entirely, Datagrok in-memory engine is designed in such way that all data-related operations, such as aggregations
or statistical computations, are performed efficiently on modern machines.

Here are a few of unique dataframe design features:

* Built to allocate little memory space, utilizes adaptive bit storage
* Operates on raw data instead of JavaScript wrapping objects
* Utilizes manual memory management with the column-based layout
* Doesn't have blocking structures, works well in multithreaded environments
* Made to work entirely in the browser, but same dataframe code also runs on Datagrok servers, e.g. for serialization
* Uses custom serialization codecs for efficiency

Datagrok [visualizations][102] sit on top of the dataframe, that's one of the key reasons for their high performance.

## `DG.Column`

Columns support the types specified in [`DG.COLUMN_TYPE`][103] enum: `STRING`, `INT`, `FLOAT`, `BOOL`,
`DATE_TIME`, `BIG_INT`, `QNUM`, `DATA_FRAME` and `OBJECT`. Find the details of handling these types
[further][105].

### Column properties

Every column has:

* a `.name`: set it at construction, change at runtime
* a `.length`: it's a read-only property
* a .`type`: one of `DG.COLUMN_TYPE`, possible to change later via `DG.DataFrame` class
* a parent `.dataFrame`, if the column is a part of some dataframe
* `.temp` and `.tags`: provide access to [special proxy classes][114] for named temporary values and named tags
  which can be associated with any column
  
`.tags` play key role in platform extensibility, allowing to introduce custom behavior and new data types.

A column which is already constructed, either by an [initialization from a list][117], or being part of some
dataframe from its creation, isn't a subject to changing its length, if taken standalone. This is due to its
memory-optimized nature. However, the columns length is changed with [adding rows][] to the embodying dataframe.
  
#### Semantic type

A property `.semType` is a string value representing a tag which associates a column data with a logical type.

An underlying raw data type of a molecule or a peptide is usually just a `string`. Semantic typing allows expanding
the platform with such logical types recognition and understanding with a help of 3-rd party Datagrok packages.
For example, [Chem][112] package adds visualizing and handling molecules, including [detecting][111] that a `string`
column is of type `Molecule`, computing molecular fingerprints, searching for substructures, ranking by similarity.

The property `.semType` is a getter to a [tag][114] named `'quality'`. This distinct getter is there due to a key
role of this tag.

Learn more about semantic types in this article: [Link][109].

### Construct a column

* The most common way is to use a method `.fromList`, explicitly specifying a type and a column name:  
  ```let col =  DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Column Name', [1, 2, 3])```
* The method `.fromStrings` recognizes and automatically assigns a type:  
  ```let col =  DG.Column.fromStrings('Column Name', ['3.14', '2.71']); // col.type === DG.COLUMN_TYPES.FLOAT```
* To create a column with [`NULL` values][116] of a pre-specified length, use `.fromType` or concrete types shortcuts:  
  ```let col =  DG.Column.fromType(DG.COLUMN_TYPE.INT, 'Name', 3); // col.get(0) === DG.INT_NULL```  
  ```let col =  DG.Column.string('Name', 5); // col.get(2) === ""```
  
The column, once constructed, may later be [added to a dataframe]().

### Manipulate column values

#### Access and modify column values

* A method `.get` is passed an index `i` to return an `i`-th value: `const value = column.get(idx)`  
* A method `.set` sets `i`-th value to `x`: `column.set(i, x)`

A pair of methods `.getString`/`.setString` work similarly, but with formatting and parsing:
* `.getString` returns a value converted to a string taking into account a column's
  [tag `format`]()
* `.setString` attempts to set an `i`-th value by converting a provided string to the corresponding
  strongly-typed value, returns `true` if a text was successfully parsed and set, otherwise `false`

```javascript
const table = grok.data.demo.demog();
const column = table.columns.byName('weight');
column.tags.format = '#.0000';
grok.shell.info(column.getString(217)); // displays '108.7208'
grok.shell.info(column.setString(15, '3.1415')); // displays 'true'
grok.shell.info(column.setString(16, 'non-number')); // displays 'false'
```

Learn more about supported formats in this article: [Link][108].

<!-- TODO: Explain `notify` -->

#### Refresh the state after modifying the column

Some column associated structures, such as [categories][137], won't be automatically updated
after the column data modification with `.setString`, `.set` and similar. This is intentional
for performance reasons. After the modification is meant to be completed, it is possible
to refresh the column `col` representation, including categories:

```
col.compact();
```

This call compacts the internal column representation. Currently, it only affects string
columns where values were modified

#### Initialize values with a function

It's handy to set all column values in a single batch. This is possile with `.init` function:

```javascript
let df = DG.DataFrame.fromCsv(
  `x, y, s
   1, 2, Point A
   3, 5, Point B
`);
df.col('x').init(2); // df.col('x').get(1) === 2
df.col('s').init('No comment'); // df.col('s').get(0) === 'No comment'
df.columns.addNewInt('counter').init((i) => i * 3); // df.col('counter').get(1) === 3
grok.shell.addTableView(df);
```

#### Access performance

To learn the typical times it takes to run various column access patterns, run
[this example]((https://public.datagrok.ai/js/samples/data-frame/performance/access)). In summary,
it is advised to explicitly get a column and its length  as separate variables before looping through:

```javascript
const t = grok.data.demo.demog(100000);
const column = t.columns.byName('age');
const rowCount = t.rowCount;
let sum = 0;
for (let i = 0; i < rowCount; i++)
  sum += column.get(i);
```

Doing a `.byName` or a `t.rowCount` call as part of the loop shall incur up to 20x overhead.
A `column.values()` iterator is also available, it is 2-3 times slower than this snippet.

If the fastest access is required for numerical column data, which usually happens in computing new values
atop a column, accessing data with a result of calling `.getRawData()` is advised:

```javascript
const table = grok.data.demo.demog(100000);
const column = table.columns.byName('age');
const array = column.getRawData();
const rowCount = column.length;
let sum = 0;
for (let i = 0; i < rowCount; i++)
  sum += array[i];
```

It is 4-6 times faster than a snippet from above with explicit `column` and `rowCount`.

`.getRawData` returns a `Float32Array` for `DG.COLUMN_TYPE.FLOAT`, `Int32Array` for `DG.COLUMN_TYPE.INT`.

### Column data types

#### Special values

In popular languages like JavaScript, Python and R, special values representing missing values, unset values,
or values for the undefined results (such as not valid numbers), are known as `None`, `NaN` (Not A Number),
`null`, or `undefined`. Dataframe's type system implements special handling of such values in
a way great for performance yet still allowing for these kinds of values to be present.

<!-- TODO: How is this handled when we work with these in server functions? -->

##### `NULL` values for numeric types

For performance reasons, Datagrok type system implements an equivalent of an empty value for numeric types
`DG.COLUMN_TYPE.INT` and `DG.COLUMN_TYPE.FLOAT`. There are [special constants][104] `DG.INT_NULL = -2147483648` and
`DG.FLOAT_NULL = 2.6789344063684636e-34`, respectively.

To set a value as a Datagrok `NULL`, pass either a JavaScript `null` or a corresponding constant:

```javascript
let col = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Name', [314, null, DG.INT_NULL, 143]);
col.set(0, null);
grok.shell.info(col.get(0)); // shows '-2147483648'
grok.shell.info(col.get(1)); // shows '-2147483648'
grok.shell.info(col.get(2)); // shows '-2147483648'
```

When checking if a value is a Datagrok `NULL`, be cautious of using `=== null`, as the value isn't a JavaScript `null`.
Instead, use a `DG.Column`'s method `.isNone(i)` to check if the `i`-th element is Datagrok `NULL`. Continuing with
the previous example:

```javascript
grok.shell.info(col.isNone(2)); // shows 'true'
grok.shell.info(col.isNone(3)); // shows 'false'
```

##### `NULL` values in table views 

* In the [table view][107], Datagrok `NULL` values are seen as empty cells
* Entering an empty string in a numeric column makes a cell valued to a Datagrok `NULL`
* If an empty string is passed to a [`.setString`][106], a Datagrok `NULL` value shall be assigned

##### `NULL` values for strings

An empty string is an equivalent of a Datagrok `NULL`:

```javascript
let col = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Name', [null, '']);
grok.shell.info(col.isNone(0)); // shows 'true'
grok.shell.info(col.isNone(1)); // shows 'true'
``` 

##### `undefined` and `NaN` values

* Passing `undefined` as a column item value is equivalent to passing a `null`
* [`NaN`](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/NaN)
in JavaScript represents a value which is not a valid number: this value isn't currently supported


<!-- TODO: How's that for other types? -->

#### String, integer, float, boolean

Enum values: `DG.COLUMN_TYPE.STRING`, `DG.COLUMN_TYPE.INT`, `DG.COLUMN_TYPE.FLOAT`, `DG.COLUMN_TYPE.BOOL`

These are regular JavaScript types. The only difference is handling `NULL`-values as described above.

#### Datetime

Enum value: `DG.COLUMN_TYPE.DATE_TIME`

In future releases, `DateTime` type becomes a wrapper around a [dayJs][118] value.

#### Bigint

Enum value: `DG.COLUMN_TYPE.BIG_INT`

The type for working with integers that do not fit into 53 bits.

`BIG_INT` won't be rendered by a [grid viewer][119] and won't become part of [aggregations][120]. It is introduced
for compatibility with `BIG_INT` types in databases.

#### Dataframe

Enum value: `DG.COLUMN_TYPE.DATA_FRAME`

[This example](https://public.datagrok.ai/js/samples/data-frame/advanced/data-frames-in-columns) demonstrates
how to store dataframes as column values.  

#### Qualified number

Enum value: `DG.COLUMN_TYPE.QNUM`

Qualified numbers, or QNums, are typically used to represent measurements when the exact value is not known,
but it is known that it is either less or greater than some value. The examples of such values are `<3.5` or `>5E-6`.

To keep the performance of this data type on par with `INTEGER` and `FLOAT`, `QNUM` was based on the 64-bit floating
point number where two least significant bits of mantissa are reserved for the qualifier. This storage structure
comes with ability to efficiently store the numbers and perform arithmetic operations on them without having to incur
costly packing/unpacking. The qualifier part is valued to one of: `1` (`LESS`), `2` (`EXACT`), or `3` (`GREATER`).
This way, `QNums` are compared using regular floating point number comparison.

There is no special internal data type, as the values are 64-bit IEEE754 floats. In most cases the isn't need to
check whether the value is qualified or not, since the result for the most operations (rendering, using for
visualization, arithmetic operations) will be the same. However, in cases it does matter the programmer
has to keep track of whether the values are qualified, and pack/unpack accordingly. This is 
achieved with [`DG.Qnum` class][] class containing helper methods:

```javascript
let col = DG.Column.qnum('col', 3);
col.set(0, DG.Qnum.greater(5));
col.set(1, DG.Qnum.exact(5));
col.set(1, DG.Qnum.less(5));
grok.shell.info(DG.Qnum.getQ(col.get(0)) == DG.QNUM_GREATER); // shows `true`
```

#### Object

Enum value: `DG.COLUMN_TYPE.OBJECT`

This may be any JavaScript object. Pass it as usual.

### Statistics

The properties `.min` and `.max` of `Column` return and cache the minimum and maximum numeric values. Using
the internal [version][128] of the column, the caches for `.min` and for `.max` are automatically invalidated
whenever the column is modified.

Computing popular statistics, such as average `.avg`, median `.med` or standard deviation `.stdev`,
is more efficient when done altogether. If that's the case you need, it's possible to get all such statistics
using a `.`, which returns an istance of [`Stats`][130] class.

Run `Stats` example: [Link][129].

### Numerical and categotical columns

Currently, columns of types `BOOL` and `STRING` are considered as categorical, s.t. they represent values from a
discrete set. Columns of types `INT` and `FLOAT` are considered as numerical, s.t. they represent values from
a continuous set. Later, this will be adjusted s.t. whether column is categorical or numerical will be automatically
detected based on its content analysis.

This property mostly affects the way columns are treated in the UI. For instance, this property affects how
[color coding][124] works.

### Work with categories

For a `string` column, the column collection `.categories` returns a list of categories in an alphabetical sort.
For instance, for a column of strings `['a', 'b', 'b', 'c']` the categories list will be `['a', 'b', 'c']`.

It is possible to introduce a custom category order, as [this example][127] demonstrates. That is useful if
the alphabetical order does not reflect a natural order, for example in `['Low', 'Medium', 'High']`.

### Versioning

A column increments an integer version number whenever its contents are changed. The `.version` property
returning this number is useful when invalidating cached computation results stored as columns or dataframes. Check
[this example][126] to see how column versioning works.

## `DG.DataFrame`

### Construct

Dataframes may be obtained through the JavaSript or TypeScript code in various ways:

* a new dataframe constructed from predefined data: [Link][122]
* a new dataframe from precreated Datagrok columns: [Link][122]
* a table already being rendered by a table view: [Link][]
* a dataframe constructed from a file in a file share: [Link][]
* a CSV file uploaded to a browser
* a dataframe returned by a script: [Link][]
* as calculated on the flight for aggregations: [Link][]

#### Construct from in-place content

A variety of in-place construction methods are available:

* Construct from array of explicitly created columns, all same length: [link](
  https://public.datagrok.ai/js/samples/data-frame/construction/create-from-arrays)
* Specify number of rows in a dataframe at creation and later adding columns with this number of items: [link](
  https://public.datagrok.ai/js/samples/data-frame/construction/create-from-columns)
* Create from a string containing a CSV: [link](
  https://public.datagrok.ai/js/samples/data-frame/construction/create-from-csv-format)
* Create from a string containing a JSON: [link](
  https://public.datagrok.ai/js/samples/data-frame/construction/create-from-json)
* Create from JavaScript objects: [link](
  https://public.datagrok.ai/js/samples/data-frame/construction/create-from-objects)
* Create from JavaScript typed arrays: [link](
  https://public.datagrok.ai/js/samples/data-frame/construction/create-from-typed-arrays)

#### Load a demo dataset

Datagrok comes with a hundred demo datasets [available as a folder][121] in the platform interface. It's possible
to load these datasets programmatically using Datagrok API methods:

* `let table = grok.data.demo.demog()` loads the first `10'000` rows of the demographic dataset we used a lot
   in this article. It's possible to specify a desired number of rows as a parameter, and use one of
   `.biosensor`/`.wells`/`.geo` functions instead of `.demog` to load corresponding data
* `let data = grok.data.demo.randomWalk(rows, cols)` generates random walk data
* `await table = grok.data.getDemoTable('geo/earthquakes.csv')` allows loading any demo dataset from the specified
  demo datasets file share by its path (the method is asynchronous)

### Access

#### Access columns

Dataframes' columns are accessed and manipulated as a collection by an instance of [`DG.ColumnList`][099] class
called `.columns`. Calling `.columns` methods is often a starting point to doing anything with a given dataframe:

```javascript
let d = grok.data.demo.demog();
let len = d.columns.length;
ageColumn = d.columns.byName('age');
ageColumn = d.columns['age'];
ageColumn = d.columns.byIndex(3);
ageColumn = d.columns[3];
ageColumn = d.columns.age;
```

Note that access by name is case-insensitive.

There are also two shortcuts to get the column by name:

```javascript
ageColumn = d.col('age');    // won't throw if column isn't found
ageColumn = d.getCol('age'); // will throw if column isn't found
```

Technically, `.columns` object is an instance of a special JavaScript wrapper around an instance of `ColumnList`.
While the `ColumnList` class itself provides for the access methods such as `.byIndex` or `.length`, this wrapper
adds support for square brackets, access by explicit column name, like `.d.columns.age`, and other syntactic
sugaring of JavaScript such as iterating columns with a `for (... of ...)` loop:

```javascript
let df = grok.data.demo.demog();
for (let col of df.columns) {
  grok.shell.info(col.name);
}
```

In addition to regular access to columns by index and name, there's a group of methods covered in
[this example][123] for retrieving columns `Iterable` by a specific property:

* With pre-specified [tags][114]:  
  ```javascript
  let demog = grok.data.demo.demog();
  demog.getCol('race').setTag('tag1', 'value1').setTag('tag2', 'value2');
  // undefined or null means that any value passes
  for (let column of demog.columns.byTags({'tag1': 'value1', 'tag2': undefined}))
    grok.shell.info(column.name);
  ```
* With [categorical or numerical][125] columns:
  ```javascript
  let demog = grok.data.demo.demog();
  for (let column of demog.columns.categorical) grok.shell.info(column.name);
  for (let column of demog.columns.numerical) grok.shell.info(column.name);
  ```
  
#### Access rows

As [`DG.ColumnList`][099] provides for column access and manipulation through a `.columns` property, its counterpart
[`DG.RowList`][098] with an instance called `.rows` allows for rows access and manipulation.

As dataframe is columnar-optimized, it is preferred to work with it in columns-then-rows rather than
rows-then-columns fashion. Typically, accessing a dataframe value is based on a column and then row,
simply by obtaining a reference to a relevant column and then accessing its elements by their indices.
This won't create any intermediate objects besides the ones already existing as part of the dataframe.

If a row-by-row access is desired, that is possible with additional rows materialization. This brings memory
and performance overheads, though it may be a convenient tool in cases where dataframes are small — no more than
a 100 of rows. An example of iterating through a dataframe row-by-row:

```javascript
const df = grok.data.demo.demog(5);
for (let row of df.rows) {
  grok.shell.info(row.idx);
}
```

To access rows' values in such iteration, use a collection of [`DG.Cell`][097] named `.cells`, which is
another wrapper around a values stored at a columns' offset:

```javascript
const df = grok.data.demo.demog(10);
for (const row of df.rows) {
  for (const cell of row.cells) {
    grok.shell.info(cell.value);
  }
}
```

<!-- TODO: Get a row at an index? -->

### Manipulate

#### Manipulate columns

To add, remove or replace columns in a dataframe:
* for an existing instance `column` of `Column`:
  * `.add(column)` adds it as a last column of the dataframe
  * `.insert(column, index)` adds it _before_ the column position which is currently at position `index`,
    starting count from `0`. The default value of `index` is `null`, which corresponds to adding the column
    at the end of the columns list
* `.addNew` adds a new empty column at the end of the column list
* `.addNew<TypeName>(name)` (`.addNewInt`, `.addNewString` etc) adds a new empty column of a [type][105]
  `<TypeName>` at the end of the column list
* `.addNewVirtual` adds a new [virtual column][135] at the end of the column list
* `.replace(oldColumn, newColumn)`  substitutes inside a dataframe an `oldColumn` passed as an object
  with a `newColumn`. The column remains avaialble by `oldColumn`, but it is no longer part of the dataframe
* 

Examples:
* run "Add columns": [Link][133]
* run "Manipulate columns": [Link][134]

#### Add a column by a formula

The `.columns` property of `DataFrame` provides for a method to create a column by a mathematical formula,
specified with a string, which may also involve any [function][132] registered within the platform.
For example, in a dataframe `df` with columns `X` and `Y` it's possible to add a new column `Z`, specified as follows:

```javascript
df.columns.addNewCalculated('Z', 'Sin({X}+${Y})/2');
```

The column type shall be deduced automatically, or may be specified explicitly as one of this method arguments.

Run "Add calculated columns" example: [Link][131].

#### Manipulate rows

While it is not generally advised to access the dataframe values row-by-row rather than column-then-offset,
it is often useful to expand or shrink dataframes row-by-row. The `.rows` property enables a set of methods for this:

* `addNew([value1, ..., valueN])` appends a row to the end of the dataframe, filling it with values
  specified in the passed array, there should be as many values as there are columns. If a null value
  is passed as an array, an empty row shall be created with the values of [Datagrok `None`][116]

* `.setValues(rowIdx, [value1, ..., valueN])` will set values in the row at the specified offset according
  to the values in the array
 
* `.insertAt(idx, count)` inserts `count` new empty rows _before_ a row `idx`

* `.removeAt(idx, count)` removes `count` rows, _starting_ from a row `idx` and further

#### Extend

#### Join

#### Aggregate

#### Pivot

### Virtual columns

Consider a person's `weight`, `height` and a derived value of a mass index `BMI`, which is computed as
`weight / height^2`. It is handy to access a `BMI` value as if it is stored in the dataframe along with other
column values. However, it would not be practical to pre-compute these values and actually store in the dataframe, and
also be inconvenient to maintain them while the dataframe is being updated.

Virtual columns is a tool to access such computed values of arbitrary type as if they are stored in the dataframe,
without physically storing them. For their construction, a callback function taking a row index and returning a scalar
value or an object is passed. The platform would call this function whenever a table cell is accessed
(for instance, when rendering a grid).

#### Virtual columns for scalar types

In the example below we are expanding a `demog` table with the two virtual columns:
* `idx`, which only depends on a row's index
* `BMI`, which depends on two values of the row

Property `.isVirtual` indicates if a column is virtual.

```js
let table = grok.data.demo.demog();
table.columns.addNewVirtual('idx',
  (i) => i + 1, DG.TYPE.INT);
table.columns.addNewVirtual('BMI',
  (i) => table.row(i).weight / Math.pow(table.row(i).height / 100, 2), DG.COLUMN_TYPE.FLOAT);
grok.shell.add(table);
grok.shell.info(table.columns['idx'].isVirtual); // shows 'true`
```

Access virtual columns in a most standard way:

```js
console.shell.info(table.get('idx', 11)); // displays '12'
```

<!-- TODO: A picture -->

#### Virtual columns for objects

In addition to scalar types, it is possible to store JavaScript objects in columns of type [Object]().
However, if such object is a proxy, domain-specific representation of row's data, most often it is not desired to
physically store such objects in the dataframe. A virtual column with a type of `DG.COLUMN_TYPE.OBJECT` (it is
a default value for a type in `.addNewVirtual`) with a callback function returning a new object is a convenient
alternative to an `OBJECT`-typed `Column`. Origin data still resides in the most efficient way in the
[compressed columnar format](#dataframe-design), at the same time developers have an option to access it using
custom, domain-specific objects which get constructed on the fly.

In the example below we expand the `demog` dataframe with objects of a class `PersonalData`:

```javascript
class PersonalData {
  constructor(age, state = '', SSN = '') {
    this.age = age;
    this.state = state;
    this.SSN = SSN;
  }
  toString() {
    return `Social data: Age ${this.age}, ${this.state} ${this.SSN}`;
  }
}

let table = grok.data.demo.demog();
table.columns.addNewVirtual('personalData', (i) => new PersonalData(table.row(i).age, 'New York'));
grok.shell.info(table.row(1).personalData.state); // will emit 'New York'
grok.shell.add(table);
```

The platform uses `.toString` method to render the virtual column contents in the [grid][].

As such objects are virtual, a direct modification of their fields won't take any effect on the column's data,
unless the object class is implemented in such way that it modifies the original dataframe (for example, in a
JavaScript property setter).

### Filter

### Syncing

### Custom value comparers

## `.tags` and `.temp` collections

[097]: https://github.com/datagrok-ai/public/blob/e444332f49d6672c0fc93ad54dcbc601d3f332a1/js-api/src/dataframe.ts#L1431 "Class Cell"
[098]: https://github.com/datagrok-ai/public/blob/e444332f49d6672c0fc93ad54dcbc601d3f332a1/js-api/src/dataframe.ts#L1330 "Class RowList"
[099]: https://github.com/datagrok-ai/public/blob/e444332f49d6672c0fc93ad54dcbc601d3f332a1/js-api/src/dataframe.ts#L136 "Class ColumnList"
[100]: #dg-datagrame "Dataframe"
[101]: #dg-column "Dataframe column"
[102]: visualize/viewers.md "Datagrok Viewers"
[103]: https://github.com/datagrok-ai/public/blob/c4b913ef931e457144f773b1d8c55430c509657e/js-api/src/const.ts#L50 "DG.COLUMN_TYPE"
[104]: https://github.com/datagrok-ai/public/blob/c4b913ef931e457144f773b1d8c55430c509657e/js-api/src/const.ts#L39 "NULL constants"
[105]: #column-data-types "Column data types"
[106]: #accessing-and-modifying-column-values "Accessing and modifying column values"
[107]: overview/table-view.md "Table Views"
[108]: discover/tags.md#format "Values of a format tag"
[109]: discover/semantic-types.md "Semantic types"
[110]: ../domains/chem/cheminformatics.md "Cheminformatics overview"
[111]: https://github.com/datagrok-ai/public/blob/1393df83ef2eea80dc2c19b5e6e541cb35d60f91/packages/Chem/detectors.js#L6 "Chem Molecule detector"
[112]: https://github.com/datagrok-ai/public/tree/master/packages/Chem "Chem package"
[113]: develop/how-to/define-semantic-type-detectors.md "Defining semantic types detectors"
[114]: #tags-and-temp-collections "Tags and Temp collections"
[115]: visualize/viewers/grid.md "Grid"
[116]: #null-values-for-numeric-types "NULL values for numeric types"
[117]: #constructing-a-column
[118]: https://day.js.org/
[119]: visualize/viewers/grid.md
[120]: link
[121]: https://public.datagrok.ai/files/demo.testjobs.files.demofiles/ "Demo files"
[122]: #construct-from-in-place-content "Construct a dataframe from in-place content"
[123]: https://public.datagrok.ai/js/samples/data-frame/find-columns "Find columns"
[124]: ../how-to/customize-grid.md#color-coding "Color coding"
[125]: #numerical-and-categotical-columns "Numerical and categorical columns"
[126]: https://public.datagrok.ai/js/samples/data-frame/modification/manipulate "Manipulating dataframes"
[127]: https://public.datagrok.ai/js/samples/data-frame/advanced/custom-category-order "Custom categories order"
[128]: #versioning
[129]: https://public.datagrok.ai/js/samples/data-frame/stats
[130]: https://github.com/datagrok-ai/public/blob/be5486c9c761947bd916202343fad0adb2deaef3/js-api/src/dataframe.ts#L1669
[131]: https://public.datagrok.ai/js/samples/data-frame/modification/calculated-columns
[132]: overview/functions.md
[133]: https://public.datagrok.ai/js/samples/data-frame/modification/add-columns
[134]: https://dev.datagrok.ai/js/samples/data-frame/modification/manipulate
[135]: #virtual-columns
[136]: #work-with-categories