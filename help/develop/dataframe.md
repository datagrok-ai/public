<!-- TITLE: Dataframe -->

# Dataframe

Dataframe is a tabular structure with strongly-typed columns of different types. A dataframe class `DG.DataFrame`,
which is used in virtually any Datagrok extension or application. It operates through a columnar in-memory data engine
which we implemented from scratch to support highly efficient operation with data in a modern browser along
with fast in-browser data visualizations. 

## Dataframe JavaScript API

The dataframe's class and its related classes are available at a usual [`DG`]() namespace.
Use classes [`DG.DataFrame`](/js-api/DataFrame.html) and [`DG.Column`](/js-api/Column.html) for data initialization
and access, along with instances of [`DG.ColumnList`](/js-api/ColumnList.html), [`DG.Row`](/js-api/Row.html),
[`DG.Cell`](/js-api/Cell.html) as related properties or return values, for data manipulation and data-related
event handling.

Since a dataframe stores data as list of columns, the functionality related to
constructing, modifying and efficiently accessing data is embodied in `DG.Column` class, whereas event handling,
visual aspects of working with dataframes, handy construction methods and row-based access is provided in `DG.DataFrame`.

## Dataframe design

A Datagrok dataframe reminds of functionally similar structures in
[Python](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html) or
[R](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/data.frame).

Compared to these, Datagrok implementation is considerably optimized in numerous ways:

* Made to work entirely in the browser
* Operates on raw data instead of the JavaScript wrapping objects
* Doesn't have blocking structures
* Works well in multithreaded environments
* Designed to allocate little space in memory
* Same dataframe code also runs on Datagrok servers (e.g. for serialization)
* While it isn't possible to control the browser entirely, the in-memory engine is designed in such way that all
  operations are performed efficiently on modern machines, involving
  * manual memory management with the column-based layout
  * adaptive bit storage
  * custom serialization codecs

Datagrok [visualizations]() are fast because they work directly
with the dataframe data instead of numerous layers of physical abstractions.

## `DG.Column`

### Column types

Dataframe's columns support the types specified in a [`DG.COLUMN_TYPE`]() enum, each described below.

#### Special values

In popular programming languages like JavaScript, Python and R, special values representing missing,
unset values, or values for the undefined results (such as not valid numbers), are called `None`, `NaN` (not a number),
`null`, `undefined`, and so forth. Datagrok's dataframe type system implements special handling of such values,
which is good for performance yet still allows for these kinds of values to be set.

##### `NaN` values

[`NaN`](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/NaN)
in JavaScript represents a value which is not a valid number. If a value of `NaN` is passed as a value
at the 

##### `null` values

 For performance reasons,
Datagrok type system doesn't support JavaScript's `NaN` as a value of
a column item. Instead, it introduces values reserved for a `None` (null) for


#### String

<!-- In particular, describe a behavior with nulls -->

#### Integer

#### Float

#### Boolean

#### Datetime

#### BigInt

#### Qualified number

#### Object

<!-- Add the type to the enum -->

Enum value: `DG.COLUMN_TYPE.BIG_INT`

The type for working with integers that do not fit into 53 bits.

<!-- TODO: Figure out how to construct -->

### Column properties

Every column has:

* a `.name`: set it at construction, change at runtime 
* a .`type`: one of `DG.COLUMN_TYPE`
* a `.length`: it's a read-only property
* a parent `.dataFrame`, if part of a dataframe
* a semantic type — `.semType`: a string value
* `.temp` and `.tags` provide access to 

### Constructing a column

* The most common way is to use a method `.fromList`, explicitly specifying a type and a column name:  
  ```let col =  DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Column Name', [1, 2, 3])```
* The method `.fromStrings` recognizes and automatically assigns a type:  
  ```let col =  DG.Column.fromStrings('Column Name', ['3.14', '2.71']); // col.type === DG.COLUMN_TYPES.FLOAT```
* To create a column with NULL-values of a pre-specified length, use `.fromType` or concrete types shortcuts:  
  ```let col =  DG.Column.fromType(DG.COLUMN_TYPE.INT, 'Name', 3); // col.get(0) === DG.INT_NULL```  
  ```let col =  DG.Column.string('Name', 5); // col.get(2) === ""```
  
The column, once constructed, may later be [added to a dataframe]().

### Accessing and modifying column values

#### Accessing and modifying items

* A method `.get` is passed an index `i` to return an `i`-th value: `const value = column.get(idx);`  
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

<!-- TODO: Explain `notify` -->

#### Initializing with a function

#### Accessing raw data

If the fastest access is required for numerical columns, which usually happens in computing new values atop a column,
accessing data with a result of calling `.getRawData()` is advised:

```javascript
const table = grok.data.demo.demog(100000);
const column = table.columns.byName('age');
const array = column.getRawData();
const rowCount = column.length;
let sum = 0;
for (let i = 0; i < rowCount; i++)
  sum += array[i];
```

The `.getRawData` returns a `Float32Array` for `DG.COLUMN_TYPE.FLOAT`, `Int32Array` for `DG.COLUMN_TYPE.INT`.

To see the typical times it takes to run various column access patterns, run
[this example]((https://public.datagrok.ai/js/samples/data-frame/performance/access))
containing the methods from above.

## `DG.DataFrame`

### Constructing a dataframe

Dataframes may be obtained through the JavaSript or TypeScript code in various ways:

* a new dataframe constructed from a list columns
* a table already being rendered by a table view
* a dataframe constructed from a file in a file share
* a CSV file uploaded to a browser
* a dataframe returned by a script
* as calculated on the flight for aggregations

#### Constructing columns

#### Construct from columns

Dataframe's data resides in a list of columns. The fastest way to construct a new dataframe in-place is to
construct a list of columns in-place (specifying their types and names also and pass it to the `DataFrame`
method `fromColumns`):

```javascript
let newDataframe = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Sequence of integers', [1, 2, 3]),
  DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Floats below one', [0.1, 0.2, 0.3]),
  DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'String sequences', ["ABC", "CDE", "EFG"])
]);
```

<!-- TODO: adapt the samples accordingly — with a better naming -->

### Manipulating dataframes

#### Extending

#### Joining

#### Aggregating

#### Pivoting

### Versioning dataframes

### Syncing dataframes

### Virtual columns

Consider a person's `weight`, `height` and a derived value of a mass index `BMI`, which is computed as
`weight / height^2`. It's handy to access a `BMI` value as if it is stored in the dataframe along with
other column values. However, if access to such computed value only happens at some rows in some occasions,
it would not be practical to pre-compute them and actually store in the dataframe. It would also be inconvenient
to maintain these computed values while the dataframe is being modified.

Virtual columns is a tool to have such computed values of arbitrary type available as if they are stored in
the dataframe without actually storing them. For their construction, a callback function taking a row index and
returning a scalar value or an object is passed. The platform would call this function whenever a table cell is accessed
(for instance, when rendering a grid).

#### Virtual columns for scalar types

In the example below we are expanding a `demog` table with the two virtual columns:
* `idx`, which only depends on a row's index
* `BMI`, which depends on two values of the row

```js
let table = grok.data.demo.demog();
table.columns.addNewVirtual('idx',
  (i) => i + 1, DG.TYPE.INT);
table.columns.addNewVirtual('BMI',
  (i) => table.row(i).weight / Math.pow(table.row(i).height / 100, 2), DG.COLUMN_TYPE.FLOAT);
grok.shell.add(table);
```

Access virtual columns in a standard way:

```js
console.shell.info(table.get('idx', 11)); // displays '12'
```

<!-- TODO: A picture -->
<!-- TODO: Accessing these values -->

#### Virtual columns for objects

In addition to scalar types, it is possible to store JavaScript objects in columns of type [Object]().
However, if such object is a proxy, domain-specific representation of row's data, most often it is not desired to
physically store such objects in the dataframe. A virtual column with a type of `DG.COLUMN_TYPE.OBJECT` (it is
a default value for a type in `.addNewVirtual`) with a callback function returning a new object is a convenient
alternative to an `Object`-typed `Column`. Origin data still resides in the most efficient way in the
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
table.columns.addNewVirtual('DATA', (i) => new PersonalData(table.row(i).age, 'New York'));
grok.shell.info(table.row(1).DATA.state); // will emit 'New York'
grok.shell.add(table);
```

The `.toString` method is used to render the column contents by the grid.

As such objects are virtual, a direct modification of their fields won't take any effect on the column's data,
unless the object class is implemented in such way that it modifies the original dataframe (for example, in a
JavaScript property setter).

### Custom value comparers

## `.tags` and `.temp`