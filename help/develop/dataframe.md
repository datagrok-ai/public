<!-- TITLE: Dataframe -->

# Dataframe

Dataframe is a tabular structure with strongly-typed columns of different types. A dataframe class `DG.DataFrame`,
which would be used in virtually any Datagrok extension or application, sits on top of a columnar in-memory data engine
which we implemented from scratch to support highly efficient operation with data in a modern browser along
with fast in-browser data visualizations. 

## Dataframe JavaScript API

The dataframe's class and its related classes are available at a usual `DG` namespace.
Use classes [`DG.DataFrame`](/js-api/DataFrame.html) and [`DG.Column`](/js-api/Column.html), along with their
properties or return values being instances of [`DG.ColumnList`](/js-api/ColumnList.html), [`DG.Row`](/js-api/Row.html),
[`DG.Cell`](/js-api/Cell.html), for data manipulation and data-related event handling.

Since a dataframe stores data as list of columns, the functionality related to
constructing, modifying and efficiently accessing data is embodied in `DG.Column` class, whereas event handling,
visual aspects of working with dataframes, handy construction methods and row-based access is provided in `DG.DataFrame`.

We review the contents of the above classes in detail through this document.  

## Dataframe design

A Datagrok dataframe reminds of functionally similar structures in
[Python](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html) or
[R](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/data.frame).

However, Datagrok's implementation is much more efficient in numerous ways:

* Made to work entirely in the browser
* Same dataframe code also runs on Datagrok servers (e.g. for serialization)
* While it isn't possible to control the browser entirely, the in-memory engine is designed in such way that all
  operations are performed efficiently on modern machines, involving
  * manual memory management with the column-based layout
  * adaptive bit storage
  * custom serialization codecs
* Operates on raw data instead of the JavaScript wrapping objects
* Doesn't have blocking structures
* Works well in multithreaded environments
* Designed to allocate little space in memory

Our [visualizations]() are fast because they work directly
with the dataframe data instead of numerous layers of physical abstractions.

## `DG.Column`

### Column types

Dataframe's columns support the types specified in a [`DG.COLUMN_TYPE`]() enum, each described below.

#### Special values

In popular programming languages like JavaScript, Python and R, common values representing missing values,
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
* a `.length`: read-only property
* a parent `.dataFrame`
* a semantic type — `.semType` (string value)

### Constructing a column

* The most common way is to use a method `.fromList`, explicitly specifying a type and a column name:  
  ```let col =  DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Column Name', [1, 2, 3])```
* The method `.fromStrings` recognizes and automatically assigns a type:  
  ```let col =  DG.Column.fromStrings('Column Name', ['3.14', '2.71']); // col.type === DG.COLUMN_TYPES.FLOAT```
* To create a column with NULL-values of a pre-specified length, use `.fromType` or concrete types shortcuts:  
  ```let col =  DG.Column.fromType(DG.COLUMN_TYPE.INT, 'Name', 3); // col.get(0) === DG.INT_NULL```  
  ```let col =  DG.Column.string('Name', 5); // col.get(2) === ""```
  
The column, once constructed, may later be added to a dataframe.

## `DG.DataFrame`

### Constructing a dataframe

Dataframes may be obtained through the JavaSript or TypeScript code in various ways:

* a new dataframe constructed from a set columns
* a table already being rendered by a table view
* a dataframe constructed from a file in a file share
* a CSV file uploaded to a browser
* a dataframe returned by a script
* calculated on the flight for aggregations

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
`weight / height^2`. It would be convenient to access a `BMI` value as if it is stored in the dataframe along with
other values. However, if access to such computed value only happens to some rows in some occasions, it would not be
practical to pre-compute them and actually store in the dataframe. It would also be inconvenient to maintain
these computed values while the dataframe is being modified.

Virtual columns is a tool to have such computed values of arbitrary type available as if they are stored in
the dataframe without actually storing them. For their construction, a callback function taking a row index and
returning a scalar value or an object is passed. The platform would call this function whenever a table cell is accessed
(for instance, when rendering a grid).

#### Virtual columns for scalar types

In the example below we are expanding a `demog` table with the two virtual columns:
* `idx`, which only depends on a row's index
* `bmi`, which depends on two values of the row

```js
let table = grok.data.demo.demog();
table.columns.addNewVirtual('idx',
  (i) => i + 1, DG.TYPE.INT);
table.columns.addNewVirtual('BMI',
  (i) => table.row(i).weight / Math.pow(table.row(i).height / 100, 2), DG.COLUMN_TYPE.FLOAT);
grok.shell.add(table);
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