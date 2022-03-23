<!-- TITLE: Parameterized queries -->

# Parameterized queries

A parameterized query is a query with one or more parameters. When the query is executed from the UI, a user is prompted
to enter parameters. It is also possible to run a query programmatically with the specified parameters
(see this [code snippet](https://public.datagrok.ai/js/samples/data-access/parameterized-query)). See how it works
in [this video](https://www.youtube.com/watch?v=sSJp5CXcYKQ&ab_channel=Datagrok).

## Creating a parameterized query

To create a parameterized query, open `Data | Databases`, right-click on a connection, and select `Add query...`. After
that, annotate parameters in the query header in using SQL/Sparql comments: `--` for SQL, and `#` for Sparql.

Here is an example of a simple SQL query that defines a `productName` input parameter:

```sql
--input: string productName
select * from products where productname = @productName
```

## Syntax

The syntax for defining query parameters is based on [scripting](../compute/scripting.md) with additions specific to
queries.

### Query parameters

| Parameter      | Description            |
|----------------|------------------------|
| `name`         | Name                   |
| `friendlyName` | Friendly name          |
| `description`  | Description            |
| `help`         | Help URL               |
| `tags`         | Tags                   |
| `input`        | An input parameter     |
| `output`       | An output parameter    |

All parameters are optional.

### Input parameters

#### `input` format template

The following format is used for query input parameters:

```sql
--input: <type> <name> = <value> {<option>: <value>; ...} [<description>]
```

#### Supported types

Is one of:

* **`int`** – integer scalar
* **`double`** – float scalar
* **`bool`** – boolean scalar
* **`string`** – string
* **`datetime`** – DateTime

<!-- TODO -->

* **`list<T>`** — a list of type `T` (currently `string` is supported)

##### Using lists in queries

To use lists as inputs to the queries, consider the following transformation:

* Specify `list<T>` as a type of a parameter `p`, `T` is the type of list elements
* Wrap the use of `p` inside the SQL query into `= ANY` operator: `= ANY(p)`, or a similar operator with an alternative
  selection of the comparison type, such as `>= ANY` or `< ANY`

For example, here is how we can transform a query to `northwind` taking a single `string`
parameter for a country:

```
--input: string country
select * from customers where country = @country
```

into a query taking a comma-separated list of countries:

```
--input: list<string> countries
select * from customers where country = ANY(@countries)
```

Learn more about using the lists feature in this video: [link](https://www.youtube.com/watch?v=meRAEF7ogtw).

#### Choices and suggestions

Options for supported types are described in the [Scripting](../compute/scripting.md) section.

| Option        | Description                                                                       |
|---------------|-----------------------------------------------------------------------------------|
| `choices`     | A comma-separated list of values, a name of the query, or the actual SQL query    |
| `suggestions` | Name of the query to be called to generate suggestion as the user types the value |

Examples:

```sql
--input: string shipCountry = France {choices: ['France', 'Italy', 'Germany']}
--input: string shipCountry = France {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCountry = France {choices: Demo:northwind:countries}
--input: string shipCountry = France {suggestions: Demo:northwind:countries}
```

#### Re-using input parameters

It's possible to re-use one or more existing input parameters as values inside parameters' `choices`
queries:

```sql
--input: string firstLetter = F
--input: string shipCountry = France {choices: Query("SELECT DISTINCT shipCountry FROM Orders WHERE shipCountry LIKE @firstLetter || '%')}
SELECT * FROM Orders WHERE (shipCountry = @shipCountry)
```

This is handy for queries with hierarchical choices, where each following parameter is dependent on the previous.

<!--
This query can be used as a "suggestion" query. It accepts exactly one parameter,
which is what a user has typed in the input box so far:

```sql
--name: country
--input: string sub
SELECT DISTINCT shipCountry FROM Orders WHERE shipCountry LIKE '%' || @sub || '%'
```
-->

#### Patterns

Sometimes users need to enter the filtering criteria as free text. On a server side, this query would be parsed and
safely transformed to a proper SQL clause. Check out [search patterns](../explore/data-search-patterns.md)
for more details.

To use the feature, the input type has to be `string`, since the user will be entering a free-text query, and the actual
data type should be put in the `pattern` option. In the dependent query, a reference
`@<patternName>(columnName)` should be used to specify a pattern that will be evaluated against the specified column. In
the example below, `@freightValue(freight)` will be transformed into `freight > 200.0` if a value for `freightValue`
is specified as `> 200.0`:

```sql
--input: string freightValue = >= 10.0 {pattern: double}
select * from Orders where @freightValue(freight)
```

A `datetime` type is supported as well:

```sql
--input: string orderDate = after 1/1/1995 {pattern: datetime}
select * from orders where @orderDate(orderDate)
```

#### Patterns summary

| Type               | Value         | Description or Example       |
|--------------------|---------------|------------------------------|
| `num, int, double` | `=`           | `= 100`                      |
|                    | `>`           | `> 1.02`                     |
|                    | `>=`          | `>= 4.1`                     |
|                    | `<`           | `< 5`                        |
|                    | `<=`          | `<= 2`                       |
|                    | `in`          | `in (1, 3, 10.2)`            |
|                    | `min-max`     | `Range: 1.5-10.0`            |
| `string`           | `contains`    | `contains ea`                |
|                    | `starts with` | `starts with R`              |
|                    | `ends with`   | `ends with w`                |
|                    | `regex`       | `regex 1(\w+)1`              |
|                    | `in`          | `in (ab, "c d", "e\\"f\\"")` |
| `datetime`         | `anytime`     |                              |
|                    | `today`       |                              |
|                    | `this week`   |                              |
|                    | `this month`  |                              |
|                    | `this year`   |                              |
|                    | `yesterday`   |                              |
|                    | `last week`   |                              |
|                    | `last month`  |                              |
|                    | `last year`   |                              |
|                    | `before`      | `before July 1984`           |
|                    | `after`       | `after March 2001`           |
|                    | `min-max`     | `Range: 1941-1945`           |

### Output parameters

A dataframe is returned by default as a query result. If you plan to obtain a value of a different data type
(for instance, in your JavaScript code), you can explicitly specify it in the output parameter. Below is an example
from [Chembl](https://github.com/datagrok-ai/public/tree/master/packages/Chembl)
package:

```sql
--output: string smiles {semType: Molecule}
```

The query with such header parameter will output a string of the semantic type `Molecule`.

## Example

Here is an example of the parameterized query applicable to the "Northwind" database:

```sql
--input: int employeeId = 5
--input: string shipVia = = 3 {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = France {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCity = starts with r {pattern: string}
--input: bool freightLess1000 = true
--input: datetime requiredDate = 1/1/1995
--input: string orderDate = after 1/1/1995 {pattern: datetime}
SELECT * FROM Orders WHERE (employeeId = @employeeId)
    AND (freight >= @freight)
    AND @shipVia(shipVia)
    AND ((freight < 1000) OR NOT @freightLess1000)
    AND (shipCountry = @shipCountry)
    AND @shipCity(shipCity)
    AND @orderDate(orderDate)
    AND (requiredDate >= @requiredDate)
```

Run it on Datagrok: [link](https://public.datagrok.ai/func/Demo.TestJobs.PostgreSQL.Orders).

When this query is run, the following dialog with the auto-generated inputs appears. Note that inputs marked as patterns
allow users to enter expressions like "> 5" for numbers,
"after 2019" for dates, and "starts with r" for strings.

![Parametrized queries](parameterized-queries.png)

Behind the scenes, Datagrok will parse the free-text query and execute a parameterized, safe, provider-specific SQL
query on the backend.

## Videos

[![Parameterized queries](../uploads/youtube/data_access.png "Open on Youtube")](https://www.youtube.com/watch?v=dKrCk38A1m8&t=1980s)

[YouTube: Datagrok database parameterized queries](https://www.youtube.com/watch?v=sSJp5CXcYKQ&ab_channel=Datagrok)

See also:

* [Data query](data-query.md)
* [Search patterns](../explore/data-search-patterns.md)
* [Function](../overview/functions/function.md)
* [Scripting](../compute/scripting.md)
* [JavaScript API Samples](https://public.datagrok.ai/js/samples/data-access/parameterized-query)
