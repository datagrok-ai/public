<!-- TITLE: Parameterized Queries -->
<!-- SUBTITLE: -->

# Parameterized Queries

A parameterized query is a query in which placeholders are used for 
parameters, and the actual parameter values are supplied at execution time. 

## Creating a parameterized query

To create a parameterized query, open **File | Connect to Data**, right-click
on the connection, and select **Add Query**. After that, annotate parameter 
information in the query header. Below is an example of a simple query that
defines a "productName" parameter:

```$sql
--input: string productName
select * from products where name == '@productName'
```

## Syntax

Syntax for query parameters is based on [Scripting](../../features/scripting.md) 
with some specifics:

### Header parameters

| Parameter   | Description                        |
|-------------|------------------------------------|
| name        | Name                               |
| description | Description                        |
| help        | Help URL                           |
| tags        | Tags                               |
| input       | Input parameter                    |

### Format template for 'input':

```
--input: <type> <name> = <value> {<option tag>:<value>; ...} [<description>]  
```

**type** - parameter type:

*   **int** \- integer scalar
*   **double** \- float scalar
*   **bool** \- boolean scalar
*   **string** \- string
*   **datetime** \- datetime

Comments style can be used '#' for Sparql.

### Options

Options for supported data typed are described in [Scripting](../../features/scripting.md) section. 

| Option      | Value                     | Description                             |
|-------------|---------------------------|-----------------------------------------|
| choices     | List separates with comma | List of choices for string parameter    |
| suggestions | Query name                | Name of query tat generates suggestions |     

"choices" option can be result of defined query: as direct query definition or reference by name to 
existing query.

Examples:
```
--input: string shipCountry = France {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCountry = France {choices: northwind:countries}
--input: string shipCountry = France {suggestions: northwind:countries}
```

Example of query for suggestions generation (one parameter is required):
```
--name: country
--input: string sub
SELECT DISTINCT shipCountry FROM Orders WHERE shipCountry LIKE '%' || @sub || '%'
```

#### Patterns

For all parameters types is supported [Search Patterns](../../features/data-search-patterns.md) feature.

| Option  | Value               | Description                                                           |
|---------|---------------------|-----------------------------------------------------------------------|
| pattern | True parameter type | See in [Search Patterns](../../features/data-search-patterns.md) section |

To enable this feature set input parameter type to "string" and add "pattern" option 
with true column data type.

### SQL syntax

#### Parameters

```
@<parameter name>
```

#### Parameters with pattern option

```
@<parameter name>(<column name>)
```

## Example

Following example can be applied to "Northwind" database "Orders" table. 

```$sql
--input: int employeeId = 5
--input: string shipVia = = 3 {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = France {choices: ["France", "Germany", "USA", "Finland"]}
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

See also:

  * [Data Query](../data-query.md)
  * [Search Patterns](../../features/data-search-patterns.md)
  * [Function](../function.md)
  * [Scripting](../../features/scripting.md)
