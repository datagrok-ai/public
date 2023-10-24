---
title: "Function annotations"
---

# Function annotations

There are various types of [functions](functions.md) such as [scripts](../../../compute/scripting.md) or
[queries](../../../access/access.md#data-query). What is common to all of them is the annotation of parameters. This is part of
the mechanism that enables universal support for functions in the platform.

A function annotation, also known as a header, is a multi-line comment that is placed above the function declaration.
It contains information about the function name, role, as well as function inputs and outputs. Inputs and outputs
can have metadata associated with them as well.

## Function

There are general parameters common to all functions, as well as parameters specific to certain function types. Not all
general parameters are required, the list of parameters depends on the function
[type](#parameter-types-and-options), [role](../../../develop/function-roles.md), and so on.

<details> <summary> Simple python script </summary> <div>

```python
#language: python
#name: GetCellNumber
#description: Calculates number of cells in the table
#tags: demo
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
count = table.shape[0] * table.shape[1]
```
</div> </details>

These are the common parameters for all functions:
* `name`: shows up in the user interface
* `description`: shows up in a function tooltip
* `tags`: comma-separated tags that you can use in search
* `help-url`: help that shows up when you click on the "?" in the function dialog
* `reference`: Reference to a research paper, Wikipedia article, Git repository, etc.
* `top-menu`: Top menu path separated with pipes (`|`), such as `Chem | Gasteiger Charges`

Some parameters are specific to script language and/or technology:
* Script
  * `language`: script language, so that Datagrok knows how to execute it
  * `environment`: [script environment](../../../compute/scripting.md#environments) (Conda environment for python, etc)
  * `sample`: path to a sample data csv file. When defined, a `*` icon appears on the ribbon panel that loads it. 
* Script-based info panels
  * `condition`: GrokScript condition that gets evaluated to decide whether to show the panel for the object
* Query
  * `connection`: Name of the db connection that the query should use

To add additional parameters, use the `meta.` prefix. They can be used for dynamically searching for
the functions of interest.


## Inputs and outputs

Each input and output definition take one line that starts with the comment, followed 
by type, name, optional default value, options, and optional description, just like here:

```
#input: string country {choices: ["USA", "Canada"]} [Country to ship from]
#input: double weight
#output: double price
```

### Parameter types and options

Datagrok supports the following types in all scripting languages: 

* Scalars: `int`, `double`, `bool`, `string`, `datetime`
* Table: `dataframe`, `column`, `column_list`
* Collections: `list` (typically of strings)
* `graphics`: typically a function output. See [example]()
* `file`: when the script is executed, contains a string with the path to a file
* `blob`: array of bytes

Some of the options apply to all parameters, while other are type-specific. 

For all parameters:

| Option      | Value  | Description                                        |
|-------------|--------|----------------------------------------------------|
| validators  | string | Comma-separated list of [validators](#validation) |
| caption     | string | Custom field caption                               |
| postfix     | string | Field postfix                                      |
| units       | string | Value unit name                                    |
| nullable    | bool   | Makes it an [optional parameter]()                 |

For `dataframe` type:

| Option      | Value       | Description                             |
|-------------|-------------|-----------------------------------------|
| columns     | numerical   | Only numerical columns will be loaded   |
| columns     | categorical | Only categorical columns will be loaded |

For `column` and `column_list` types

| Option     | Value                           | Description                                                                 |
|------------|---------------------------------|-----------------------------------------------------------------------------|
| type       | numerical,categorical,dateTime  | In a dialog, only numerical columns will be shown                           |
| format     | MM/dd/yyyy                      | Datetime format, for dateTime columns and datetime type only                |
| allowNulls | true/false                      | Adds validation of the missing values presence                              |
| action     | join("table parameter name")    | Joins result to the specified table, for output parameters only             |
| action     | replace("table parameter name") | Replaces result with columns in specified table, for output parameters only |

For `string` type

| Option                       | Value                                                                             | Description                                        |
|------------------------------|-----------------------------------------------------------------------------------|----------------------------------------------------|
| [choices](choices)           | Comma-separated list of values, or a function name that returns a list of strings | Makes it a combo box                               |
| [suggestions](#autocomplete) | Name of a function that returns a list of strings to autocomplete user input      | Autocomplete gives you options as you type |

For `numeric` types

| Option | Description                    |
|--------|--------------------------------|
| min    | Minimum value to be validated. |                                                                                |
| max    | Maximum value to be validated. | 

### Initial values and optional parameters

Proper handling of the empty parameters requires special efforts when building a SQL query,
and passing empty parameters to the function that does not expect it is a major source of errors. 
To deal with it, by default each parameter is required, but you can specify initial values
and make it optional (nullable).

* *Initial value* gets shown in the dialog when you execute the function. If you remove the value, an empty value
  is passed to the function. But before the function is executed, the input is validated, and you will get
  an error if a required parameter is not specified.
* *Optional parameter*: make a parameter optional by specifying the `nullable: true` option.

For example, to create an optional string parameter with the initial value:

```
--input: string shipCountry = "France" { nullable: true }
SELECT * FROM customers where shipCountry = @shipCountry
```

### Filter patterns

**Filter pattern** allows you to use free-text conditions like "this week" for dates, or ">50" for numbers.

To use search patterns, set the input type to `string`, and set the `pattern` option
to the type of the column you are filtering. Then, reference the patter in the query
as `@patternName(columnName)`, just like we did here for the "freight" column:

```
--input: string freightValue = >= 10.0 {pattern: double}
select * from Orders where @freightValue(freight)
```
Different inputs would produce differently structured SQL (also dependent on the database). 

| Input | SQL                                                                | Description                     |
|-------|--------------------------------------------------------------------|---------------------------------|
|       | select * from orders <br/> where 1 = 1                             | No input => no filter           |
| >3    | select * from orders <br/> where freight > 3                       | Using column name to filter     |
| 10-20 | select * from orders <br/> where (freight >= 10 and freight <= 20) | Have to do multiple comparisons |


In this example, the `freightValue` input parameter is defined as a string with a default value of `>= 10.0`. 
The `pattern` _option_ specifies that the actual data type is a `double`. In the query, a reference to
`@freightValue(freight)` specifies the _pattern_ that will be evaluated against the "freight" column.

Here's a list of all supported search patterns:

<details>
<summary> Patterns </summary>

| Type               | Value         | Description or example       |
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

To learn more, see [search patterns](../../../explore/search-filter-select/data-search-patterns.md).
</details>


### Choices

Use `choices` to make input a combo box, and restrict the selection to the defined set of values. 
You can define `choices` using a comma-separated list of values, a name of another query, 
or by writing an actual SQL query. Here's an example of how to define `choices` for a `shipCountry` 
input parameter using all three methods:

```sql
--input: string shipCountry = "France" {choices: ['France', 'Italy', 'Germany']}
--input: string shipCountry = "France" {choices: Demo:northwind:countries}
--input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
```

### Validation

Use validators to assure that the parameters are valid before calling the function. 

You can specify 
Named validators:

* containsMissingValues
* columnName - checks if table contains the column name
* columnIsNumerical
* columnIsCategorical.

Any function can be a validator if it returns `null` when validation is right, or string with error message otherwise.

### Autocomplete

Use the `suggestions` option to enable autocomplete, and specify the name of a function that 
accepts one string parameter, and returns a list of strings (or a dataframe with one string column)
to be used as suggestions as the user types the value. 

Here are two SQL functions 
(from the [Chembl](https://github.com/datagrok-ai/public/blob/master/packages/Chembl/queries/cartridge.sql#L63) package), 
where the UI for the "Structures by Organism" query uses the "organismsAutocomplete" function
to complete user input:

<details> <summary> Autocomplete function </summary> <div>

```sql
--name: organismsAutocomplete
--input: string sub
select distinct organism from target_dictionary
where organism ilike '%' || @sub || '%'
limit 50
```

```sql
--name: StructuresByOrganism
--input: string organism = "Shigella" {suggestions: Chembl:organismsAutocomplete}
SELECT md.chembl_id AS compound_chembl_id,
cs.canonical_smiles,
act.standard_type,
act.standard_value,
act.standard_units,
td.chembl_id AS target_chembl_id,
td.organism,   td.pref_name
FROM target_dictionary td
  JOIN assays a ON td.tid = a.tid
  JOIN activities act ON a.assay_id = act.assay_id
  JOIN molecule_dictionary md ON act.molregno = md.molregno
  JOIN compound_structures cs ON md.molregno   = cs.molregno
  JOIN organism_class oc ON td.tax_id = oc.tax_id
    AND td.organism = @organism
    AND oc.L1 = 'Bacteria'
```
</div></details>

![](autocomplete.gif)


## Examples

<details>
<summary> TypeScript function </summary>
<div>

```ts
//name: Len
//description: Calculates the length of a string
//input: string s
//output: int n
export function getLength(s: string): number {
  return s.length;
}
```

</div>
</details>

<details>
<summary> Python script </summary>
<div>

```python
#name: Template
#description: Calculates number of cells in the table
#language: python
#tags: template, demo
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
count = table.shape[0] * table.shape[1]
```

</div>
</details>

<details>
<summary> Query </summary>
<div>

```sql
--name: protein classification
--connection: chembl
select * from protein_classification;
--end
```

</div>
</details>

<details>
<summary> Complex annotation example </summary>
<div>

```python
#input: dataframe t1 {columns:numerical} [first input data table]
#input: dataframe t2 {columns:numerical} [second input data table]
#input: column x {type:numerical; table:t1} [x axis column name]
#input: column y {type:numerical} [y axis column name]
#input: column date {type:datetime; format:mm/dd/yyyy} [date column name]
#input: column_list numdata {type:numerical; table:t1} [numerical columns names]
#input: int numcomp = 2 {range:2-7} [number of components]
#input: bool center = true [number of components]
#input: string type = high {choices: ["high", "low"]} [type of filter]
#output: dataframe result {action:join(t1)} [pca components]
#output: graphics scatter [scatter plot]
```

</div>
</details>

See also:

* [Functions](functions.md)
* [Function parameters enhancement](func-params-enhancement.md)
