---
title: "Parameter annotation"
---

There are various types of [functions](functions.md) such as [scripts](../../../compute/scripting.md) or
[queries](../../../access/access.md#data-query). What is common to all of them is the annotation of parameters. This is part of
the mechanism that enables universal support for functions in the platform.

A function annotation, also known as a header, is a multi-line comment that is placed above the function declaration.
It contains information about the function name, role, as well as function inputs and outputs. Inputs and outputs
can have metadata associated with them as well.

## Function

There are general parameters common to all functions, as well as parameters specific to certain function types. Not all
general parameters are required, the list of parameters depends on the function
[type](#parameters-specific-to-function-type), [role](../../../develop/function-roles.md), and so on.

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

### Supported types

Datagrok supports the following types in all scripting languages: 

* Scalars: `int`, `double`, `bool`, `string`, `datetime`
* Table: `dataframe`, `column`, `column_list`
* Collections: `list` (typically of strings)
* `graphics`: typically a function output. See [example]()
* `file`: when the script is executed, contains a string with the path to a file
* `blob`: array of bytes


### Options

#### Common

| Option     | Value                     | Description              |
|------------|---------------------------|--------------------------|
| validators | List separated with comma | List of named validators |
| caption    | Text string               | Custom field caption     |
| postfix    | Text string               | Field postfix            |
| units      | Same as postfix           | Value unit name          |

Named validators:

* containsMissingValues
* columnName - checks if table contains the column name
* columnIsNumerical
* columnIsCategorical.

Any function can be a validator if it returns `null` when validation is right, or string with error message otherwise.

#### For "dataframe" type

| Option      | Value       | Description                             |
|-------------|-------------|-----------------------------------------|
| columns     | numerical   | Only numerical columns will be loaded   |
| columns     | categorical | Only categorical columns will be loaded |

#### For "column" and "column_list" types

| Option     | Value                           | Description                                                                 |
|------------|---------------------------------|-----------------------------------------------------------------------------|
| type       | numerical                       | In a dialog, only numerical columns will be shown                           |
| type       | categorical                     | In a dialog, only categorical columns will be shown                         |
| type       | dateTime                        | In a dialog, only dateTime columns will be shown                            |
| format     | MM/dd/yyyy                      | Datetime format, for dateTime columns and datetime type only                |
| allowNulls | true/false                      | Adds validation of the missing values presence                              |
| action     | join("table parameter name")    | Joins result to the specified table, for output parameters only             |
| action     | replace("table parameter name") | Replaces result with columns in specified table, for output parameters only |

#### For "string" type

| Option      | Value                                                                                      | Description                              |
|-------------|--------------------------------------------------------------------------------------------|------------------------------------------|
| choices     | A comma-separated list of values, or a function name that returns a list of strings        | List of choices for string parameter     |
| suggestions | Name of a function that returns a list of strings to be used as suggestions as the user types the value | List of suggestions for string parameter |

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
