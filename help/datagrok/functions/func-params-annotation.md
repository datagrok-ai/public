---
title: "Parameter annotation"
---

There are various types of [functions](functions.md) such as [scripts](../../compute/scripting.md) or
[queries](../../access/data-query.md). What is common to all of them is the annotation of parameters. This is part of
the mechanism that enables universal support for functions in the platform.

A function annotation is a multi-line comment that is placed above the function declaration and contains information
about the function name, role, input and output parameters, and so on. The parameter annotation is part of the function
annotation.

## Header parameters

There are general parameters common to all functions, as well as parameters specific to certain function types. Not all
general parameters are required, the list of parameters depends on the function
[type](#parameters-specific-to-function-type), [role](../../develop/function-roles.md), and so on.

### Common parameters

| Parameter   | Description                                                                         |
|-------------|-------------------------------------------------------------------------------------|
| name        | Name                                                                                |
| description | Description                                                                         |
| tags        | Tags (see also `DG.TAGS`)                                                           |
| input       | Input parameter                                                                     |
| output      | Output parameter                                                                    |
| help-url    | Datagrok's Wiki URL                                                                 |
| reference   | Reference to a research paper, Wikipedia article, Git repository, etc.              |
| top-menu    | Top menu path separated with pipes ("\|")                                           |

### Parameters specific to function type

| Type       | Parameter   | Description                                                                                             |
|------------|-------------|---------------------------------------------------------------------------------------------------------|
| Script     | language    | Script language (see the list of [supported languages](../../compute/scripting.md#supported-languages)) |
| Script     | sample      | Name of a sample file                                                                                   |
| Script     | environment | Environment name                                                                                        |
| Script, info panel | condition   | Function applicability conditions                                                               |
| Query      | connection  | Connection name                                                                                         |

:::tip Tip

You can also add custom parameters using the `meta.` prefix. For example, `meta.icon` or `meta.ext`.

:::

### Inputs and outputs

The syntax for defining input and output parameters looks like this:

```
<comment symbol><direction>: <type> <name> = <value> {<option tag>:<value>; ...} [<description>]
```

1. **direction** - parameter direction:
    1. input
    2. output
2. **type** - parameter type:
    * **int** - integer (scalar)
    * **double** - float (scalar)
    * **bool** - boolean (scalar)
    * **string** - string (scalar)
    * **dataframe** - dataframe
    * **column** - column from selected table
    * **column_list** - list of columns from selected table
    * **datetime** - datetime
    * **graphics** - graphics
    * **file** - file, variable in the function body contains path to file
    * **blob** - blob, variable in the function body contains path to binary file
    * **list<T>** - list of objects of type T (support depends on function type)
3. **name** - parameter name that will be used in the function body. Optional for a graphical output.
4. **value** - default value. For scalar inputs, it corresponds to a value; for graphics outputs, it is used as index of
   graphical element. Optional.
5. **options** - the list of options as key-value pairs. It is used to enhance the UI of a function call. Optional.
6. **description** - brief description for a parameter that appears in UI. Optional.

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
* columnIsCategorical

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
