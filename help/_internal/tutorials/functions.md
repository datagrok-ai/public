<!-- TITLE: Use Cases: Functions -->
<!-- SUBTITLE: -->

# Use cases: Functions

Owner: Andrey

Goal: demo advanced capabilities for data processing

Features: functions, function browser, query editor, console, macros (recording and reusing), auto UI parameter builder,
help, auto-completion, invocation of UI commands/scripts,

## Functions overview

Everything that accepts something and returns something is a function in Datagrok. DataQueries accept parameters and return data, 
package functions also take arguments and optionally return data or another parameters.

You can execute functions via UI or [programmatically](https://public.datagrok.ai/js/samples/functions/calling), combine them in any possible pipelines.

Datagrok can build UI for any function automatically using parameters annotation:
```python
#name: BMI
#description: Body Mass Index
#language: python
#input: double height
#input: double weight
#output: double bmi

bmi = weight / height**2
```

## Console

Press `~` to show console. You can see, that each function call is reflected in the Datagrok console. You can copy and re-run every line.
Hit `tab` to auto-complete function.

## Queries

Go to Data -> Datasets and create a new query. Note, that query has a name, and it's located in your namespace. 
After you run it, you see in console something like this: 
```
> AlexAprm:TestQuery()
  TestQuery: TestQuery (1 rows, 1 columns)
```

You can reuse it for Grok-Script or in JS. See more in [joining data section](#joining-data)

## Files

Go to Data -> Files, and open any file. 

## Joining Data

You can reuse all functions in JS. Create script in Functions -> Scripts -> New JS Script
```javascript
let queryData = await grok.functions.call('AlexAprm:TestQuery' , {});
let filesData = await grok.functions.call('OpenFile' , {fullName: 'Demo:Files/some_data.csv'});
```
Note, `Demo:Files` is namespace-qualified file-connection name. Just like `AlexAprm:TestQuery`.
Then, you can join them
```javascript
//name: Join data
//language: javascript
//output: dataframe result
 
let queryData = await grok.functions.call('AlexAprm:TestQuery' , {});
let filesData = await grok.functions.call('OpenFile' , {fullName: 'Demo:Files/some_data.csv'});
let tj = grok.data.joinTables(queryData, filesData, 
  ['id1', 'id2'], 
  ['id3', 'id4'], 
  ['v1', 'v2'], 
  ['v3', 'v4'], DG.JOIN_TYPE.INNER, false);
```
See also:
- [Scripts](scripting.md)
- [Join Tables](https://public.datagrok.ai/js/samples/data-frame/join-link/join-tables.js)

You can also use Grok Script as a script language:
```
//name: Join data
//language: grok
//output: dataframe result
 
queryData = AlexAprm:TestQuery();
filesData = OpenFile('Demo:Files/some_data.csv');
result = JoinTables(queryData, filesData, ["ID"], ["ID"], ["ID", "D1", "D2"], ["ID", "F1", "F2"])
```