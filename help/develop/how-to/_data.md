<!-- TITLE: Efficiently work with data -->

# Efficiently work with data

As Datagrok is centralized around data, it's pivotal to know key means and best practices for working with data
programmatically on the platform. This document addresses various topics related to that.

Datagrok dataframes are highly optimized. They are implemented with our proprietary technology allowing to efficiently
work with huge datasets in the browser. Essentially, it is a columnar in-memory database that was engineered from
scratch and optimized for the purpose of exploratory data analysis, interactive visualizations, and machine learning.

Note that Datagrok dataframes live and operate entirely inside the browser, but not on our
[compute server](../admin/infrastructure.md#compute-components). However, it's possible to pass dataframes to scripts (in Python, R and others)
which run on the server, and get dataframes in return.

You get dataframes within your application in various ways. Dataframe may be a table rendered by a table view, a new
dataframe constructed from a set columns, a dataframe constructed from a file in a file share, a CSV file uploaded to a
browser, a dataframe returned by a script, and so forth. You can add and delete rows and columns of the dataframe, view
it in a table view and let viewer be attached to it to render. Dataframes can also be calculated on the flight for
aggregations.

Dataframes are comprised of [columns](). Columns may be used as Datagrok functions arguments. For columns, it's possible
to get its underlying dataframe. In return, columns are comprised of cells, and it's possible to get a cell's underlying
column. There is also a diverse [system of events](https://datagrok.ai/js-api/classes/dg.DataFrame) one can subscribe on
a dataframe.

### Data types

---

Use [DataFrame](https://datagrok.ai/js-api/classes/dg.DataFrame), [Column](https://datagrok.ai/js-api/classes/dg.Column)
, [ColumnList](https://datagrok.ai/js-api/classes/dg.ColumnList), and [Row](https://datagrok.ai/js-api/classes/dg.Row) classes for table manipulation.

```javascript
demog = grok.testData('demog', 5000);
demog.cols.remove('sex');
foo = demog.cols.addNew('foo', 'int');
demog.rows.removeAt(1, 3);
demog.rows.insertAt(2, 2);
demog.rows.addNew(['Spiderman', 'studyX', 'NYC', 32, 'Spider', 'Net', new Date(2020), 180, 80, 666]);
demog.rows.addNew().subj = 'Iron Man';

// alternative ways of setting values
foo.set(1, 777);
demog.set('age', 1, 44);

```

## Best practices

<!--
- use special data structures for high-performance tasks
- formatting in visualizations
- tags in dataframes, detectors
-->

### Formatting

### Tags

When annotating dataframes and columns with [metadata](../../discover/tags.md), take a moment to consider whether you
want a tag to be shown to the user. If metadata is attached to facilitate other parts of the data workflow, or you find
it difficult to read the value stored in a tag, it makes sense to discard it from tooltips where people can see it. To
do this, simply start a tag name with the `.` prefix.

### Detectors

Sematic type [detectors](define-semantic-type-detectors.md) are meant to be lightweight, simple, and efficient. Keep
this fact in mind when designing these functions. If columns you are interested in generally contain a lot of unique
values, use a special sampling method to run your checks on a random subset of column categories. Empty values deserve
special attention: make sure you don't match them with your semantic type as it will be confusing to find it assigned to
a column consisting solely of nulls.
