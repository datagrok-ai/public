<!-- TITLE: Link tables -->
<!-- SUBTITLE: -->

# Link tables

Links two tables based on the key columns. The link type determines which records should be synchronized between the
datasets (current record, filter, selection, and mouse-over record), which gives us the following combinations:

* current row to row
* current row to selection
* current row to filter
* mouse-over row to selection
* mouse-over row to filter
* filter to filter
* filter to selection
* selection to filter
* selection to selection

For example, "selection to filter" means that all records of the second table, except for those that correspond to the
selected rows in the main table, will be filtered out. This is due to the fact that selection in the first table
controls the filter in the second one.

See also:

* [JavaScript API Samples: Linking Tables](https://public.datagrok.ai/js/samples/data-frame/link-tables)
* [Join tables](../transform/join-tables.md)
