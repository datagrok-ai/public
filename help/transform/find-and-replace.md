<!-- TITLE: Find and replace -->
<!-- SUBTITLE: -->

# Find and replace

This is the typical "Find and Replace" dialog that you see in every text editor, except that it supports specifying the
columns to run against, and [search patterns](../explore/data-search-patterns.md)
for matching non-textual columns. It means that you can search for 'this week' in the datetime column, and replace it
with a value that the datetime column understands, such as 'Nov 7, 2000' or '
7/11/2000'.

`Match case`, `Match whole word`, and `Use regular expressions` fields affect only string columns.

When a replace command is executed, it is logged into [console](../overview/navigation.md#console).

See also

* [Search](../explore/data-search.md)
* [Search patterns](../explore/data-search-patterns.md)
