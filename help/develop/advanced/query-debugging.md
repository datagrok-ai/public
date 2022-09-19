<!-- TITLE: Query debugging -->
<!-- SUBTITLE: -->

# Query debugging

To debug and troubleshoot the database queries, use the built-in tool
that allows you to see the following information about the query:

- Column/Rows - the number of columns and rows in the requested table
- Query - the resulting query string with substituted parameters
- Params array - the array of query parameters
- Column - the information about column name, type, type name, precision, and scale
- Java type - the type which is used in the Java for this column
- Execution time by steps - the time benchmark for each computational step
- Execution time - the time benchmark for the all computational steps
- Memory - the used and free JVM memory

You can find this information by doing the following:

- Open the Datagrok inspector by pressing Alt+I (or Tools | Dev | Inspector)
- Enable 'Debug Queries' on the 'Debug' pane. You might close the inspector now.
- Execute the query you want to profile
- Open your requested table in the property panel by clicking on the table header
- Click the 'DataQuery' in the property panel
- Click the 'Last call'
- Expand the 'Log' menu item
