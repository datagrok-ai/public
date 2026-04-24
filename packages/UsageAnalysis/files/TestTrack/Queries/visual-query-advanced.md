1. Create a new visual query (with any parameter)
1. Run it
1. Set the custom name (e.g. `test_visual_query`)
1. Save it
1. Share it
1. Add a post-process, e.g. on line 7, add: 
`grok.shell.info(result.rowCount);`
1. Add a layout (viewers, color coding, change format, row size)
1. Save it
1. Close all
1. Click the query to preview - verify its name, new layout and that post-process executes
1. Edit the query: 

   1. change some settings on the Query tab (e.g. add or remove Order by)
   1. change the layout
   1. save
1. Close all
1. Click the query to preview - verify its name, new layout and that post-process executes
1. Run the query
1. On the Toolbox, change the parameter value and refresh
1. Add some more viewers
1. Delete some rows and columns
1. Refresh with Enrich on - layout should not change, deleted rows and columns should be restored
1. Save the project
1. Close all
1. Open the saved project