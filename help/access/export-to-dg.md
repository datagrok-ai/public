<!-- TITLE: Export Data to Datagrok -->
<!-- SUBTITLE: -->

# Export Data to Datagrok

Say, you have a system that works with tabular data, and you want to implement "Open in Datagrok" feature.
Datagrok has an endpoint to store file in user's storage. All you need to do is to make single POST request with CSV data.

```http request
method: POST
uri: $DATAGROK_URI/api/files/$LOGIN/uploads/any/path/you/want/file.csv

CSV_DATA
```
- set $LOGIN to User login to own the file.
- note, that you can pick any/path/you/want to keep file in user's upload folder.
- after you perform the request, Datagrok replies with URL you need to open in browser.

## Layout

If you want to apply a layout automatically:

- Save the layout [manually](../overview/table-view.md) or [programmatically](../develop/how-to/layouts.md)
- Go to Manage - Layouts, find your saved layout. Navigate to Property Panel, hit `Links...`
- In the new window, grab ID or Grok name of the layout. If layout name is too generic - just rename it and try again
- Pass ID or Name as `?layout` parameter to the request.