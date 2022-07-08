<!-- TITLE: Upload data -->
<!-- SUBTITLE: -->

# Dataframe

To upload a CSV file to the user storage, make the following POST request with the CSV content as body:

```http request
method: POST
uri: $DATAGROK_URI/api/files/$LOGIN/uploads/any/path/you/want/file.csv

CSV_DATA
```

* Set $LOGIN to the User login to own the file.
* Use any path after /uploads/.

The server replies to this request with the URL of the uploaded project.

## Layout

To apply a layout, use the `?layout` parameter in the following way:

1. Save the layout [manually](../../overview/table-view.md) or [programmatically](layouts.md)
2. Go to **Manage / Layouts**, find your saved layout. Navigate to Property Panel, hit `Links...`
3. In the new window, copy either ID or fully qualified name of the layout. You might want to rename the layout if the
   name is ambiguous.
4. Pass ID or Name as a `?layout` parameter to the request.

# Example

```curl
curl --location --request POST 'http://localhost:8080/api/files/alex.aprm/uploads/from_excel/test.csv?layout=alexaprm.superlayout' \
--header 'Content-Type: text/csv' \
--data-raw 'a,b
1,2'
```
