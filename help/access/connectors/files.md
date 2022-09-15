<!-- TITLE: Files -->
<!-- SUBTITLE: -->

# Files

File shares are often used for storing and exchanging data, and Datagrok provides first-class support for them. Once a
file share is mounted as a network drive on a server and registered within the platform, its content gets automatically
[indexed](../../access/files-indexer.md) and can be browsed.

## Connection parameters

```json
{
  "parameters": {
    "dir": "",
    "indexFiles": true
  },
  "credentials": {
    "parameters": {
      "login": "",
      "password": ""
    }
  }
}
```

See also:

* [File Manager](link)
* [Data connection](../data-connection.md)
* [File Indexer](../../access/files-indexer.md)
