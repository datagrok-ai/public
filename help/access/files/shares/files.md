---
title: "Files"
description: Mount a network drive or server folder as a Datagrok file share with automatic indexing.
keywords:
  - network drive
  - network share
  - windows share
  - linux file share
  - mounted drive
  - file indexing
  - login and password
---

Once you mount a file share as a network drive on a server and register it within the platform, its content gets
automatically indexed and can be browsed.

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

* [File shares](../files.md)
* [Data connection](../../access.md#data-connection)
