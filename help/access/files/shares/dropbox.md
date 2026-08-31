---
title: "Dropbox"
description: Connect Datagrok to Dropbox storage as a file share.
keywords:
  - dropbox storage
  - dropbox file share
  - cloud storage connector
  - file sharing service
---

Provides access to [Dropbox](https://www.dropbox.com) storage as [file share](../files.md).

## Connection parameters

```json
{
  "parameters": {
    "clientId": "",
    "dir": ""
  },
  "credentials": {
    "parameters": {
      "password": ""
    }
  }
}
```

`clientId` is the app key of the [Dropbox app](https://www.dropbox.com/developers/apps) that
authorizes access. Register your own app and paste its app key here — Datagrok does not ship one.

See also:

* [File shares](../files.md)
* [Dropbox](https://www.dropbox.com)
