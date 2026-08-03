---
title: "Git"
description: Connect Datagrok to a Git repository as a file share, including private repositories.
keywords:
  - git repository
  - github
  - version control
  - private repo
  - source code repository
  - file indexing
---

Provides access to [Git](https://git-scm.com/) repository as [file share](files.md). Works with any Git repository.

## Connection parameters

```json
{
  "url": "https://github.com/datagrok-ai/public.git",
  "dir": "packages/ApiSamples",
  "indexFiles": true,
  "privateRepo": false
}
```

See also:

* [Data connection](../../access.md#data-connection)
* [Git](https://git-scm.com/)
* [Files](files.md)
