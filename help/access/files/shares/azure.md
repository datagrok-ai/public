---
title: "Azure Blob"
---

Provides access to the
[Azure](https://azure.microsoft.com/en-us) as
[file share](../files.md).

The connector supports two authorization methods: [shared tokens](https://learn.microsoft.com/en-us/azure/storage/common/storage-account-keys-manage?tabs=azure-portal#view-account-access-keys) and [SAS](https://learn.microsoft.com/en-us/azure/ai-services/translator/document-translation/how-to-guides/create-sas-tokens?tabs=Containers#create-sas-tokens-in-the-azure-portal). Shared tokens work at the account level, while SAS works at the container level. SAS also has lifetimes and a set of permissions. Reading and listing are the minimal rights for read-only usage. You can also add writing and delete rights.

## Connection parameters

````json
{
  "parameters": {
    "account": "",
    "container": "",
    "use SAS": false, // select between "shared token" and SAS
  },
  "credentials": {
    "parameters": {
      "shared token": "",
      "SAS": "",
    }
  }
}
````

See also:

* [File shares](../files.md)
* [Azure](https://azure.microsoft.com/en-us)
