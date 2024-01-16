---
title: "Azure Blob"
---

Provides access to the
[Azure](https://azure.microsoft.com/en-us) as
[file share](../../files/files.mdx).

It supports two different authorization methods: [shared tokens](https://learn.microsoft.com/en-us/azure/storage/common/storage-account-keys-manage?tabs=azure-portal#view-account-access-keys) and [SAS](https://learn.microsoft.com/en-us/azure/ai-services/translator/document-translation/how-to-guides/create-sas-tokens?tabs=Containers#create-sas-tokens-in-the-azure-portal). Shared tokens are used on account level, SAS works on container level. SAS also has lifetimes and set of permissions. Reading and listing is the minimal set of rights for RO-usage (for RW, it is possible toadd writing and deleting rights as well).

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

* [File shares](../files.mdx)
* [Azure](https://azure.microsoft.com/en-us)
