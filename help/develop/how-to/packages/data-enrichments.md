---
title: "Publishing enrichments with Datagrok plugins"
---

## Overview

Datagrok plugins can include [enrichments](../../../access/databases/databases.md#data-enrichment) that let users join additional database fields to a table with one click.

Enrichments are defined as JSON files and can be bundled with plugins.

> Note: PowerPack plugin should be installed to use the Enrichments feature.

## Where to put enrichments

Place enrichment JSON files in the plugin’s `enrichments/` folder:

```
my-plugin/
└─ enrichments/
   └─ my-enrichment.json
```

Datagrok loads these files automatically when the plugin is published.

## Enrichment JSON Format

Each file defines one enrichment:

```json
{
  "name": "Company info",
  "connection": "Dbtests:Northwind",
  "keySchema": "public",
  "keyTable": "orders",
  "keyColumn": "customerid",
  "fields": [
    "public.customers.companyname",
    "public.customers.contacttitle"
  ],
  "joins": [
    {
      "leftTableName": "public.orders",
      "rightTableName": "public.customers",
      "joinType": "left",
      "leftTableKeys": ["customerid"],
      "rightTableKeys": ["customerid"]
    }
  ]
}
```

### Fields

| Field        | Description                                 |
|--------------|---------------------------------------------|
| `name`       | Display name shown to users                 |
| `connection` | Database connection nqName                  |
| `keySchema`  | Schema of the key table                     |
| `keyTable`   | Key table                                   |
| `keyColumn`  | Key column                                  |
| `fields`     | Columns to retrieve (`schema.table.column`) |
| `joins`      | Table join definitions                      |

### Joins

Each join specifies how tables are connected:

* `leftTableName`, `rightTableName`: fully qualified table names
* `joinType`: left, inner, right
* `leftTableKeys`, `rightTableKeys`: join columns

## Conclusion

After publishing, the enrichment appears in the column context panel and can be applied with one click.


See also:
- [Enrichments](../../../access/databases/databases.md#data-enrichment)
- [Packages](../../develop.md#packages)
- [Creating a package](./create-package.md)
