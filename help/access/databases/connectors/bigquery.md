---
title: "BigQuery"
---

Provides access to [Google BigQuery](https://cloud.google.com/bigquery/) database using SQL queries via a JDBC driver.

## Connection parameters

```json
{
  "parameters": {
    "connString": "",
    "projectId": ""
  },
  "credentials": {
    "parameters": {
      "login": "",
      "password": ""
    }
  }
}
```

## Authentication

BigQuery connections support two authentication methods:

| Method              | Use case                                                                      |
|---------------------|-------------------------------------------------------------------------------|
| **Service Account** | Shared identity via a Google service-account JSON key                         |
| **OAuth**           | Per-user identity via your identity provider ([details](oauth-connectors.md)) |

With **OAuth**, each user signs in with their own Google or Azure AD identity
and queries run under that identity — row-level security and dataset-level
IAM apply per user. Tokens refresh silently in the background after the
initial consent.

See also:

* [Data connection](../../access.md#data-connection)
* [OAuth authentication for connectors](oauth-connectors.md)
* [BigQuery](https://cloud.google.com/bigquery/)
