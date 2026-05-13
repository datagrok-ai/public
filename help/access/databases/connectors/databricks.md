---
title: "Databricks"
---

Enables direct access to [Databricks](https://www.databricks.com/)
lakehouse platform through JDBC connectivity. This connector supports querying data
warehouses, data lakes, and analytics workloads hosted on **Databricks**.

## Connection parameters

```json
{
  "server": "",
  "port": "",
  "db": "",
  "connString": ""
}
```

## Authentication

Databricks connections support three authentication methods:

| Method                       | Use case                                                                             |
|------------------------------|--------------------------------------------------------------------------------------|
| **Personal Access Token**    | Shared identity via a Databricks PAT                                                 |
| **OAuth Client Credentials** | Shared service-principal identity via the M2M client-credentials grant               |
| **OAuth**                    | Per-user identity via your identity provider ([details](oauth-connectors.md))        |

With **OAuth**, each user signs in with their own Azure AD or OIDC identity
and queries run under that identity. On Azure AD, the token is passed
directly to the Databricks JDBC driver. On other OIDC providers, Datagrok
exchanges the token at the workspace's `/oidc/v1/token` endpoint using a
configured federation policy.

See also:

* [Data connection](../../access.md#data-connection)
* [OAuth authentication for connectors](oauth-connectors.md)
* [Databricks](https://www.databricks.com/)
