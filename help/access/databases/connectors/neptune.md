---
title: "Neptune"
description: Connect Datagrok to the Amazon Neptune graph database service using SQL over JDBC.
keywords:
  - amazon neptune
  - aws neptune database
  - neptune jdbc driver
  - connect to neptune
  - graph database service
  - service region
---

Provides access to the [Neptune](https://aws.amazon.com/neptune/) database service using SQL queries via a JDBC driver.
This connector is served by the optional **grok_connect_extended** service. If it does not appear
under **Data > Databases**, ask your administrator to enable it (Helm: `grokConnectExtended.enabled: true`;
Docker Compose: start the `grok_connect_extended` profile).


## Connection parameters

```json
{
  "server": "",
  "port": "",
  "serviceRegion": "",
  "connString": ""
}
```

See also:

* [Data connection](../../access.md#data-connection)
* [Neptune](https://aws.amazon.com/neptune/)
