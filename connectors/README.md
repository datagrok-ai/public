# Datagrok Database Connectors

Out of the box, the platform comes with the data connectors for [30+ popular databases](../help/access/data-connection.md#connectors), and the list is constantly growing.
In addition to that, it is possible to develop your own data connectors, and seamlessly integrate them
into the platform.

## How it works

A server works with the "Grok Connect" REST endpoints registered with the platform. The endpoint's methods are:
* `getConnectors` - returns all connectors that the endpoint supports (one for database type)
* `getSchema(connection)` - if applicable, returns database schema for the given connection
* `testConnection(connection)` - tests the connection
* `execute(query)` - executes the specified query
* `queryTable(structuredQuery)` - executes a structured query

At startup, the server asks each registered endpoint for the list of supported connectors, and 
creates a global list of supported connectors. The client asks the server for the available connectors,
and populates the UI accordingly. Later on, when client makes a request to query a database, this request
gets accepted by a server, and then routed to the corresponding database connector.

## How to add a new database connector 

Create a [DataProvider](https://github.com/datagrok-ai/public/blob/5c9a8df6b7f1494ae5f666bd2aaf5c6d55bc4dee/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java) 
subclass. If you are adding a JDBC provider, subclass it from the 
[JdbcDataProvider](https://github.com/datagrok-ai/public/blob/5c9a8df6b7f1494ae5f666bd2aaf5c6d55bc4dee/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java)
Then, override the required methods depending on the specifics of the database/provider, and add it to 
the list of supported providers (see `DataProvider.Providers`).

For testing and debug purposes, a command-line 
["GrokConnectTest"](https://github.com/datagrok-ai/public/tree/master/connectors/grok_connect/src/test/java/grok_connect) 
application could be quite useful. 