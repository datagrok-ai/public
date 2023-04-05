# Create custom connectors

Datagrok connectors are plugins that provide a way for Datagrok to interact with external data sources. Datagrok platform supports [30+ popular databases](supported-connectors.md) right out of the box, with the list of supported databases continuously growing. You can also create your own connectors and seamlessly integrate them into the platform.

## How it works

The Grok Connect REST endpoints registered with the platform enable a server to work with data. The endpoints have the following methods:

* `getConnectors`: returns all connectors supported by the endpoint (one for each database type)
* `getSchema(connection)`: if applicable, returns database schema for the given connection
* `testConnection(connection)`: tests the connection
* `execute(query)`: executes the specified query
* `queryTable(structuredQuery)`: executes a structured query.

At startup, the server requests the list of supported connectors from each registered endpoint and creates a global list of supported connectors. The client requests the available connectors from the server, and the UI populates accordingly. When a client requests to query a database, the request is accepted by the server and then routed to the corresponding database connector.

## Adding a new connector

To add a new connector, follow these steps:

1. Create a [DataProvider](https://github.com/datagrok-ai/public/blob/5c9a8df6b7f1494ae5f666bd2aaf5c6d55bc4dee/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java) subclass. If you are adding a JDBC provider, subclass it from the [JdbcDataProvider](https://github.com/datagrok-ai/public/blob/5c9a8df6b7f1494ae5f666bd2aaf5c6d55bc4dee/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java).
1. Override the required methods, depending on the specifics of the database/provider.
1. Add it to the `Providers` list (list of supported providers) in the `connectors_info/DataProvider.java`.

## Testing a new connector

You have two options for testing your database connector:

1. Test with `GrokConnectShell`:
  *. Open `connectors/examples/query.json` (or another file that you need to use in the below configuration) and enter the query information: set the query itself, datasource equal to `descriptor.type` in your provider, credentials to server and server address, port, and the database name.
  * Configure your next run: use `Application` and these program arguments:
    * To save your output in the `connectors/examples/output.csv`, use `--q examples/query.json --o examples/output.csv`
    * To print your output, use `--q examples/query.json`.
  * Run `GrokConnectShell.main()` and check the output.
1. Connect Datagrok Docker image to the grok_connect.