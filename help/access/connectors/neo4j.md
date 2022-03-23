<!-- TITLE: Neo4j -->

# Neo4j

This is a [connector](../data-connection.md#connectors) that provides access to the [Neo4j](https://neo4j.com/) graph
database via JDBC driver. Allows to query Neo4j using [Cypher](https://neo4j.com/developer/cypher-query-language)
language, and use results in dashboards, data augmentation panels, or via the [JS API](../../develop/js-api.md).

## Creating a connection

To create a new Neo4j connection, open the "Databases" pane (Open | Databases), right-click on Neo4j, and choose "Add
connection...". Then, enter the connection details, click TEST to make sure you've entered everything correctly, and
click SAVE. Note that by default, the connection will only be accessible to you, so you will have to share it in order
for others to use it. Check out [creating a connection](../data-connection.md#creating-a-new-connection) for more
details.

Once a connection is created, it's time to create a query.

## Creating a query

Right-click on the connection, select "Add query", and enter the Cypher query. Press F5 to run query (this is useful for
debugging); once you are satisfied with the results, press SAVE to save it on the server.

Note that you can introduce [parameters](../parameterized-queries.md) to the query. This is an incredibly powerful
concept that allows Cypher queries to be used as
[functions](../../overview/functions/function.md), or become info panels.

See [query editor](../data-query-view.md) for details.

![Add query](../../uploads/gifs/query-add.gif "Add query")

## Using a query

A query is a [function](../../overview/functions/function.md), and therefore could be used for multiple purposes in many
different contexts. Here are some of them:

* Executed manually
    * From the [query editor](../data-query-view.md)
    * From the "Functions" pane
    * As part of expression for [calculated columns](../../transform/add-new-column.md)
* As an [info panel](../../discover/info-panels.md)
* From [JS API](../../develop/js-api.md)

See also:

* [Data connection](../data-connection.md)
* [Neo4j](https://neo4j.com/)
