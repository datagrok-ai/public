# Access databases

Databases can be accessed from a section `Databases` under the `Folder` 📁 icon on the left.

## Database connections

Datagrok provides a way to connect to pretty much any database out of the box. We support more than 30 of them, and its
simple to add a new connector with some custom development.

Adding a new connection is as simple as adding a [file share](). In a similar way, hit a
`New Connection` command on the left, choose the type of the provider, and then depending on the provider different
login credentials are required. Login credentials are stored in a separate secured vault in an encrypted way which can
only be decrypted with your login credentials — we thought about it a lot.

Also, all the connections to the databases as well as connections to file systems are also subject of privileges and
security checks. Each connection and sharable object has a panel called "Sharing". Expand it to see who has privileges.
It is easy to share it with a particular group or an individual and provide a level of access needed.

## Exploring a database

Once a connection to the database is created, there are multiple ways to connect, explore and work with that database.
Let's look at what we could do on the exploration side.

### Preview the table data

Under each connection we see a Node.js called `Tables`. It all depends on the provider, where some providers also
support `Schemas`. The simplest way is to simply click on the table and expand contents of it. By default, it shows the
first 50 rows, so it isn't going to take a lot of time to load the data preview. We find it a very useful tool for just
checking what is in the database. The content is in our dataframes and data spreadsheet, so that you can already do some
basic data profiling already at that data preview. If the data contains some molecules or other user-defined types, we
would see these data types rendered in the preview.

### Explore columns

In the below section properties for each column are displayed. In addition to some general information, such as names,
types and additional metadata, there is also a way to quickly inspect it with the basic descriptive statistics.

## Query the database visually

### Querying by aggregation

For example, let us do a couple of visual queries against the `Orders` table of the celebrated
`Northwind` database. Let us right-click on the `Orders` table to which we navigated and select a `Visual Query`
command. What this particular tool does is it lets you create an aggregation query against a particular table visually.

While the functionality for joining tables is not yet available and comes later, it is often possible to add a view to
the actual database and then interrogate this view from the `Visual Query`
tool.

Building an aggregate query is easy. Start with the `Measures` section. For example, we select an average of `freight`,
and the result appears instantaneously. It is a nice way to explore datasets which don't fit in the browser's memory, or
something where you know in advance what you are looking for.

It is possible to change a measure of the aggregation by right-clicking on it and selecting the measure in interest,
such as choosing a `sum` instead of an `avg`. Same applies to the column by which the aggregation is being computed.

To group by different columns, use the `Rows` section of the dialog. For example, let us group by `Country` and then
by `City`.

### Data pivot

On of the popular features when people query data is the ability to pivot it, essentially put values in columns. We
support this feature with the `Columns` field. As an example, let us select the ...

### Querying by joining

Another popular way of querying the database is the one when you start with particular table and then join particular
attributes using left join which are related to that table.

For example, ...

Right-click on the table and choose `Build Query`. Our platform figures out the schema of the database, and starting
from that table it adds all the tables that could be reached by following the foreign keys.