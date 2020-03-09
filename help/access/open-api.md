<!-- TITLE: OpenAPI -->
<!-- SUBTITLE: -->

# OpenAPI

[OpenAPI](https://swagger.io/docs/specification/about/), also known as swagger,
is a popular format that describes the structure of the server APIs so that machines 
can read the document and use the service.

Grok platform integrates with OpenAPI really well. Once a swagger file is imported
(you can simply drag-and-drop yaml file into the app), its content gets translated
to Grok connections, queries, and [functions](../overview/functions/function.md). All of them can be combined and used
in data jobs, calculations, [info panels](../discover/info-panels.md), executed from [console](../features/console.md), etc.   

Let's take a look at the example of AirNow, a site that tells you how clean or polluted
your outdoor air is. Once its swagger file is imported (which takes a second), 
the following connection and queries appear in Grok:

![AirNow connection](../uploads/features/open-api-airnow-connection.png "AirNow")

Click on the query to bring up its details in the property panel on the right:

![AirNow query](../uploads/features/open-api-airnow-query.png "AirNow")

Note that it became a regular [Grok query](data-query.md), which means all the features of the platform
are now applicable to it - in fact, end users won't even have to know where the data is coming
from, or what technology is used there. All data governance capabilities, such as data lineage,
history, and security, can be used. Query results can be automatically transformed, and there
is even a chat associated with it. But of course, the most important thing is that
it can be executed right from there. Let's fill in the required zip code,
hit RUN, and voila - the query is executed and results are shown in the table.

![AirNow results](../uploads/features/open-api-airnow-results.png "AirNow")

Is it important to understand that by using the OpenAPI Grok connector, it is possible to not
only retrieve data from external services, but do pretty much everything the service lets you
to do. For instance, you can automatically provision virtual machines via the Amazon AWS swagger,
create Jira tickets, or purchase items on eBay. 

See also:
* [Public datasets](public-datasets.md)
* [OpenAPI collection](https://apis.guru/browse-apis/)
* [Data query](data-query.md)
