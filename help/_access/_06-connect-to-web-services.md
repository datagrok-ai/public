# Connect to web services

In addition to all the data sources like files or databases, it's possible to use web services. We can ingest them in
multiple different ways. Of course it is always possible to write a custom JavaScript / TypeScript code
using [Datagrok API](../develop/js-api) to do exactly what is needed. However, for web services that provide an OpenAPI decription we
have a solution that works out of the box.

OpenAPI, sometimes referred to as Swagger, is an open specification that defines methods, or functions, that a
particular web server exposes, along with the parameters and other information such as how do you authorize against it,
and so on. The specification is pretty complicated and is evolving. We've completely supported OpenAPI 3.0. Not all
cases are done, but we do most of the usually used cases, and there's a way to intercept whatever you need and customize
where necessary.

Here is a collection of different OpenAPI services that we currently have in our repository. None of them were created
manually. They were introduced by importing OpenAPI / Swagger files.

This is a YAML-file which describes what the web service is doing. Let's simply drag-and-drop it into the platform. The
platform is smart enough to understand what the file has. After some processing we see that the connector has been
imported and appeared in the property panel. Alternatively we can go to `Web Services` and find the newly created
connection there.

Once expanded, we can see all the methods supported by the service. Each is a function just like a [query](../access/data-query) or
a [script](../compute/scripting). We can run it and the UI is formed automatically. Hint: it's possible to see the previous runs of the
function using a `Watch` icon.

They are subject to acess control and priveleges as well.

The result of the query to the web service comes as a dataframe table. For example, the actual output of the service may
be JSON or XML, but the output is a nice table. We do a lot atomatically underneath, such as unnesting of JSON structure
in order to represent results in a tabular manner. It doesn't work for all types of responses. There are a lot of
heuristics to determine the proper level of unnesting.

In general, a transformation to the query output may be applied. There's an ability to intercept the output of the
Swagger call.
