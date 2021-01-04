<!-- TITLE: OpenAPI -->
<!-- SUBTITLE: -->

# OpenAPI

[OpenAPI](https://swagger.io/docs/specification/about/), also known as swagger,
is a popular format that describes the structure of the server APIs so that machines 
can read the document and use the service.

Simple example of contents of a Swagger file in yaml format (Datagrok also supports json format for Swagger files) with a minimum set of attributes might look like this:

```
swagger: '2.0'
info:
  description: 'AirNow'
  title: AirNow
host: airnowapi.org
basePath: /aq
schemes:
  - http
paths:
  /observation/latLong/historical/:
    get:
      summary: Historical Observation By Latitude and Longitude
      operationId: historicalObservationByLatitudeAndLongitude
      produces:
        - text/csv
        - application/xml
        - application/json
      parameters:
        - name: date
          in: query
          required: false
          description: Date of forecast. If date is omitted, the current forecast is returned.
          type: string
          format: date-time
          grok-datetime-format: yyyy-MM-ddT00-0000
        - $ref: '#/parameters/latitude'
        - $ref: '#/parameters/longitude'
        - $ref: '#/parameters/distance'
      responses:
        '200':
          description: successful operation
          schema:
            type: array
            items:
              $ref: '#/definitions/Observation'
        '400':
          description: Invalid status value
grok-datetime-format: yyyy-MM-dd
parameters:
  distance:
     name: distance
     in: query
     required: false
     description: If no reporting area is associated with the specified Zip Code, return a forecast from a nearby reporting area within this distance (in miles).
     type: integer
     format: int32
  latitude:
     name: latitude
     in: query
     description: Latitude in decimal degrees.
     required: true
     type: number
     format: float
  longitude:
     name: longitude
     in: query
     required: true
     description: Longitude in decimal degrees.
     type: number
     format: float
securityDefinitions:
  api_key:
    type: apiKey
    name: API_KEY
    in: query

```

Grok platform integrates with OpenAPI really well. Once a swagger file is imported
(you can simply drag-and-drop yaml or json file into the app), its content gets translated
to Grok connections, queries, and [functions](../overview/functions/function.md). All of them can be combined and used
in data jobs, calculations, [info panels](../discover/info-panels.md), executed from [console](../overview/navigation.md#console), etc.   

After  swagger file is uploaded to the platform, based on what is described in it, new [data connection](data-connection.md) will be created in Datagrok.

You can find this connection in Connections Tree under the source "Web". Also in the user interface of Datagrok there is a special view "Web Services", which displays only connections for OpenAPI. (Located in the "Data" section on the left sidebar next to "Databases")

Also, you can find all  created connections to OpenAPI in the Data Connections browser. ("Manage" section on side bar and then "Connections" view)

Let's consider how Datagrok interprets individual attributes of the Swagger file:

| In Swagger File | In Datagrok                                          |
|-----------------|------------------------------------------------------|
| title           | [Data connection](data-connection.md) name           |
| description     | [Data connection](data-connection.md) description    |
| paths           | [Data query](data-query.md) is created for each path |
| summary         | [Data query](data-query.md) name                     |


In some cases, the standard attributes of Swagger format are not enough to import file to Datagrok.

For example, for the correct interpretation of datetime format in the corresponding parameters, an additional ```grok-datetime-format``` field inside the Swagger file is used, which is not part of the standard Swagger format, but is used only for more correct import into Datagrok in rare cases.
As we can see from the Swagger file example above, this field can be used at the level of entire Swagger file and be redefined at the level of individual parameters within individual paths.


Here it is defined at the level of the entire Swagger file:
```
swagger: '2.0'
info:
  description: 'AirNow'
  title: AirNow
host: airnowapi.org
...

grok-datetime-format: yyyy-MM-dd
parameters:
...
```


And here it is overridden at the level of an individual path (for "data" parameter):
```
swagger: '2.0'
info:
  description: 'AirNow'
  title: AirNow
host: airnowapi.org
...
paths:
  /observation/latLong/historical/:
    get:
      summary: Historical Observation By Latitude and Longitude
      operationId: historicalObservationByLatitudeAndLongitude
      produces:
        - text/csv
        - application/xml
        - application/json
      parameters:
        - name: date
          in: query
          required: false
          description: Date of forecast. If date is omitted, the current forecast is returned.
          type: string
          format: date-time
          grok-datetime-format: yyyy-MM-ddT00-0000
 ...
```

## Parameters

Swagger format supports the direct definition of parameters within each path and definition of parameters at the entire file level and their further use in any path.

Datagrok perfectly understands both options or their combination. 
In our AirNow example, the ```date``` parameter is defined inside the path and the rest of parameters are placed to ```parameters:``` section at the level of entire file.

It is important to mention that Datagrok can work both with parameterized queries with and with queries without parameters.

Also in Datagrok there are no restrictions on the use of parameters in different places, such as:

- path parameters, such as /users/{id}

- query parameters, such as /users?role=admin

- header parameters, such as X-MyHeader: Value

## Credentials 

Your access parameters to one or another OpenAPI service are not defined inside the Swagger file. 

Swagger file only describes the type and access parameters. 
The ```"securityDefinitions:``` section is used to determine the type and access parameters. 

In the AirNow example, the api_key type is used with a parameter named "API_KEY" (external service itself provides the parameter naming requirement):

```
securityDefinitions:
  api_key:
    type: apiKey
    name: API_KEY
    in: query

```

Datagrok works great with all types of secret access that the Swagger format supports:

- Basic authentication

- API key (as a header or a query string parameter)

- OAuth 2 common flows (authorization code, implicit, resource owner password credentials, client credentials)

After you have imported the Swagger file into Datagrok, you need to specify your access in the settings of corresponding [Data connection](data-connection.md).
To do this, find your connection in the connections tree (or in connection browser, special "Web Services" view, etc.) and open settings dialog for it:

![AirNow connection](../uploads/features/swagger-security-definitions.png "AirNow")

There are times when services do not require secret access. In situations like this, just don't specify the appropriate block ```securityDefinitions:``` inside the Swagger file and Datagrok will work fine with this.
An example of Swagger file without security definitions can be seen in our public repository, just click [here](https://github.com/datagrok-ai/public/blob/master/packages/Swaggers/swaggers/countries.yaml). 

## Example

Let's see what happens after we drag-and-drop the Swagger for AirNow, a site that tells you how clean or polluted your outdoor air is.
Once its Swagger file is imported (which takes a second), the following connection and queries appear in Grok:


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

## Packages

Swagger files can be imported into Datagrok not only from local storage by drag-and-drop.
These files can be also stored in external [packages](../develop/develop.md). 


You can check out one of these packages [here](https://github.com/datagrok-ai/public/tree/master/packages/Swaggers). 
This package contains our Swagger demo files. It contains enough Swagger examples for different services in different formats (yaml, json). 

In this variant  of storing Swagger files, they will be imported to the platform (new [data connections](data-connection.md) will appear) at the same time as the corresponding package is published.


## Videos

<iframe width="560" height="315" src="https://www.youtube.com/embed/dKrCk38A1m8?start=3121" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

See also:

  * [Public Datasets](public-datasets.md)
  * [OpenAPI Collection](https://apis.guru/browse-apis/)
  * [Data Query](data-query.md)
