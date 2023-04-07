---
title: "OpenAPI"
---

[OpenAPI](https://swagger.io/docs/specification/about/), also known as Swagger, is a popular format that describes the structure of the server APIs allowing machines to read the document and use the service. Datagrok seamlessly integrates with OpenAPI, making it easy to connect to webservices and execute queries using the platform's features.

## Connecting to a webservice

To upload an OpenAPI document, drag and drop a YAML or JSON file into Datagrok. Datagrok detects different attributes of a Swagger file and automatically creates a data connection and associated queries:

| Swagger file    | Datagrok                                             |
|-----------------|------------------------------------------------------|
| title           | Data connection name                                 |
| description     | Data connection description                          |
| paths           | Data query is created for each path                  |
| summary         | Data query name                                      |
| parameters      | Data query parameters                                |

When the file import is complete, the data connection, along with all associated queries, appears in the **Webservices Manager** (**Data** > **Webservices**).

:::tip

If you can't find a Swagger JSON/YAML file link, inspect the API service's [Swagger-UI browser](https://swagger.io/tools/swagger-ui/) browser page. Use Chrome Developer Tools' "Network" tab. Usually, the file link contains `/v1` or `/v2`. Alternatively, add `/api/v2/api-docs` to the original service `<URI>`. This is the default location for Swagger JSON data in Swagger-UI browsers.

:::

### Parameters

Swagger format allows defining parameters both within individual paths and at the file level, with an option to use them in other paths. Datagrok effectively works with both approaches (or their combination) for both the standard Swagger attributes and any custom parameters that you add.

For example, consider the following AirNow YAML file, where the `date` parameter is defined within the path, and the remaining parameters are placed to `parameters:` section at the file level:

```yaml
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
    description: |
      If no reporting area is associated with the specified Zip Code, 
      return a forecast from a nearby reporting area within this distance (in miles).
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

To add a custom `grok-datetime-format parameter`, you can either apply it to the entire Swagger file or define it within specific paths:

Custom parameter applied to the entire Swagger file:

```yaml
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

Custom parameter applied to an individual path:

```yaml
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

Datagrok can work with both parameterized queries and queries without parameters. Furthermore, there are no restrictions in Datagrok on using parameters in different locations, including:

* path parameters, such as `/users/{id}`
* query parameters, such as `/users?role=admin`
* header parameters, such as `X-MyHeader: Value`.

### Editing Swagger files before import

To edit the original Swagger file provided by the service or enhance the file with simpler queries not present in the original Swagger file, we recommend using [Postman](https://www.postman.com/). You can import a Swagger JSON/YAML file into Postman for introspection, manipulation, and pruning using the "Import" button. If you need to remove some Swagger items, do it directly in Datagrok after uploading or importing it.

Usually, a Swagger file from the API service's Swagger UI works well with both Datagrok and Postman. Import issues, if any, usually arise due to Swagger version differences or parser inconsistencies. First, ensure the Swagger name is included under the `"info"` > `"title"` section. If import errors persist, consider these steps:

<details>
<summary> Issue: Datagrok loads the Swagger file successfully but `basePath` or `host` are missing along with the Swagger icon </summary>

Add this section to the file:

```
 "schemes": [
   "https",
   "http"
 ]
```

</details>

<details>
<summary> Issue: Postman or Datagrok can't open the file </summary>

Solution 1. Change `"swagger": "2.0"` to `"openapi": "2.0"`

Solution 2. If `"version"` isn't present in the original file, add a `"version"` section with an arbitrary version to the `"info"` section:

```
 "version": "1.0.0"
```

</details>

### Credentials

The Swagger file doesn't define access parameters for an OpenAPI service directly. Instead, it describes the types and access parameters in the `securityDefinitions:` section. In the AirNow example, the `api_key` type is used with the `"API_KEY"` parameter, following the provider's naming requirement:

```yaml
securityDefinitions:
  api_key:
    type: apiKey
    name: API_KEY
    in: query
```

:::note

In OpenAPI 3.x, `securityDefinitions` were renamed to `securitySchemes` and moved inside `components`. Datagrok supports both major versions (OpenAPI 2.x and OpenAPI 3.x).

:::

To specify access credentials for a newly created connection, you need to edit its setting within Datagrok. To do this, in the **Webservices Manager**, right-click a connection and select **Edit...** from its context menu. In the **Edit Connection** dialog, select the appropriate **Security** type from the dropdown and enter the **ApiKey**.

![AirNow connection](../uploads/features/swagger-security-definitions.png "AirNow")

:::note

Datagrok supports all types of secret access for Swagger, such as:

* basic authentication
* API key (as a header or a query string parameter), and
* OAuth 2 common flows (authorization code, implicit, resource owner password credentials, client credentials).

:::

When webservices don't require secret access, omit the `securityDefinitions:` block inside the Swagger file. See this [sample Swagger file without the security definitions](https://github.com/datagrok-ai/public/blob/master/packages/Swaggers/swaggers/countries.yaml)
.

### Packages

In addition to uploading Swagger files via drag-and-drop, Datagrok supports import of these files from [external packages](../develop/develop.md). Here's one such [package](https://github.com/datagrok-ai/public/tree/master/packages/Swaggers), which contains our Swagger demo files. The package includes many Swagger examples for different services in various formats (YAML, JSON).

When Swagger files are stored in this manner, they are imported to Datagrok (and new data connections appear) simultaneously as the corresponding package is published.

## Webservices Manager

The **Webservices Manager** provides a convenient interface for hierarchical browsing and managing connections and associated queries. To access the context actions for a connection or a query, right-click it. If you don't see a certain action, it may be due to insufficient permissions. In such cases, contact your Datagrok administrator for assistance.

Here's an example.

When you import the Swagger file for AirNow, a site that provides information about outdoor air quality, the following connection and queries appear in Datagrok:

![AirNow connection](../uploads/features/open-api-airnow-connection.png "AirNow")

The queries are now standard [Datagrok queries](access.md/#data-query), with all platform features available. For example, you can run queries using a convenient [UI for entering parameters](databases.md/#running-queries) or edit queries with [built-in editors](databases.md/#querying-data). Data governance and automation capabilities also apply.

Whenever you click a connection or query in the **Webservices Manager**, the [**Context Panel**](../datagrok/navigation.md#context-panel) to the right displays the object's properties and context actions. For example, when you click a query, the **Context Panel** lets you view the query's details, run, edit, or share it, and access other relevant information and options:

![AirNow query](../uploads/features/open-api-airnow-query.png "AirNow")

Running a query opens it in Datagrok as an [interactive grid](visualize/viewers/grid.md):

![AirNow results](../uploads/features/open-api-airnow-results.png "AirNow")

You can now explore the query results further, create [query views and dynamic dashboards](databases.md/#creating-views-for-query-results) and [share them with others](databases.md/#sharing-query-results).

:::tip

When using the OpenAPI connector in Datagrok, you're not limited to only retrieving data from external services; you can also perform any action that the service allows. For instance, you can automatically provision virtual machines using the Amazon AWS Swagger, create Jira tickets, or even make purchases on eBay.

:::

## Resources

[![Open Api](../uploads/youtube/data_access.png "Open on Youtube")](https://www.youtube.com/watch?v=dKrCk38A1m8\&t=3121s)
