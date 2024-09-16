# Caching function results

Caching is a process that stores multiple copies of data or files in a temporary storage location so that future requests for that data are served up faster than accessing the original source. 
Datagrok can cache the results of [functions](../../datagrok/concepts/functions/functions.md) on both client and server sides. This feature can drastically enhance the user experience and increase network efficiency.
Caching is especially useful for functions that are immutable and produce the same output for the same input.

## Client-side cache

Datagrok uses [IndexedDb](https://www.w3.org/TR/IndexedDB/) for client-side cache, and this type of cache has several limitations:

* Function output parameters should be scalar (integer, float or string) or of the following types - [dataframe](../../datagrok/concepts/table.md), [graphics or datetime](../../datagrok/concepts/functions/func-params-annotation.md).
* The function results have a cache limit of 100 MB. This means that the total size of the cached data for different input parameters of the same function cannot exceed 100 MB.
* Maximum record count in the **Client-side cache** is 100000.

To enable Client-side cache:

1. Go to **Settings** > **Cache**.
2. Check **Client-side cache** switch.
3. Click **Apply** button.

> Note: You can clear client-side cache by clicking **Clear client cache**.

## Server-side cache

Datagrok can cache the results of functions on the server. There are no restrictions on the size or type of parameters for caching. 
Both **Server-side cache** and **Client-side cache** can work together and enhance each other.

To enable Server-side cache:

1. Go to **Settings** > **Cache**.
2. Check **Server-side cache** switch.
3. Click **Apply** button.

## Using cache

You can cache results of almost all types of functions. But it really shines with [DataQueries](../../access/access.md#data-query) and [Scripts](../../compute/scripting/scripting.mdx).

To cache results of a particular function:

1. Add `meta.cache` annotation parameter in a function header. This parameter can take several values:
   * `client` - to use only **Client-side cache**
   * `server` - to use only **Server-side cache**
   * `all` - to use both types of cache.

2. Optionally add `meta.cache.invalidateOn` to specify when the cache is invalidated. This parameter accepts valid [cron expressions](https://www.adminschoice.com/crontab-quick-reference). If not specified, the cache will never be updated.

For example:

```
#name: Example
#language: python
#meta.cache: client
#meta.cache.invalidateOn: 0 * * * *
#input: string table [Data table]
#output: int result

...
```

To cache results of all queries under the specific connection:

1. Right-click the connection, select **Edit...** and check **Cache Results**. This activates client and server-side caches.
2. Optionally fill in the **Invalidate On** field or check **Cache Schema**: to cache database schema.

>Note: Caching the results of the database query is only a good idea when data stored in the database remains static and never changes.

Also, you can configure the cache for connection programmatically in the connection JSON definition. For example:

```json
{
  "name": "Northwind",
  "parameters": {
    "server": "dev.datagrok.ai",
    "port": 23306,
    "db": "Northwind",
    "cacheSchema": false,
    "cacheResults": true,
    "ssl": false,
    "connString": ""
  }
}
```
