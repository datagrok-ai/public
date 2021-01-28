<!-- TITLE: Exercises -->
<!-- SUBTITLE: -->

# Exercises

This is a set of programming exercises designed to make developers proficient with the
Datagrok platform. The exercises are organized as progressive steps, with tasks
of increasing complexity built on top of the previously completed steps.

During this course, we will be building support for
handling DNA nucleotide sequences. Let that not scare you, think of them as regular
strings that can only contain characters G, A, C, and T (and now you know the origins
of the "Gattaca" movie name). We will start with writing stand-alone functions, then 
automatically recognizing nucleotide sequences in the imported data, and then 
going all the way to custom visualizations, querying relational databases, 
predictive models, integration with the external utilities, data augmentation, and 
custom applications. 

Table of contents
  * [Setting up the environment](#setting-up-the-environment)
  * [Semantic types](#semantic-types)
  * [Querying databases](#querying-databases)

## Setting up the environment

Prerequisites: basic JavaScript knowledge

1. Install the necessary tools (Node.js, npm, webpack, datagrok-tools) following [these instructions](develop.md#getting-started)
2. Get a dev key for https://dev.datagrok.ai (you will work with this server) and add it by running `grok config`
3. Create a default package called `<Name>Sequence` using datagrok-tools: `grok create <Name>Sequence`
4. Upload it to the server: `grok publish dev --rebuild` (see other options [here](develop.md#deployment-modes))
5. Launch the platform and run the package's `test` function using different methods: 
    * via the [Functions](https://dev.datagrok.ai/functions?q=test) view
    * via the [Packages](https://dev.datagrok.ai/packages?) menu (find your package and run `test` from the 'Content' pane in the property panel)
    * via the console: press `~` and type `<Name>Sequence:test()` (note that the function's namespace corresponds to the package name)

## Semantic types

Prerequisites: basic JavaScript knowledge

Details: [How to Create a Semantic Type Detector](how-to/semantic-type-detector.md), [How to Add an Info Panel](how-to/add-info-panel.md)

You will learn: how to write semantic type detectors, how to develop context-specific data augmentation.  

1. Create a `complement` function that takes a nucleotide string and returns its complement.
   Essentially, change each character to the complementary one: A<=>T, G<=>C. 
   Run it and check whether everything works fine. 
2. Now, let's specify that this function is meant to accept not any string, but nucleotides only,
   and returns a nucleotide string as well. In order to do that, let's annotate both input and output parameters with the `dna_nucleotide` semantic type (add `{semType: dna_nucleotide}` after the parameter name). At this point, `dna_nucleotide` string does not have any meaning, but we will connect the dots later.
3. Define a `detectNucleotides` semantic type detector function as part of the special
   `detectors.js` file. It should check whether a column is a string column, and whether
   each string represents a nucleotide (hint: for best performance, 
   don't iterate over all column values, instead iterate on `column.categories`). When everything is done
   correctly, the `detectors.js` file will get loaded by the platform automatically, and the
   `detectNucleotides` function will be executed against every column in a newly added table.
4. Test your implementation by opening the following CSV file (or go Data | Text, and paste it there):
   ```
   sequence, id
   GATTACA, 1997
   ATTCGGA, 1984
   TTTAGGC, 2021 
   ```
   Hover over the `sequence` column header after the data is imported — if everything is done correctly,
   you will see 'quality: dna_nucleotide' in the bottom of the tooltip. Alternatively, you can 
   find this information if you click on the column and expand the 'Details' pane in the property panel on the right.
5. Now transform the previously created `complement` function into an info panel: tag it with `panel` and `widgets` tags
   and change the output type to `widget` (see an example [here](how-to/add-info-panel.md#functions)). This will instruct the platform to use the `complement` function for providing additional information for string values
   of the `dna_nucleotide` semantic type. To test it, simply open our test file, click on any cell
   in the `sequence` column, and find the `complement` property in the panel on the right.

## Querying databases

Prerequisites: basic SQL knowledge

Details: [Connecting to Databases](https://www.youtube.com/watch?v=dKrCk38A1m8&t=1048s), [How to Access Data](how-to/access-data.md)

In this exercise, we will work with a `northwind` Postgres database (in case the name sounds 
familiar, this is a demo database that Microsoft often uses for showcasing its technology).
The database is already deployed and is accessible from our server.

1. Navigate to the `Data | Databases | Postgres | northwind | orders` table
2. Make this table current by clicking on it, and explore its property panels on the right. The 
   `Content` pane should be showing first 50 rows of that table.  
3. Right-click on the table, and choose `New SQL Query...`
4. Execute the query and make sure it returns results.   
5. Modify the query to accept a `country` parameter, and return the sum of freights for the specified 
   country, grouped by `customerid`. Here is one way to do it:
   ```sql
   --input: string country
   select customerid, sum(freight)
   from public.orders
   where shipcountry = @country
   group by customerid
   ```
6. Make sure the query executes fine, when prompted, enter one of the countries in the input box (such as "USA").
   Run it the second time, notice that previously entered parameters could be quickly reused by clicking
   on the watch icon in the left bottom corner.
7. Give this query `ordersByCountry` name, and save it.
8. Use different ways to execute it:
    * Right-click on `Data | Databases | Postgres | northwind | ordersByCountry`, select `Run` from the menu
    * Click on `Data | Databases | Postgres | northwind | ordersByCountry`, expand the `Run` pane on the right, enter the country name and run it
    * Open console, see the results of the previous invocations. Copy-paste the corresponding command and 
      run it from the console.
9. Now, let's add this query to our package. Create a connection by running `grok add connection <name>`, then, as instructed [here](how-to/access-data.md#creating-queries), create the '.sql' file under the `queries` folder, and paste our query there. Give it a name by adding the `--name: ordersByCountry` line on top of it.
10. Deploy the package, launch the platform, find the query in the package, and run it.
11. Create a JavaScript function that has no parameters and returns a dataframe with 
    the results of the `ordersByCountry('USA')` call:
    ```javascript
    //name: getOrders
    //output: dataframe df
    export async function getOrders() {
      return await grok.data.query(`${packageName}:${queryName}`, { country: 'USA'});
    }
    ```
    There is another way to pass a country name to the query: you can provide a default value for the input parameter (see examples in the article [Parameterized Queries](../access/parameterized-queries.md)).

   
