<!-- TITLE: &#8204;Exercises -->
<!-- SUBTITLE: -->

# Exercises

These programming exercises are designed to help developers get proficient with the Datagrok platform. The exercises are
organized as progressive steps, with tasks of increasing complexity built on top of the previously completed steps.

During this course, we will be building support for handling DNA nucleotide sequences. Let that not scare you, think of
them as regular strings that can only contain characters `G`, `A`, `C`, and `T` (and now you know the origins of the
"Gattaca" movie name). We will start with writing standalone functions, then automatically recognizing nucleotide
sequences in the imported data, and then going all the way to custom visualizations, relational databases querying,
predictive models, integration with the external utilities, data augmentation, and custom applications.

## Table of contents

* [Setting up the environment](#setting-up-the-environment)
* [Semantic-types](#exercise-1-semantic-types)
* [Scripting and functions](#exercise-2-scripting-and-functions)
  * [Scripting with server functions](#scripting-with-server-functions)
  * [Modifying dataframes with scripts](#modifying-dataframes-with-scripts)
  * [Scripting with client functions](#scripting-with-client-functions)

<!---
* [Composing functions](#composing-functions)
* [Composing a JavaScript and a Python function](#composing-a-javascript-and-a-python-function)
* [Composing two JavaScript functions](#composing-two-javascript-functions)
* Composing functions in the package
--->

* [Querying databases](#exercise-3-querying-databases)
* [Creating a scripting viewer](#exercise-4-creating-a-scripting-viewer)
* [Transforming dataframes](#exercise-5-transforming-dataframes)
* [Custom cell renderers with 3-rd party JS libraries](#exercise-6-custom-cell-renderers-with-3-rd-party-js-libraries)
* [Accessing Web services with OpenAPI](#exercise-7-accessing-web-services-with-openapi)
* [Creating an info panel with a REST web service](#exercise-8-creating-an-info-panel-with-a-rest-web-service)
* [Enhancing Datagrok with dialog-based functions](#exercise-9-enhancing-datagrok-with-dialog-based-functions)

<!---
* Creating an application
* Accessing Web services in JavaScript with REST
* Creating a custom JavaScript viewer
* Optimizing a custom cell renderer with a cache
* Extending Datagrok with info panels
* Customize packages with properties
* Persisting user sessions and tables
* Using WebAssembly with Datagrok functions
* Webpack packages with WebAssembly
* Using Web Workers for background computations
--->

## Setting up the environment

*Prerequisites:* basic TypeScript or JavaScript knowledge.
*Useful links:*

* [Datagrok tools](https://www.npmjs.com/package/datagrok-tools)
* [Naming conventions](https://datagrok.ai/help/develop/develop#naming-conventions)

1. Install the necessary tools (Node.js, npm, webpack, datagrok-tools) following
   [these instructions](../set-up-environment.md)
2. Create a branch from master at [GitHub](https://github.com/datagrok-ai/public/branches) or using your IDE; Use your
   credentials as a name of the branch
3. Get a dev key for [Dev Server](https://dev.datagrok.ai) (you will work with this server) and add it by
   running `grok config`. Open [https://dev.datagrok.ai/u](https://dev.datagrok.ai/u), click on `Developer key`, copy
   the `grok` command and execute it to add the key to your config
4. Create a default package in your branch [called](https://datagrok.ai/help/develop/develop#naming-conventions)
   `<yourFirstName>-sequence` using datagrok-tools:
   `grok create <yourFirstName>-sequence` with specifying the `--ts` option to create a package with TypeScript
   configuration
   (if you are new to TypeScript, you can specify the `--js` option);
   Note that detectors.js file should be in JavaScript anyway.
   Also you can add `--eslint` option to add eslint checker feature to the package
5. Run `npm install` to link the dependencies mentioned in `package.json` file of your package
6. Upload it to the server: run `webpack` and `grok publish dev` (see other
   options [here](../develop.md#deployment-modes))
7. Launch the platform and run the package's `info()` function using different methods:

* via the [Functions](https://dev.datagrok.ai/functions?q=info) view
* via the [Packages](https://dev.datagrok.ai/packages?) menu (find your package, click on it and run `info()`
  from the `Functions` pane in the property panel on the left)
* via the [console](../../datagrok/navigation.md#console): press `~` key anywhere inside Datagrok, the Console will
  appear to the right; execute `<loginName>Sequence:info()` there. The identifier used as package name (before ':') will
  be obtained by transformation kebab style of folder name to camel style, or can be specified directly with
  attribute `friendlyName` in `package.json` file.

As a result of the function execution you should see an info notification with url of package's webRoot.

## Exercise 1: Semantic types

*Prerequisites:* basic TypeScript or JavaScript knowledge.

Details: [How to define semantic type detectors](../how-to/define-semantic-type-detectors.md),
[How to add an info panel](../how-to/add-info-panel.md).

You will learn: how to write semantic type detectors, how to develop context-specific data augmentation.

1. Create a `complement` function in `src/package.ts` which takes a nucleotide string and returns its complement:

    ```javascript
    //name: complement
    //input: string nucleotides
    //output: string result
    export function complement(nucleotides): /*type*/ {
        // your code goes here
    }
    ```

   Note that comments on the top of the function declaration are crucial for running it on the platform. They determine
   the function name, the input and output types. Essentially, change each character to the complementary one: `A <=> T`
   , `G <=> C`. Run it and check whether everything works fine.

2. Now, let's specify that this function is meant to accept not any string, but nucleotides only, and to return a
   nucleotide string as well. In order to do that, let's annotate both input and output parameters with
   the `dna_nucleotide` semantic type:

    ```javascript
    //input: string nucleotides {semType: dna_nucleotide}
    //output: string result {semType: dna_nucleotide}
   ```

   At this point, `dna_nucleotide` string does not have any meaning, but we will connect the dots later.

3. <a name="detectors"></a> Define a `detectNucleotides` semantic type detector function as part of the special
   `detectors.js` file.

    ```javascript
   class <yourFirstName>SequencePackageDetectors extends DG.Package {

     //tags: semTypeDetector
     //input: column col
     //output: string semType
     detectNucleotides(col) {
         // your code goes here
     }
   }
   ```

   It should check whether a column is a string column, and whether each string represents a nucleotide. If condition is
   met, it should return `"dna_nucleotide"` string.

   For best performance, don't iterate over all column values, instead
   iterate [on `column.categories`](https://datagrok.ai/help/develop/advanced/data-frame#work-with-categories)
   . Full Datagrok Column type API could be found [here](https://datagrok.ai/js-api/classes/dg.Column).

4. Upload your package to `dev.datagrok.ai` using `grok publish dev`
   command. When everything is done correctly, the `detectors.js` file will get loaded by the platform automatically,
   and the
   `detectNucleotides` function will be executed against every column in a newly added table.

5. Reload the `dev.datagrok.ai` page to use updated version of your package.

6. Test your implementation by opening the following CSV or TXT file (or go to `📁 (Data) | Text`
   and paste it there). Make sure you click on DONE (this will trigger semantic types detection):

   ```
   sequence, id
   GATTACA, 1997
   ATTCGGA, 1984
   TTTAGGC, 2021
   ```

   Hover over the `sequence` column header after the data is imported — if everything is done correctly, you will
   see `quality: dna_nucleotide` in the bottom of the tooltip:
   ![exercises-semantic-tooltip](exercises-semantic-tooltip.png)
   Alternatively, you can find this information if you click on the column and expand the 'Details' pane in the property
   panel on the right.
7. Now let’s put the result of the previously created `complement` function into an [info panel](../how-to/add-info-panel.md):
   Create function `complementWidget` and add special comments block to allow Datagrok system recognise it and upload
   properly (see an example [here][014]).

   ```javascript
    //name: complementWidget
    //tags: panel, widgets
    //input: string nucleotides {semType: dna_nucleotide}
    //output: widget result
    //condition: true
   ```

   The `panel` and `widgets` tags and output type `widget` allows Datagrok to determine how the result of
   `complementWidget` function will appear in the system. Listed above block of comments will
   instruct the platform to use the `complementWidget` function for providing additional information for string values
   of the `dna_nucleotide` semantic type. To test it, simply open our test file, click on any cell in the `sequence`
   column, and find the `complementWidget` property in the panel on the right as it is shown on screenshot:
   ![exercises-complement-data-panel](exercises-complement-data-panel.png)

## Exercise 2: Scripting and functions

### Scripting with server functions

*Prerequisites:* basic Python knowledge.

*Details:* [Scripting](../../compute/scripting.md), [Dev Meeting 1 | First-class functions][015]

*You will learn:* how to create and invoke Datagrok scripts in data science languages like R and Python.

In this exercise, we will count occurrences of a given subsequence in a nucleotide sequence, using Python.

1. Open Datagrok and navigate to `Functions | Scripts | Actions | New Python Script`.
2. Observe a default script created for you. All script attributes are specified in the beginning in comments. There we
   have the script name, language, one input value of type `dataframe`, and one output value of type `int`. The script
   simply computes number of cells in the dataframe.
   [Dataframe](../how-to/build-an-app.md) is a high-performance, easy to use tabular structure with strongly-typed
   columns of different types (supported types are: `string`, `bool`, `int`
   , `bigint`,
   `double`, `qnum` and `datetime`). In this exercise, we only see a dataframe as is in the default script; there is
   another exercise to learn manipulating dataframes in JavaScript.
3. Run the script to get a hint for creating an input file. An attribute `#sample: cars.csv`
   is responsible for it. To open a default input file `cars`, click the `Star` icon in the top menu.
4. Run the script again and proceed to the Datagrok's console. As in Quake, it's available by pressing a `~` button
   anywhere inside Datagrok. In the console, you would see the script execution result. Just one line above the result
   you could see the console's command to execute the script. Enter it again to the console to get the same result
   (but in console you should specify script with namespace prefix as `<yourLogin>:<script_name>`).
5. Let's modify the script to solve the task of counting sequence occurrences. Add a new preamble:
   (use any `#description` you like). Spaces are not allowed between '#' and attribute name:

    ```python
    # name: CountSubsequencePython
    # language: python
    # input: string sequence
    # input: string subsequence
    # output: int count
    ```

   In the body, implement a Python function counting all occurrences of a given `subsequence` in a `sequence`. Return
   a `count` the same way as in the default script from p. 2.

6. Run the script function, provide input values in the dialog and get to the console to see the result. Now run the
   script function again through the console completely, passing different arguments values:
   ```<yourLogin>:CountSubsequencePython('ATGATC', 'A')```. You can find your login inside the profile page between name
   and email (under avatar), or in the profile URL:
   `https://dev.datagrok.ai/u/<yourLogin>/summary`.
7. Let's apply `CountSubsequencePython` to the input dataframe using Datagrok UI. Open a table — say, let's go
   for `sars-cov-2.csv`. Navigate to `Data | Files` and open
   `Demo Files / bio / sars-cov-2.csv`. Navigate to a menu item `Edit | Add New Column...`
   and click it. Type in your expression using the function you've just previously created:
   ![exercises-add-new-column](exercises-add-new-column.png)
   Observe how the `Preview Result Columns` change while you are modifying the expression. There, notice a
   namespace `<yourLogin>` as part of a qualified function name `<yourLogin>:<functionName>`,
   `JDoe:CountSubsequencePython` in this case. Namespaces are used through Datagrok very commonly. In general, there shall
   be no case where you would call a function without specifying a namespace. Datagrok namespaces originate from the
   names of packages, projects, and users, and always qualify a resource name, be it a package, a function, a connection
   or a query. Now hit "Ok" and have the new column inserted to the dataframe.

<!-- TODO: Update the dialog -->

### Modifying dataframes with scripts

*Prerequisites:* basic Python knowledge.

*You will learn:* how to manipulate tables, which we usually call dataframes, using a server scripting language, expand
dataframes with newly computed values, and modify the dataframes.

In the previous exercise we learnt a fast method to apply a function to a table and produce a new column in it. Another
means to introduce new columns to the dataframes is to programmatically manipulate dataframes right in scripts. Let's
repeat what we've achieved in the last point of the previous exercise, now with more scripting.

1. Let's create a different kind of our `CountSubsequencePython` function, now called `CountSubsequencePythonDataframe`.
   While the original function could only operate on a single row, the new function shall operate on the entire
   dataframe. To start with, the function's Datagrok signature should look as follows:

    ```python
    # name: CountSubsequencePythonDataframe
    # language: python
    # input: dataframe sequences
    # input: column columnName
    # input: string subsequence = "acc"
    # output: dataframe result {action:join(sequences)}
    ```

   This function takes as an input a dataframe with a column containing nucleotide sequences, named as a value of
   `columnName`, a nucleotide subsequence `subsequence` being sought, and outputs an input dataframe with a new column
   *appended* to it, containing numbers of subsequence occurrences in the nucleotide sequences. Say, for a table on the
   left the following table on the right should be produced for a subsequence `acc` being sought:

    <table>
    <tr><td>

   | GenBank    | ID         |
   |------------|------------|
   | MT079845.1 | ctacaagaga |
   | MT079851.1 | attaaaggtt |
   | MT326187.1 | gttctctaaa |

    </td><td>

   | GenBank    | ID         | N(acc) |
   |------------|------------|--------|
   | MT079845.1 | ctaccagaga | 1      |
   | MT079851.1 | attaaaggtt | 0      |
   | MT326187.1 | gttctctacc | 1      |

    </td></tr>
    </table>

2. Implement a function `CountSubsequencePythonDataframe`. Assume the `result` is a Python dataframe with just this one
   column `columnName`. After the `result` column is computed and returned from the server to the client, based on
   the `join` instruction, `result` will be  *appended* to the existing input dataframe `sequences`. As this is
   performed purely on the client, we save the bandwidth without needing to return a copy of a dataframe which we
   already passed to the server.

    * Use Pandas dataframes as `pd` to access the input dataframe and create an output dataframe
    * You don't need to import `pandas`, Datagrok does this automatically: to each Python script it adds a preamble with
      most popular imports (`os`, `io`, `json`, `pandas as pd`, `requests`
      , `datetime`, `timedelta`)
    * Note that the column `columnName` is just a string with a column name passed to a script, not an actual column content

3. Run the function with a "Play" button on top of the function window. The dialog will prompt you to select a
   dataframe. Navigate to a "Data" view (first button on the left sidebar) and open a file with nucleotide sequences
   (say, `Demo Files / bio / sars-cov-2.csv` available at public.datagrok.ai). Go back to the `Run Function` dialog to
   select the opened dataframe.

4. Now choose a column with nucleotide sequences from the dropdown. Notice how the list of columns is automatically
   formed for the selected dataframe. Finally, run the function to get the resulting dataframe.

5. As for modifying the dataframes in general, just consider removing the `{action:join}` option and do whatever is
   needed to the output dataframe `result` before the end line of the script. This will return exactly the `result`
   dataframe after all modifications.

6. Consider that the function may have several outputs. In case you return two dataframes, both will appear in the
   Datagrok interface. There's also a special syntax in Datagrok JS API to call functions which return several
   parameters, we'll review this in one of the following exercises.

### Scripting with client functions

*Prerequisites:* basic JavaScript knowledge.

*You will learn:* how to create and invoke Datagrok JavaScript scripts.

1. Go to `Functions | Scripts` and hit `New JavaScript Script`.
2. Implement the function `CountSubsequenceJS` in JavaScript, which does the same as
   [`CountSubsequencePython`]. Follow the same conventions on the parameters in the comments block and returning a
   result via a variable.
3. Run `CountSubsequenceJS` using the `Play` button; using the console. From same console, run `CountSubsequencePython`
   yet again. You can notice that both Python and JS versions of our function, implemented as scripts, are homogeneous
   functions in Datagrok. It's also possible to call them in a uniform
   fashion [using our JavaScript API](../../compute/scripting.md#running-a-script).
4. Don't forget to save these two scripts. We would re-use parts of them in the following exercises.

The difference between the two scripts is that the first, `CountSubsequencePython`, runs on our server by
a [compute virtual machine](../admin/infrastructure.md#compute-components), whereas the second, `CountSubsequenceJS`,
runs directly in the browser. To run `CountSubsequencePython`, Datagrok passes the script arguments over the network and
fetches back the result to the browser.

<!-- Redo with a simpler stance

## Composing functions

### Composing a JavaScript and a Python function

*You will learn:* how to invoke arbitrary Datagrok functions in JavaScript and augment tables.

1. Open Datagrok and navigate to `Functions | Scripts | New Python Script`.

3. Let's create a wrapping function `CountSubsequenceTableAugment` in JavaScript which would do all the technical
   job. We would use it to give the new column with counts a proper name, and to *augment* the input dataframe
   with the newly computed column. Make it look like this:

    ```javascript
    //name: CountSubsequenceTablePythonAugment
    //language: javascript
    //input: dataframe df
    //input: column colName
    //input: string subseq = ATG

    grok.functions.call(
      "<yourLogin>:CountSubsequenceTablePython", {
        'inputDf': df,
        'inputColName': colName,
        'outputColName': `N(${subseq})`,
        'subseq': subseq
      }).then((resultDf) => {
        df.columns.insert(resultDf.columns.byIndex(0));
      });
    ```

   Follow how we match input parameters to their values in `grok.functions.call` using a JSON object.
   In the `.then` [continuation]() we manipulate the original dataframe by inserting a new column which
   we get as a result of the script's execution.

4. Let' prepare a visual layout before running our script. Navigate to `Data | Files` and open
   `Demo Files / bio / sars-cov-2.csv`. Then click on `Windows` and activate `Menu`, this will let
   you move windows around. Stick the `Sequences` table to the side of the `CountSubsequenceTable`
   window and leave enough space in the table to see the new column coming when running the script.

5. Run the script and check the new column appears in the grid.

6. Add `CountSubsequenceTablePythonAugment` and `CountSubsequenceTablePython` as part of the package
   `<yourFirstName>-sequence` prepared in ["Semantic types"](#exercise-1-semantic-types) exercise. Deploy the package,
   reload Datagrok and run `CountSubsequenceTablePythonAugment` from the package.

7. Notice that we don't need the entire input dataframe to run the script, just one column. Optimize the 2 scripts
   of this exercise to only take this one column as an input, and thus optimizing the data network roundtrip.

8. In this example, we've manually and programmatically expanded the original dataframe with one new column.
   As this is a common action in many scenarios, the same result in Datagrok can be achieved with setting
   `{action:join(inputDf)}` on the output dataframe:
   ```
   #output: dataframe adfwithanewcolumn {action:join(inputdf)}
   ```
   This is useful not only for brevity, but also for saving bandwidth: the resulting dataframe will be joined with
   the origin dataframe on the client side, only the necessary new data shall be passed over to the client.
   For example, this is how to convert an existing column into a new one and add it to the original dataframe:

   ```
   #language: python
   #input: dataframe inputdf
   #output: dataframe outputdf {action:join(inputdf)}

   outputDf = inputDf[['someColumn']]
   outputDf['newColumn'] = 2 * outputDf['someColumn']
   ```
   Note that the `someColumn` column shall pop up as a new column in the resulting dataframe, as well,
   under a name `someColumn (2)`. If you don't want this, you should drop this column in the
   `outputDf`.

### Composing two JavaScript functions

8. Let's repeat this augmentation for a JavaScript function. First, delete the newly created `N(ATG)` column
   by clicking with a right mouse button on it and selecting `Remove`.

9. Create a function `CountSubsequenceTableJS`, which has the following structure:
    ```javascript
    //name: CountSubsequenceTableJS
    //language: javascript
    //input: dataframe inputDf
    //input: string inputColName
    //input: string outputColName
    //input: string subseq
    //output: dataframe outputDf

    let newCol = DG.Column.fromType(
      DG.COLUMN_TYPE.INT, outputColName, inputDf.rowCount);
    for (let i = 0; i < newCol.length; ++i) {
      const seq = inputDf.get(inputColName, i);
      const count = ...; // your subsequence counting here
      newCol.set(i, count);
    }
    outputDf = DG.DataFrame.fromColumns([newCol]);
    ```
   Grasp a number of techniques with Datagrok JS API:
     * creating a new column of a given type
     * iterating through a dataframe by a row index
     * creating a dataframe from an array of columns

10. In `CountSubsequenceTablePythonAugment`, replace `CountSubsequenceTablePython` to `CountSubsequenceTableJS`,
    rename itself to `CountSubsequenceTableJSAugment` and run it. Check that the exact same new column is produced
    as it was for a Python version.

11. Add `CountSubsequenceTableJS` and `CountSubsequenceTableJSAugment` as part of the
    package `<yourFirstName>-sequence` prepared in ["Semantic types"](#exercise-1-semantic-types) exercise.
    Deploy the package, reload Datagrok and run `CountSubsequenceTableJSAugment` from the package.

12. In contrast to `CountSubsequenceTablePythonAugment` running in the browser and `CountSubsequenceTablePython`
    running on the server, now both `CountSubsequenceTableJSAugment` and `CountSubsequenceTableJS` run
    in the same browser tab. Also, the dataframe-typed arguments passed from one JS function to another
    are just references to one JS object. Would it make sense to optimize the 2 functions in the same fashion
    as in p.6 of the previous exercise? If not, show a suitable optimization.

In these two exercises, we could try another approach. Instead of passing the column to
and forming the new column in the function being called, we could just populate the new column
in a loop inside the `CountSubsequenceTable` by calls to a function operating on row data, such as
the one created in ["Scripting and functions"](#scripting-and-functions):
`CountSubsequencePython(seq, subseq)`.

However, this incurs a substantial overhead in a Python version.
For a table of 10 rows we'd have to call the server scripting 10 times with a network roundtrip
to deliver parameters to [compute virtual machine](). This will be less overhead for a JavaScript version,
as all will happen in the browser, but still not as optimal as the case where we do just one call to a
nested script.

-->

## Exercise 3: Querying databases

*Prerequisites:* basic SQL knowledge

*Details:* [Connecting to Databases](https://www.youtube.com/watch?v=dKrCk38A1m8&t=1048s),
[How to Access Data](../how-to/access-data.md)

*Note:* Editing an existing data query requires the respective access permission. You might need to request one.

In this exercise, we will work with a `northwind` PostgreSQL database (in case the name sounds familiar, this is a demo
database that Microsoft often uses for showcasing its technology). The database is already deployed and is accessible
from our server.

1. Navigate to the `Data | Databases | PostgreSQL | northwind | Schemas | public | orders` table
2. Make this table current by left-clicking on it, and explore its property panels on the right. The
   `Content` pane should be showing first 50 rows of that table.
3. Right-click on the table, and choose `New SQL Query...`
4. Execute the query and make sure it returns results.
5. Modify the query to accept a `country` parameter, and return the sum of freights for the specified country, grouped
   by `customerid`:

   ```sql
   --input: string country
   select customerid, sum(freight)
   from public.orders
   where shipcountry = @country
   group by customerid
   ```

6. Run the query, enter one of the countries in the input box (such as `USA`, without quotation marks or apostrophes).
   Run it the second time, notice that previously entered parameters could be quickly reused by clicking on the watch
   icon in the left bottom corner of the dialog window.
7. Rename this query from your name to `ordersByCountry`, and save it.
8. Try different ways to execute it:

    * Right-click on `Data | Databases | PostgreSQL | northwind | ordersByCountry`, select `Run` from the context menu,
      enter the country name, and run it
    * Click on `Data | Databases | PostgreSQL | northwind | ordersByCountry`, expand the `Run` pane on the right, enter the
      country name and run it
    * Open console by pressing `~` key, see the results of the previous invocations. Copy-paste the corresponding command
      and run it from the console.

9. Now, let's add this query to our package. Create a connection by running `grok add connection <yourFirstName>`, then,
   as instructed [here](../how-to/access-data.md#creating-queries), create the '.sql' file under the `queries`
   folder, and paste our query there. Give it a name by adding the `--name: ordersByCountry` line on top of it.
10. Deploy the package, launch the platform, find the query in the package, and run it.
11. Create a JavaScript function (in `src/package.js`) that has no parameters and returns a dataframe with the results
    of the `ordersByCountry('USA')` call:

    ```javascript
    //name: getOrders
    //output: dataframe df
    export async function getOrders() {
      return await grok.data.query(`${packageName}:${queryName}`, { country: 'USA'});
    }
    ```

    There is another way to pass a country name to the query: you can provide a default value for the input parameter
    (see examples in the article [Parameterized Queries](../../access/parameterized-queries.md)).

## Exercise 4: Creating a scripting viewer

*Prerequisites:* basic Python knowledge, [matplotlib](https://matplotlib.org/) or a similar library

*Details:* [Scripting](../../compute/scripting.md)
, [Scripting Viewer](../../visualize/viewers/scripting-viewer.md),
[Creating a scripting viewer (video)](https://www.youtube.com/embed/jHRpOnhBAz4).

*Amino acids counting task.* In this exercise, we'd use a Python script to generate a histogram
(a distribution plot) for amino acids occurring in a column of nucleotide sequences. Amino acids are simply triples of
nucleotides from which the nucleotide DNA sequence is made of. These are also called triplets, or codon-codes. As there
are 4 letters `G`, `A`, `T`, `C`, there are 4 to power 3 protein amino acids: `GTA`, `AGC`, `TTT`, and so forth.

We don't know at which starting point each nucleotide sequence was cut: it could either be a junction of two triplets,
or one-third or two-third of a triplet. Therefore, we'd count in our statistics for all three possible cuts, starting
the reading frame off at offsets 0, 1, and 2 from the beginning of the nucleotide sequence.

Say, we are given a sequence `TTTAATTACAGACCTGAA`. We start to count triplets *without overlap* from an offset 0 first,
getting: `TTT`, `AAT`, `TAC`, `...`, `GAA`. Then we move to the offset 1, getting: `TTA`, `...`, `CTG`. Lastly, we move
to the offset 2 and get: `TAA`, `...`, `TGA`. In the histogram we'd count for all these triplets.

First, let's explore how scripting viewer works.

1. Open a `demog` demo file with demographic data. It is located at `Data | Files | Demo Files | demog.csv`.
   `Data` corresponds to the first button from the top of the Datagrok sidebar. Make sure the table view with the data
   appears.
2. Activate the top menu from the sidebar, using a `Windows | Menu` switch.
3. In this menu, hit `Add | Scripting Viewers | Python | Scatter Plot`.
4. See that the viewer appeared on the right, telling though it is "Unable to plot with current settings".
5. Proceed to the viewer properties by hitting on the gear icon in the viewer's title.
6. Make sure the chosen values for "Data" are `HEIGHT` for `X`, `WEIGHT` for `Y`, and `AGE`
   for `Color`. After checking this you should see a nice scatter plot for `WEIGHT` and `HEIGHT`
   with the color corresponding to `AGE`:
   ![exercises-scripting-viewer](exercises-scripting-viewer.png)
7. In the property panel, proceed to modify the value of the "Script" field by clicking on a "..."
   icon in the text field.
8. The Python code you see is what renders the scatter plot form p.6 on the Datagrok server. Let's walkthrough this
   code.
    * The script takes as inputs the original dataframe and the three columns. Remember form p.6 there were selectors
      for `X`, `Y`, and `Color` in the property panel. In fact, these three property names are declared with the
      notation `<propertyName>ColumnName` in the names of the three `#input` columns.
    * The script produces an `#output` of type `graphics`. It is important the graphics appear in the end of the Python
      script. This is exactly what happens with the `plt.show()` in the last line of the script.

9. Modify the name of `colorColumnName` to a `temperatureColumnName`, hit `Apply` in the bottom of the window, and check
   what happens to the `Color` field in the property panel.
10. Add another input parameter to the script with a name `Title`. Hit `Apply` and check what appears in the property
    panel.
11. Add another input column to the script with a name `SEX`. Hit `Apply` and check what appears in the property panel.
12. Now there's all you need to create a Python scripting viewer for our amino acid histogram task. Open a demo file
    with nucleotide sequences. It is located at `Data | Files | Demo Files | bio | sars-cov-2.csv`.
    `Data` corresponds to the first button from the top on the Datagrok sidebar.
13. In the top menu you've activated at p.2, hit `Add | Scripting Viewers | New Scripting Viewer`.
14. Follow what you've learned in the points 1 to 11 to create a scripting viewer taking a column of strings, expecting
    to have nucleotide sequences in them, and plotting a [Matplotlib's histogram][016] with all amino acid triplets
    occurred within all of these sequences. As you may notice, `numpy` and `matplotlib` are already available for your
    Python scripting in Datagrok. Reuse them to finish this exercise.

## Exercise 5: Transforming dataframes

*Prerequisites:* exercises ["Setting up the environment"](#setting-up-the-environment),
["Semantic types"](#exercise-1-semantic-types).

*You will learn:* how to join and union dataframes using the knowledge of semantic types, and display the result.

1. Make sure the [prerequisites](#setting-up-the-environment) are prepared on your machine, including the package
   called `<yourFirstName>-sequence` Assure the package carries a relevant semantic type detector from the
   exercise ["Semantic Types"](#exercise-1-semantic-types).
2. Add a function to the package as follows:

   ```javascript
   //name: fuzzyJoin
   //input: dataframe df1
   //input: dataframe df2
   //input: int N
   ...
   ```

3. Implement a `fuzzyJoin` function which takes two dataframes `df1` and `df2`, and does the following:

    * takes a first column in `df1` which has a semantic type of `dna_nucleotide`, let's say it is `col1`
    * takes a first column in `df2` which has a semantic type of `dna_nucleotide`, let's say it is `col2`
    * creates a dataframe `df` out of `df1` and `df2` in the following way:
      * the content of `df2` goes after `df1`, and all columns of `df1` and `df2` are preserved — this is a UNION operation
        for dataframes, as in SQL; use the dataframe's [`.append`](https://public.datagrok.ai/js/samples/data-frame/append)
        method
      * a new column `Counts` appears in `df`, which contains:
        * for each row `R` from `df1`, `R.counts` is a number of matches of all the subsequences in `R.col1` of length `N`
          in *all* the sequences of `col2`
        * symmetrically, same for each row from `df2` — consider this as a fuzzy, programmatic JOIN of the two dataframes;
          use[`df.columns.addNew`](https://public.datagrok.ai/js/samples/data-frame/modification/manipulate)
          , [`col.set(i, value)`](https://public.datagrok.ai/js/samples/data-frame/advanced/data-frames-in-columns)
          on a newly created column
    * displays `df` with [`grok.shell.addTableView`](https://public.datagrok.ai/js/samples/data-frame/test-tables)

4. Deploy the package with `webpack` and `grok publish dev`. Unlike with the first exercise, where the package was built
   on the Datagrok server, in this one we locally build the package before sending it. In addition, webpack output helps
   find some syntactic errors in JavaScript.
5. Launch the platform, open the two files from `"Demo files"`: `sars-cov-2.csv` and `a-h1n1.csv`, and run the
   package's `fuzzyJoin` function using one of the methods you've learned. The result for N=3 should look similar to:
   ![exercises-transforming-dataframes](exercises-transforming-dataframes.png)
6. Read more about joining dataframes through the case reviewed at our
   [Community Forum](https://community.datagrok.ai/t/table-to-table-augmentation/493/4), and with
   [a sample](https://public.datagrok.ai/js/samples/data-frame/join-link/join-tables).

<!--- TODO: add linked dataframes demo here --->

## Exercise 6: Custom cell renderers with 3-rd party js libraries

*You will learn:* reuse 3-rd party JavaScript libraries in your Datagrok packages; render cells by semantic types.

*Prerequisites:* exercises ["Setting up the environment"](#setting-up-the-environment),
["Semantic types"](#exercise-1-semantic-types).

1. Navigate into the folder with your `<yourFirstName>-sequence` package created in
   ["Setting up the environment"](#setting-up-the-environment).
2. Let's add a custom cell renderer for a *nucleotide sequence box* to represent our sequences in high density on the
   screen. We need to render each nucleotide sequence with a monospace font in small letter sizing, fitting into a
   rectangular cell area and adding ellipsis to the end of the string if it won't fit. This is a basis for a very useful
   nucleotide sequence representation in bioscience applications. Let's use a 3-rd party JavaScript
   library `fusioncharts-smartlabel` to compute the text fit. Add it to your package by navigating in its folder and
   calling:
   `npm install fusioncharts-smartlabel --save`
   The `--save` key updates `package.json` to add this library to your package dependencies.
3. Add a class to `src/package.js` for the new cell renderer:

    * use `fusioncharts-smartlabel` to break the original sequence in the current cell into lines which fit into a cell's
      canvas rectangle; learn [here][017] how to do it, consider `SmartLabel.textToLines(...).lines`
      as a target array of lines to render
    * Datagrok [grid](../../visualize/viewers/grid.md) is rendered through an
      [HTML5 Canvas](https://en.wikipedia.org/wiki/Canvas_element). The grid's canvas is `g.canvas`. Iterate through the
      resulting lines and bring them to a `g.canvas` in the `render` method with `g.canvas.getContext("2d").fillText`; learn
      more about HTML Canvas if it's new for you
    * Hint: pay attention to managing `line-height` both at computing the box and rendering text lines

    ```javascript
    class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
      get name() { return 'Nucleotide cell renderer'; }
      get cellType() { return 'dna_nucleotide'; }
      render(g, x, y, w, h, gridCell, cellStyle) {
        let seq = gridCell.cell.value;
        const sl = new SmartLabel('id', true);
        sl.setStyle({/* ... */});
        // ...
        let ctx = g.canvas.getContext("2d");
        ctx.font = '11px courier';
        // ...
        const lines = labelObj.lines;
        for (let i = 0; i < lines.length; i++)
          ctx.fillText(/* ... */);
      }
    }
    ```

4. Add the below to `src/package.js` to make the new cell renderer part of the package:

   ```javascript
    //name: nucleotideBoxCellRenderer
    //tags: cellRenderer
    //meta.cellType: dna_nucleotide
    //output: grid_cell_renderer result
    export function nucleotideBoxCellRenderer() {
      return new NucleotideBoxCellRenderer();
    }
    ```

5. Deploy the package as usual with `grok publish dev`. In [Datagrok](https://public.datagrok.ai), navigate to a file
   with nucleotide sequences from `"Demo files"`, such as `sars-cov-2.csv`. Verify you get the desired result, it should
   look similar to this:
   ![exercises-custom-cell-renderer](exercises-custom-cell-renderer.png)
   Change the "Sequence" column width and rows heights with a mouse to see how things adujst.
6. (*) Implement a colored nucleotide sequence box where backgrounds of `A`, `G`, `C`, `T` vary. Choose one of the
   popular coloring conventions, following [this link](https://www.biostars.org/p/171056/).

## Exercise 7: Accessing web services with OpenAPI

*Details:* [OpenAPI access](../../access/open-api.md)

Web services often provide their API specs in an [OpenAPI (Swagger)](../../access/open-api.md) format in a JSON or a
yaml file. Because OpenAPI spec file is standardized, the API may now be directly loaded and later queried. Datagrok
provides for connecting to API data sources and fetching API querying results as dataframes. In this lesson we will
connect to the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/) and fetch some nucleotide data regarding
coronavirus.

1. Obtain a ENA's Swagger JSON file for the [ENA Browser](https://www.ebi.ac.uk/ena/browser),
   following [this link](https://www.ebi.ac.uk/ena/browser/api/v2/api-docs). It usually takes you some effort to reach
   the JSON Swagger at the API Browser link. The rule of thumb there is to use the browser's "Network" tab in the
   Developers Console (available by F12), and identify there a resource called "
   api-docs". Usually, when there's an [API Navigator][018], there's also a JSON Swagger. It's also possible to
   understand the API through its [Swagger API tester][019]. Follow recommendations [here][020]. In particular, modify
   the connection's `Name` to `ENA`, `Url` to
   `https://www.ebi.ac.uk/ena/browser/api/`. Save this file to a desktop with a name, say, `ENA.json`.
2. Load the Swagger into Datarok by drag-and-dropping it into the platform window.
3. Check the connection is valid with the `Test` button, and hit `Ok` to save the edit.
4. In the expanded view of the `ENA` connection, locate `Perform a text search and download data in XML format`
   and hit `Run` or double-click it.
5. Enter the parameter values: set `Query` to `coronavirus`, `Result` to `assembly`. Hit `Ok`. As a result, you'd find a
   table, which was prepared from the received XML file by Datagrok.
6. Close the table, locate the saved query in the list and run it.
7. Bring the connection to the package:

    * Put the Swagger file in a `swaggers` folder of the package. Make sure the swagger contains the correct `basePath`
      and `host`, in some Swaggers it isn't always the case.
    * Add the following function to `package.js`:

    ```
    //name: testENASwagger
    export async function testENASwagger() {
      let data = await grok.data.query('<yourFirstName>sequence:PerformATextSearchAndDownloadDataInXMLFormat',
        {'query': 'coronavirus', 'result': 'assembly'});
      grok.shell.addTableView(data);
    }
    ```

    * Note how the Swagger's query name translates into a package query name.
    * You can obtain this query name with the Datagrok UI. Click on the query of interest,
      `"Perform a text search and download data in XML format"` in our case, and find a `Links...`
      section. Click it and copy a function's name from the URI.
    * Deploy the package and make sure `testENASwagger` function works in Datagrok.

We provide a handful of demo Swaggers, check their source JSON files [here][021] and see in action in Datagrok at
the [`Web Services`](https://public.datagrok.ai/webservices) section of the Datagrok UI.

## Exercise 8: Creating an info panel with a REST web service

We will use the ENA REST API to output sequences and associated data in the info panel, based on the ENA sequence ID
contained in a currently selected grid cell.

1. Searching through [the ENA archive](https://www.ebi.ac.uk/ena/browser/text-search?query=coronavirus), you may notice
   the sequences' IDs have a format of `[A-Z]{2}[0-9]{6}` (two capital letters + six digits). Go to
   the [detectors file](#exercise-1-semantic-types) of your package and add a detector which recognizes a string of this
   form:

   ```javascript
   //input: string str
   //output: bool result
   isPotentialENAId(str) {
     // returns true, if name is of the form [A-Z]{2}[0-9]{6}
   }
   ```

2. Use [`fetchProxy`](../how-to/access-data.md#rest-endpoints) to get a sequence for the potential corresponding ENA ID
   in fasta format. For example, this GET fetches the sequence for the `ID=AA046425`:
   [`https://www.ebi.ac.uk/ena/browser/api/fasta/AA046425`](https://www.ebi.ac.uk/ena/browser/api/fasta/AA046425)
   Use the following structure for the into panel function in your `src/package.js`:

   ```javascript
    //name: ENA Sequence
    //tags: panel, widgets
    //input: string cellText {semType: ENA}
    //output: widget result
    //condition: isPotentialENAId(cellText)
    export async function enaSequence(cellText) {
      const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${cellText}`;
      const fasta = await (await grok.dapi.fetchProxy(url)).text();
      return new DG.Widget(ui.box(
        // ... the widget controls are composed here
      ));
    }
   ```

   Incorporate
   a [`textInput`](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/components/accordion.js)
   control to display a sequence in a scrollable fashion. Add a caption to that text area to display an ENA's name for
   this sequence, which also comes in the fasta file. Use a
   [`splitV`](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/layouts/splitters.js)
   control to nicely locate the caption at the top and the text area at the bottom.

`fetchProxy` mimics the regular `fetch` method of ECMAScript, but solves
a [CORS](https://developer.mozilla.org/en-US/docs/Web/HTTP/CORS) limitation of JavaScript. In this panel, you'd query
the external domain from your web page, whereas CORS prevents you from querying anything outside a reach of your web
page's domain. Thus Datagrok provides a proxy facility in the neat `fetchProxy` wrapper.

## Exercise 9: Enhancing Datagrok with dialog-based functions

In the previous exercises we've learned how the Datagrok function inputs are offered in a dialog window automatically
once you run the function. In this exercise we find how to expand these dialogs with the behaviour beyond simple
arguments-to-inputs mapping.

So, in some previous exercises we've used the files `a-h1n1.csv` and `sars-cov-2.csv` which we prepared in advance.
These files contain ENA sequence ID along with the first 60 letters of the sequence by ID. Let's construct a
dialog-based function which forms such files automatically by a given search input. The search topic may
be `coronavirus`, `influenza` etc.

1. This `GET` query performs a text search in the EMBL database, returning a `limit` first results (`10` in this case):
   `https://www.ebi.ac.uk/ena/browser/api/embl/textsearch?result=sequence&query=coronavirus&limit=10`
   By the way, you could discover this API endpoint via a Swagger API navigator [at this link][018]. Let's assume the
   result we want is always of type `sequence`. Create a function `_fetchENASequence` which takes as parameters
   a `query` ad a `limit` and returns a dataframe with two string columns `ID` and `Sequence`. Use this structure for
   dataframe construction:

   ```javascript
    df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', [ /* a list of IDs you've parsed from a ENA output */ ]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', [ /* corresponding list of sequences */ ])
    ]);
    ```

   The output from `ebi.ac.uk` is a raw text, and you need to parse it to get the desired pieces. Trim the sequence so
   that isn't longer than 60 characters. Use your previous knowledge about [`fetchProxy`][022] to do the `GET` query.

2. Make a function `formENADataTable` which constructs a dialog giving the user a two-step process for constructing a
   dataframe with ENA sequence data in it.

    * First, the user can type in the query (`coronavirus` is the default setting) and see the first 10 results in the grid
      right in this window after clicking the "Search" button. Consider this as a preview before the actual dataframe is
      produced.
    * Second, when the user is happy with what's in the preview, he/she proceeds to the "Ok" button to get the actual
      dataframe with the ENA data on the screen in the Datagrok's grid view. This table shall consist of the number of rows
      the user chooses (`100` set as a default).

    Here is the code scaffold for the `formENADataTable` function:

    ```javascript
    let grid = DG.Viewer.grid(df);
    let limitInput = ui.intInput('How many rows: ', 100);
    let queryInput = ui.stringInput('Query: ', 'coronavirus');
    let button = ui.button('Preview');
    ui.dialog('Create sequences table')
      .add(ui.splitV([
        ui.splitH([
          ui.span([queryInput.root]),
          button
        ]),
        ui.div([grid]),
        ui.div([limitInput])
      ]))
      .onOK(() => {
        /* Handle table creation */
        // Display the resulting table
        grok.shell.addTableView(df);
      })
      .show();
    ```

    Re-use twice the `_fetchENASequence` function you've prepared previously.

3. In this first version we fetched `60` characters for a sequence. Add a new text field called `Sequece length`
   to let the user specify this trim length, set it `60` as a default.

4. Make your function set a proper [semantic type](#exercise-1-semantic-types) for the `Sequence` column.

5. (*) You may notice the sequences you get in this order are not too different. Add more diversity to these tables. For
   example, you can use the `offset` parameter of the `GET` query.

<!---
Search for a keyword to form a table with limits
https://www.ebi.ac.uk/ena/browser/api/

## Persisting user sessions and tables

Saving the search parameters

## Creating an application

A simple keyword search in the ENA database (with navigation)
--->

[014]: ../how-to/add-info-panel.md#functions "How to add an info panel"

[015]: https://youtu.be/p7_qOU_IzLM?t=724 "Dev Meeting 1: Getting Started – First-class functions"

[`CountSubsequencePython`]: #scripting-with-server-functions "Scripting with server functions"

[016]: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist.html "Matplotlib Histogram"

[017]: https://medium.com/@priyanjit.dey/text-wrapping-and-ellipsis-overflow-a-platform-independent-solution-30fb737ff609 "Text Wrapping and Ellipsis"

[018]: https://www.ebi.ac.uk/ena/browser/api "ENA API Navigator"

<!--- TODO: Verify 019 --->

[019]: https://www.ebi.ac.uk/ena/browser/api "Swagger API Tester"

[020]: ../../access/open-api.md#troubleshooting "OpenAPI connections troubleshooting"

[021]: https://github.com/datagrok-ai/public/tree/master/packages/Swaggers/swaggers "Datagrok Swaggers samples"

[022]: #exercise-8-creating-an-info-panel-with-a-rest-web-service "Creating an info panel with a REST web service"
