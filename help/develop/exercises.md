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
  * [Scripting](#scripting)
  * [Querying databases](#querying-databases)
  * [Creating a scripting viewer](#creating-a-scripting-viewer)

## Setting up the environment

Prerequisites: basic JavaScript knowledge

1. Install the necessary tools (Node.js, npm, webpack, datagrok-tools) following [these instructions](develop.md#getting-started)
2. Get a dev key for https://dev.datagrok.ai (you will work with this server) and add it by running `grok config`
3. Create a default package called `<name>-sequence` using datagrok-tools: `grok create <name>-sequence`
4. Upload it to the server: `grok publish dev --rebuild` (see other options [here](develop.md#deployment-modes))
5. Launch the platform and run the package's `test` function using different methods: 
    * via the [Functions](https://dev.datagrok.ai/functions?q=test) view
    * via the [Packages](https://dev.datagrok.ai/packages?) menu (find your package and run `test` from the 'Content' pane in the property panel)
    * via the console: press `~` and type `<Name>Sequence:test()` (note that the function's namespace corresponds to the package name)

## Semantic types

Prerequisites: basic JavaScript knowledge

Details: [How to Create a Semantic Type Detector](how-to/semantic-type-detector.md),
[How to Add an Info Panel](how-to/add-info-panel.md)

You will learn: how to write semantic type detectors, how to develop context-specific data augmentation.  

1. Create a `complement` function that takes a nucleotide string and returns its complement.
   Essentially, change each character to the complementary one: A<=>T, G<=>C. 
   Run it and check whether everything works fine. 
2. Now, let's specify that this function is meant to accept not any string, but nucleotides only,
   and returns a nucleotide string as well. In order to do that, let's annotate both input and output parameters
   with the `dna_nucleotide` semantic type (add `{semType: dna_nucleotide}` after the parameter name).
   At this point, `dna_nucleotide` string does not have any meaning, but we will connect the dots later.
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
   and change the output type to `widget` (see an example [here](how-to/add-info-panel.md#functions)).
   This will instruct the platform to use the `complement` function for providing additional information for string values
   of the `dna_nucleotide` semantic type. To test it, simply open our test file, click on any cell
   in the `sequence` column, and find the `complement` property in the panel on the right.
   
## Scripting

Prerequisites: basic Python knowledge

Details: [Scripting](develop/scripting.md), [Dev Meeting 1 | First-class functions](https://youtu.be/p7_qOU_IzLM?t=724)

_You will learn:_ how to create and invoke Datagrok scripts in data science languages like R and Python.

In this exercise, we will count occurrences of a given subsequence in a nucleotide sequence, using Python.

1. Open Datagrok and navigate to `Functions | Scripts | New Python Script`.
2. Observe a default script created for you. All script attributes are specified in the beginning in comments.
   There we have the script name, language, one input value of type [`dataframe`](),
   and one output value of type `int`. The script simply computes number of cells in the dataframe.  
   [Dataframe](develop/how-to/build-an-app.md) is a high-performance, easy to use tabular structure
   with strongly-typed columns of different types (supported types are: `string`, `bool`, `int`, `bigint`,
   `double`, `qnum` and `datetime`). In this exercise, we only see a dataframe as is in the default script;
   there is [another exercise]() to learn manipulating dataframes in JavaScript.
3. Run the script to get a hint for creating an input file. An attribute `#sample: cars.csv`
   is responsible for it. To open a default input file `cars`, click the `Star` icon in the top menu.
4. Run the script again and proceed to the Datagrok's console. As in Quake, it's available
   by pressing a `~` button anywhere inside Datagrok. In the console, you would see the script
   execution result. Just one line above the result you could see the console's command to execute
   the script. Enter it again to the console to get the same result.
5. Let's modify the script to solve the task of counting sequence occurrences. Add a new preamble:
   (use any `#description` you like):
   ```python
   #name: CountSubsequencePython
   #language: python
   #input: string sequence
   #input: string subsequence
   #output: int count
   ```
   In the body, implement a Python function counting all occurrences of a given `subsequence` in a `sequence`.
   Return a `count` the same way as in the default script from p. 2.
6. Run the script function, provide input values in the dialog and get to the console to see the result.
   Now run the script function again through the console completely, passing different arguments values.
7. Go back to `Functions | Scripts` and hit `New JavaScript Script`.
8. Implement the function which does the same as `CountSubsequencePython`, now in JavaScript, named
   `CountSubsequenceJS`. Follow the same conventions on the parameters in the comments block
   and returning a result via a variable.
9. Run `CountSubsequenceJS` using the `Play` button; using the console.
   You can notice that both Python and JS versions of our function are homogenous functions
   in Datagrok. It's also possible to call them in a uniform fashion
   [using our JavaScript API](develop/scripting.md#running-a-script).
10. Don't forget to save these two scripts. We would re-use them in other exercises.

The difference between the two scripts is that the first, `CountSubsequencePython`, runs on
our server by a [compute virtual machine](develop/admin/architecture.md#compute-virtual-machine),
whereas the second, `CountSubsequenceJS`, runs directly in the browser. To run `CountSubsequencePython`,
Datagrok passes the script arguments over the network and fetches back the result to the browser.

## Querying databases

Prerequisites: basic SQL knowledge

Details: [Connecting to Databases](https://www.youtube.com/watch?v=dKrCk38A1m8&t=1048s),
[How to Access Data](how-to/access-data.md)

_Note:_ Editing an existing data query requires the respective access permission. You might need to request one.

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

## Creating a Scripting Viewer

Prerequisites: basic Python knowledge, [matplotlib](https://matplotlib.org/) or a similar library

Details: [Scripting](develop/scripting.md), [Scripting Viewer](visualize/viewers/scripting-viewer.md),
[Creating a scripting viewer (video)](https://www.youtube.com/embed/jHRpOnhBAz4).

*Amino acids counting task.* In this exercise, we'd use a Python script to generate a histogram
(a distribution plot) for amino acids occurring in a column of nucleotide sequences. Amino acids are simply triples
of nucleotides from which the nucleotide DNA sequence is made of. These are also called triplets, or codon-codes.
As there are 4 letters `G`, `A`, `T`, `C`, there are 4 to power 3 protein amino acids: `GTA`, `AGC`, `TTT`, and so forth.

We don't know at which starting point each nucleotide sequence was cut: it could either be a junction of two triplets,
or one-third or two-third of a triplet. Therefore, we'd count in our statistics for all three possible cuts, starting
the reading frame off at offsets 0, 1, and 2 from the beginning of the nucleotide sequence.

Say, we are given a sequence `TTTAATTACAGACCTGAA`. We start to count triplets _without overlap_ from an offset 0 first,
getting: `TTT`, `AAT`, `TAC`, `...`, `GAA`. Then we move to the offset 1, getting: `TTA`, `...`, `CTG`. Lastly,
we move to the offset 2 and get: `TAA`, `...`, `TGA`. In the histogram we'd count for all these triplets.

First, let's explore how scripting viewer works.

1. Open a `demog` demo file with demographic data. It is located at `Data | Files | Demo Files | demog.csv`.
`Data` corresponds to the first button from the top on the activity bar at the left of the Datagrok window.
Make sure the table view with the data appears.
2. Activate the top menu from the activity bar, using a `Windows | Menu` switch.
3. In this menu, hit `Add | Scripting Viewers | Python | Scatter Plot`.
4. See that the viewer appeared on the right, telling though it is "Unable to plot with current settings".
5. Proceed to the viewer properties by hitting on the gear icon in the viewer's title.
6. Make sure the chosen values for "Data" are `HEIGHT` for `X`, `WEIGHT` for `Y`, and `AGE` for `Color`.
After checking this you should see a nice scatter plot for `WEIGHT` and `HEIGHT` with the color corresponding to `AGE`:
![](exercises-scripting-viewer.png)
7. In the property panel, proceed to modify the value of the "Script" field by clicking on a "..." icon in the text field.
8. The Python code you see is what renders the scatter plot form p.6 on the Datagrok server. Let's walk through this code.
   * The script takes as inputs the original dataframe and the three columns. Remember form p.6 there were
     selectors for `X`, `Y`, and `Color` in the property panel. In fact, these three property names are
     declared with the notation `<NAME>ColumnName` in the names of the three `#input` columns.
   * The script produces an `#output` of type `graphics`. It is important the graphics appear in the end
     of the Python script. This is exactly what happens with the `plt.show()` in the last line of the script.
9. Modify the name of `colorColumnName` to a `temperatureColumnName`, hit `Apply` in the bottom of the window,
   and check what happens to the `Color` field in the property panel.
10. Add another input parameter to the script with a name `Title`. Hit `Apply` and check what appears in the 
    property panel.
11. Add another input column to the script with a name `SEX`. Hit `Apply` and check what appears in the property panel.
12. Now there's all you need to create a Python scripting viewer for our amino acid histogram task.
    Open a demo file with nucleotide sequences. It is located at `Data | Files | Demo Files | bio | sequences.csv`.
`Data` corresponds to the first button from the top on the activity bar at the left of the Datagrok window.
13. In the top menu you've activated at p.2, hit `Add | Scripting Viewers | New Scripting Viewer`.
14. Follow what you've learned in the points 1 to 11 to create a scripting viewer taking a column of strings,
    expecting to have nucleotide sequences in them, and plotting a Matplotlib's [histogram](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist.html) with all amino acid triplets
    occurred within all of these sequences.  
    As you may notice, `numpy` and `matplotlib` are already available for your Python scripting in Datagrok.
Reuse them to finish this exercise.