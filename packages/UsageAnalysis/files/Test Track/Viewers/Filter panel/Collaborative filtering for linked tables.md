1. Run the JS script
```JS
//name: Link Test
//language: javascript

let df1 = await grok.data.files.openTable('System:DemoFiles/SPGI.csv');
let df2 = await grok.data.files.openTable('System:DemoFiles/SPGI-linked1.csv');
let df3 = await grok.data.files.openTable('System:DemoFiles/SPGI-linked2.csv');

grok.data.linkTables(df3, df2, ['Sample Name', 'link column 1', 'link column 2', 'link column 3'], ['Sample Name', 'link column 1', 'link column 2', 'link column 3'], [DG.SYNC_TYPE.FILTER_TO_FILTER]);
grok.data.linkTables(df1, df2, ['Id'], ['Concept Id'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);
grok.shell.addTableView(df1);
grok.shell.addTableView(df2);
grok.shell.addTableView(df3);
 
```
2. Go to the **SPGI** view, select 5 rows on the top of the table
2. Switch to **SPGI-linked1** view - it should contain 9 filtered rows
3. Switch to **SPGI-linked2** view
3. Open the **Filter Panel**
3. For the `link column 3` filter, select the `v ii` category
3. Switch to **SPGI-linked1** view - it should contain 5 filtered rows
4. Open the **Filter Panel** 
4. For the `PAMPA Classification` filter, select the `Inconclusive` category - the number of the filtered rows should be 2

