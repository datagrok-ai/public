1. Run the JS script
```JS
//name: Link Test
//language: javascript
//input: dataframe df1 [Data table 1] 
//input: dataframe df2 [Data table 2]
//input: dataframe df3 [Data table 3]

grok.data.linkTables(df3, df2, ['Sample Name', 'link column 1', 'link column 2', 'link column 3'], ['Sample Name', 'link column 1', 'link column 2', 'link column 3'], [DG.SYNC_TYPE.FILTER_TO_FILTER]);
grok.data.linkTables(df1, df2, ['Id'], ['Concept Id'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);
grok.shell.addTableView(df1);
grok.shell.addTableView(df2); 
grok.shell.addTableView(df3); 
```
1. Open:
   - `Df1 = SPGI`
   - `Df2 = SPGI-linked1`
   - `Df3 = SPGI-linked2`
2. In **SPGI**, select several rows (from the beginning of the table)
3. Switch to **SPGI-linked2**, open the **Filter Panel**, and apply some filters
4. Switch to **SPGI-linked1** and apply filters using its **Filter Panel**

**Expected Result**: **SPGI-linked1** should display only the rows that meet all of the following conditions:
- selected in **SPGI**,
- filtered through **SPGI-linked2**,
- filtered again using the **SPGI-linked1** Filter Panel