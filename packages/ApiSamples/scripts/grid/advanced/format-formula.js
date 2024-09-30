// Ability to render values from other columns (useful for hyperlinks, etc)
// Not that other than formatting, the column would still
// work just like the regular column when sorting, editing, etc

const df = grok.data.demo.demog(5);
df.col('race').tags['%formatFormula'] = '${sex} text [${age}](https://google.com/search?q=${race})';

df.plot.grid.root