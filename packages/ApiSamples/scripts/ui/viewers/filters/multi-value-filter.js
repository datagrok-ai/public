// Filters automatically switch to the multi-value mode
// when the DG.TAGS.MULTI_VALUE_SEPARATOR tag is set on a column

let languages = `Country,Languages
Belgium,"Dutch
French
German"
Burundi,"French
Kirundi
English"
Cameroon, "English
French"
Canada, "English
French"
Italy, "Italian
French
German
Slovene
Corsican"
Rwanda, "English
French"
Netherlands, Dutch
Switzerland, "French
German
Italian"`;

let t = DG.DataFrame.fromCsv(languages);
t.col('languages').meta.multiValueSeparator = '\n';

grok.shell.addTableView(t).addViewer(DG.VIEWER.FILTERS);