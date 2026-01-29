// Creating a DataFrame from a CSV string and checking the formats identified for the columns

let table = DG.DataFrame.fromCsv(
  `country, continent, year,  lifeExp,  pop,       gdpPercap
Angola,   Africa,    1992,  40.647,   8735988,   2627.8457031
Austria,  Europe,    1992,  76.04,    7914969,   27042.0195313
Austria,  Europe,    1992,  76.04,    7914969,   27042.0195313
Bolivia,  Americas,  1992,  59.957,   6893451,   2961.699707
Bulgaria, Europe,    1992,  71.19,    8658506,   6302.6235352
Poland,   Europe,    1992,  70.99,    38370697,  7738.8813477
Zambia,   Africa,    1992,  46.1,     8381163,   1210.8846436
`);

let view = grok.shell.addTableView(table);
grok.shell.info(`Format detected on 'gdpPercap': ${table.getCol('gdpPercap').meta.format}`);
grok.shell.info(`Format detected on 'lifeExp': ${table.getCol('lifeExp').meta.format}`);
grok.shell.info(`Max significant digits detected for 'lifeExp': ${table.getCol('lifeExp').meta.sourcePrecision}`);
