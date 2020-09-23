// Aggregations
let demog = grok.data.demo.demog();
demog.selection.init((i) => i % 2);

let avgAgesByRace = demog
    .groupBy(['race', 'sex'])
    .whereRowMask(demog.selection)
    .avg('age')
    .aggregate();

grok.shell.addTableView(avgAgesByRace);