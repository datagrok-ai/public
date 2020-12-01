// Aggregations
let demog = grok.data.demo.demog();
demog.selection.init((i) => i % 2);

let agesByRace = demog
    .groupBy(['race', 'sex'])
    .whereRowMask(demog.selection)
    .avg('age')
    .add('med', 'age')
    .aggregate();

grok.shell.addTableView(agesByRace);