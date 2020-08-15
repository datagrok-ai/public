// Aggregations
let demog = grok.data.demo.demog();

let avgAgesByRace = demog
    .groupBy(['race', 'sex'])
    .avg('age')
    .aggregate();

grok.shell.addTableView(avgAgesByRace);