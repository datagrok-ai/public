// Aggregations

let demog = grok.testData('demog', 5000);

let avgAgesByRace = demog
    .groupBy(['race', 'sex'])
    .avg('age')
    .aggregate();

grok.addTableView(avgAgesByRace);