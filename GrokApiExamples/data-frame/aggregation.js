// Aggregations

let demog = gr.testData('demog', 5000);

let avgAgesByRace = demog
    .groupBy(['race', 'sex'])
    .avg('age')
    .aggregate();

gr.addTableView(avgAgesByRace);