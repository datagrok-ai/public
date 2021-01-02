// Row filtering

let agesByRace = grok.data.demo.demog()
  .groupBy(['race', 'sex'])
  .where('race = Asian and started after 1/1/2003')
  .avg('age')
  .min('started')
  .add('med', 'age')
  .aggregate();

grok.shell.addTableView(agesByRace);