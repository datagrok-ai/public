// Row filtering

let agesByRace = grok.data.demo.demog()
  .groupBy(['race', 'sex'])
  .where('race = Asian and started after 2/2/2000')
  .avg('age')
  .avg('started')
  .add('med', 'age')
  .aggregate();

grok.shell.addTableView(agesByRace);