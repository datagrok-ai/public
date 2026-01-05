// Shows groups of DataFrames by race and sex

let demog = grok.data.demo.demog();

let agesByRaceGroups = demog
  .groupBy(['race', 'sex'])
  .getGroups();

Object.keys(agesByRaceGroups).forEach((key, ind) => grok.shell.addTableView(agesByRaceGroups[key]));