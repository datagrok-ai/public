// An example of using Dart's future as JavaScript's promises
grok.data.openTable("System:AppData/Samples/cars.csv")
  .then(t => grok.shell.info(t.name));
