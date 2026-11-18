// Doing what a double-click on an object does in the platform: a file opens, a project opens, a script opens
// its editor. Nothing happens for an object the platform has no default action for.

const [file] = (await grok.dapi.files.list('System:DemoFiles/', false, 'csv'))
  .filter((f) => f.fileName === 'demog.csv');
DG.ObjectHandler.runDefaultAction(file);
