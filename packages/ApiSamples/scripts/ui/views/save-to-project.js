//api: DG.ViewBase.saveToProject, DG.ViewBase.saveToProjectByDefault
// Views that can reopen from their state are saved with the project by default.
// Keep a scratch view out of the project: the project then holds only the table and its view.

const tv = grok.shell.addTableView(grok.data.demo.demog());
await grok.functions.call('CmdNewMarkupView');
const notes = grok.shell.v;
notes.name = 'Scratch notes';
grok.shell.info(`Saved by default: ${notes.saveToProjectByDefault}`);
notes.saveToProject = false;
grok.shell.v = tv;
grok.shell.info(`In the project: ${grok.shell.project.children.map((c) => c.name).join(', ')}`);
