//api: DG.ViewBase.saveToProject, DG.ViewBase.saveToProjectByDefault
// A view that can reopen from its state is offered in the Save project dialog, left out until the user
// includes it. Include one from code, as the tab's "Save with project" item does.

const tv = grok.shell.addTableView(grok.data.demo.demog());
await grok.functions.call('CmdNewMarkupView');
const notes = grok.shell.v;
notes.name = 'Notes';
grok.shell.info(`Saved without asking: ${notes.saveToProjectByDefault}`);
notes.saveToProject = true;
grok.shell.v = tv;
grok.shell.info(`Saved with the project: ${notes.saveToProject}`);
