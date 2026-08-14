// Publish computed tables as a dashboard via the platform's standard
// Save-project dialog — without putting them in the workspace first.
// Data sync defaults on for tables carrying a creation script (.script tag);
// pass the TableView showing a table to ship its layout with the project.

const df = grok.data.demo.demog(100);
df.name = 'My demog snapshot';

// Optional: a (possibly detached) view whose layout should reopen with the table.
const tv = DG.TableView.create(df, false);

const project = await DG.Project.showSaveDialog({
  tables: [df],
  views: [tv],           // views[i] shows tables[i]; pass null to skip
  // layouts: [state],   // a ViewLayout.viewState string per table — ships a layout without a live view
  // project: savedProjectOrId,  // re-publish INTO a previously saved project (updates, not a new one)
  name: 'My dashboard',
  description: 'Published from a script',
});

grok.shell.info(project == null ? 'Cancelled' : `Saved: ${project.friendlyName}`);
