// File browser widget
grok.shell.newView('Files', [ui.fileBrowser().root]);

// Accepts a full-qualified name (see [Entity.nqName]) in the path parameter
// and opens a directory in the files tree, if it is specified, e.g., 'System:AppData/Chem/tests'
const testConnection = await grok.dapi.connections.filter('shortName = "Home"').first();
const fb = ui.fileBrowser({path: testConnection.nqName});
grok.shell.newView(`${testConnection.name} files`, [fb.root]);
