//api: ui.pickFileFromFiles
// The platform's 'Select a file' dialog over the file shares: double-click a file, or select one and press OK.
// The promise resolves with its FileInfo (fullPath names it on the share), or null when the dialog is closed without a pick.
// {multiple: true} takes every file selected with Ctrl+click or Shift+click; {connection} browses one share alone.

const file = await ui.pickFileFromFiles();
if (file !== null)
  grok.shell.info(`${file.fullPath}: ${(await file.readAsString()).length} characters`);
else
  grok.shell.info('No file picked');

const demo = await grok.dapi.connections.filter('shortName = "DemoFiles"').first();
const files = await ui.pickFileFromFiles({multiple: true, connection: demo});
grok.shell.info(files === null ? 'No files picked' : files.map((f) => f.fullPath).join(', '));
