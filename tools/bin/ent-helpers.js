const app = (appName) => `
The application ${appName} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/build-an-app
See application examples at https://public.datagrok.ai/apps`;

const connection = (connectionName) => `
The connection ${connectionName} has been added successfully
Read more at https://datagrok.ai/help/access/data-connection,
https://datagrok.ai/help/develop/how-to/access-data#connections
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Chembl`;

const detector = (semType) => `
The detector for ${semType} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors`;

const func = (funcName, isPanel = false) => {
  return `\nThe function ${funcName} has been added successfully\n` +
  'Read more at https://datagrok.ai/help/overview/functions/function' +
  (isPanel ? ',\nhttps://datagrok.ai/help/develop/how-to/add-info-panel\n' : '\n') +
  'See examples at https://public.datagrok.ai/functions' +
  (isPanel ? ',\nhttps://public.datagrok.ai/js/samples/functions/info-panels/info-panels' : '');
}

const query = (queryName) => `
The query ${queryName} has been added successfully
Read more at https://datagrok.ai/help/access/data-query,
https://datagrok.ai/help/develop/how-to/access-data#queries
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Chembl`;

const script = (scriptName, packageName) => `
The script ${scriptName} has been created. To call it from a JavaScript file, use:

  await grok.functions.call('${packageName}:${scriptName}', { params });
  
Read more at https://datagrok.ai/help/develop/scripting
See examples at https://public.datagrok.ai/scripts,
https://public.datagrok.ai/js/samples/scripting/scripting`;

const view = (viewName) => `
The view ${viewName} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/custom-views
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Notebooks`;

const viewer = (viewerName) => `
The viewer ${viewerName} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/develop-custom-viewer
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Viewers,
https://public.datagrok.ai/js/samples/functions/custom-viewers/viewers`;

const package = (name, ts = false) => `
Successfully created package \`${name}\`${ts ? '' : '\nConsider TypeScript as a language for package development, to start over, run `grok create` with the `--ts` flag'}
Likely next steps: \`grok add\` to add functionality, \`grok publish\` to upload the package`;

module.exports = {
  app,
  connection,
  detector,
  func,
  query,
  script,
  view,
  viewer,
  package,
};
