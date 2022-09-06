const app = (appName: string) => `
The application ${appName} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/build-an-app
See application examples at https://public.datagrok.ai/apps`;

const connection = (connectionName: string) => `
The connection ${connectionName} has been added successfully
Read more at https://datagrok.ai/help/access/data-connection,
https://datagrok.ai/help/develop/how-to/access-data#connections
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Chembl`;

const detector = (semType: string) => `
The detector for ${semType} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors`;

const func = (funcName: string, isPanel: boolean = false) => {
  return `\nThe function ${funcName} has been added successfully\n` +
  'Read more at https://datagrok.ai/help/datagrok/functions/function' +
  (isPanel ? ',\nhttps://datagrok.ai/help/develop/how-to/add-info-panel\n' : '\n') +
  'See examples at https://public.datagrok.ai/functions' +
  (isPanel ? ',\nhttps://public.datagrok.ai/js/samples/functions/info-panels/info-panels' : '');
}

const query = (queryName: string) => `
The query ${queryName} has been added successfully
Read more at https://datagrok.ai/help/access/data-query,
https://datagrok.ai/help/develop/how-to/access-data#queries
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Chembl`;

const script = (scriptName: string, packageName: string) => `
The script ${scriptName} has been created. To call it from a JavaScript file, use:

  await grok.functions.call('${packageName}:${scriptName}', { params });
  
Read more at https://datagrok.ai/help/compute/scripting
See examples at https://public.datagrok.ai/scripts,
https://public.datagrok.ai/js/samples/scripting/scripting`;

const view = (viewName: string) => `
The view ${viewName} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/custom-views
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Notebooks`;

const viewer = (viewerName: string) => `
The viewer ${viewerName} has been added successfully
Read more at https://datagrok.ai/help/develop/how-to/develop-custom-viewer
See examples at https://github.com/datagrok-ai/public/tree/master/packages/Viewers,
https://public.datagrok.ai/js/samples/functions/custom-viewers/viewers`;

const _package = (ts: boolean) =>
ts ? '' : 'Consider TypeScript as a language for package development, to start over, run `grok create` with the `--ts` flag\n' +
'Likely next steps: `grok add` to add functionality, `grok publish` to upload the package';

export const help = {
  app,
  connection,
  detector,
  func,
  query,
  script,
  view,
  viewer,
  package: _package,
};
