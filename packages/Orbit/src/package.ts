/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Orbit
//tags: app
export async function app() {
  // ?c.dataSource=S3&c.parameters.bucket=datagrok-orbit&c.parameters.region=us-east-2&folder=202208269918

  let search: string[] = [];
  const params = new URLSearchParams(window.location.search);
  params.forEach((value, key) => {
    if (key.startsWith('c.')) {
      search.push(`${key.substring(2)} = "${value}"`);
    }
  });
  // find a root connection by params. datasource, bucket, region
  let connections = await grok.dapi.connections.filter(`shares.connection = null and ${search.join(' and ')}`).list();

  if (connections.length > 1) {
    grok.shell.error('Ambiguous connection parameters');
    return;
  }
  if (connections.length == 0) {
    grok.shell.error('Connection not found');
    return;
  }
  let connection = connections[0];
  let folder = params.get('folder');

  // share the subfolder
  let newConnection = await grok.dapi.connections.shareFolder(connection, folder!)
  console.log(newConnection);

  // open files view to subfolder with shared connection filter
  let view = DG.FilesView.create({'t': `id = "${newConnection.id}"`, 'path': newConnection.nqName});
  grok.shell.addView(view);
}
