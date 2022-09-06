/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import * as alationApi from './alation-api';
import * as utils from './utils';
import * as view from './view';
import * as types from './types';

export const _package = new DG.Package();
let _p: DG.Package;

export function setPackage(p: DG.Package): void {
  _p = p;
}

export function getPackage(): DG.Package {
  return _p ?? _package;
}

export async function getBaseURL() {
  const properties = await getPackage().getProperties() as {[key: string]: any};
  let baseUrl = properties['Base URL'] as string;
  if (!baseUrl)
    throw new Error('PackagePropertyError: Base URL is not set!');
  baseUrl = baseUrl.endsWith('/') ? baseUrl : `${baseUrl}/`;
  return baseUrl;
}

export async function getUserGroup() {
  const properties = await getPackage().getProperties() as {[key: string]: any};
  const userGroupName = properties['User group'] as string;
  if (!userGroupName)
    throw new Error('PackagePropertyError: User group is not set!');
  return userGroupName;
}

//name: Alation
//tags: app
//top-menu: Admin | Alation @Toolbox Data | Alation
export async function Alation() {
  const progressIndicator = DG.TaskBarProgressIndicator.create('Loading Alation...');

  const treeHost = ui.box();

  const descriptionHost = ui.panel(undefined, 'alation-description');
  const title = ui.box(ui.h1('Data Sources'), {style: {maxHeight: '30px', margin: '5px'}});
  const rightPanelHost = ui.splitV([title, treeHost], {style:{maxWidth: '300px'}});
  
  const host = ui.splitH([rightPanelHost, descriptionHost]);
  
  let v = grok.shell.newView('Alation Browser', [host]);
  v.box = true;
  
  await utils.retrieveKeys();
  
  const dataSourcesList = await alationApi.getDataSources();
  const tree = view.createTree(dataSourcesList, 'data-source');
  treeHost.append(tree.root);

  grok.events.onContextMenu.subscribe((args: any) => {
    const obj = args.args.item.value as types.table & types.query;
    const contextMenu = args.args.menu;
    if (obj.ds_id && obj.table_type) {
      contextMenu.item('Open table', async () => {
        const progressIndicator = DG.TaskBarProgressIndicator.create('Opening table...');
        try {
          await view.connectToDb(obj.ds_id, async (conn: DG.DataConnection) => {await view.getTable(conn, obj);});
        } catch {
          grok.shell.error('Couldn\'t retrieve table');
        }
        progressIndicator.close();
      });
    } else if (obj.datasource_id && obj.content) {
      contextMenu.item('Run query', async () => {
        const progressIndicator = DG.TaskBarProgressIndicator.create('Running query...');
        try {
          await view.connectToDb(obj.datasource_id,
            async (conn: DG.DataConnection) => {await view.runQuery(conn, obj);});
        } catch {
          grok.shell.error('Couldn\'t retrieve table');
        }
        progressIndicator.close();
      });
    }
  });

  progressIndicator.close();
}
