/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import * as alationApi from './alation-api';
import * as utils from './utils';
import * as view from './view';

export const _package = new DG.Package();
let baseUrl: string;

export async function getBaseURL() {
  const properties = await _package.getProperties() as {[key: string]: any};
  baseUrl = properties['Base URL'] as string;
  if (!baseUrl)
    throw new Error('PackagePropertyError: Base URL is not set!');
  baseUrl = baseUrl.endsWith('/') ? baseUrl : `${baseUrl}/`;
  return baseUrl;
}

//name: Alation
//tags: app
//top-menu: Admin | Alation @Toolbox Data | Alation
export async function Alation() {
  const progressIndicator = DG.TaskBarProgressIndicator.create('Loading Alation...');

  const treeHost = ui.div();
  const descriptionHost = ui.div(undefined, 'alation-description');
  const rightPanelHost = ui.divV([ui.h1('Data Sources'), treeHost]);
  descriptionHost.style.maxWidth = '50%';
  descriptionHost.style.margin = '10px';
  grok.shell.newView('Alation Browser', [ui.divH([rightPanelHost, descriptionHost])]);

  await utils.retrieveKeys();

  const dataSourcesList = await alationApi.getDataSources();
  const tree = view.createTree(dataSourcesList, 'data-source');
  $(treeHost).append([tree.root]);

  grok.events.onContextMenu.subscribe((args: any) => {
    const tableObject = args.args.item.value;
    if (typeof tableObject.table_type == 'undefined')
      return;
    const name = (tableObject.title || tableObject.name).trim() || `Unnamed table id ${tableObject.id}`;
    const contextMenu = args.args.menu;
    contextMenu.item('Open table', async () => view.connectToDb(tableObject, name));
  });

  progressIndicator.close();
}
