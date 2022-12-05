/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import * as types from './types';
import * as constants from './const';
import * as alationApi from './alation-api';
import * as utils from './utils';
import {getBaseURL} from './package';

let isSettingDescription = false;

export function createTree(
  objects: types.baseEntity[], objectType: types.specialType, treeRootNode?: DG.TreeViewGroup): DG.TreeViewGroup {
  objects = utils.filterDuplicates(objects);
  treeRootNode ??= ui.tree();
  const iconClass = objectType === 'data-source' ? 'svg-data-connection' :
    objectType == 'query' ? 'svg-data-query' : 'svg-database-tables';
  const iconElement = `<i class="grok-icon svg-icon ${iconClass}"></i>`;

  if (objectType == 'query') {
    (objects as types.query[]).forEach((queryObject) => {
      const name = queryObject.title.trim() || `Unnamed query id ${queryObject.id}`;
      const item = treeRootNode!.item(name, queryObject);
      $(item.root).children().first().before(iconElement);

      item.root.addEventListener('dblclick', async () => {
        const progressIndicator = DG.TaskBarProgressIndicator.create('Opening table...');
        try {
          await connectToDb(queryObject.datasource_id, async (c: DG.DataConnection) => runQuery(c, queryObject));
        } catch {
          grok.shell.error('Couldn\'t retrieve table');
        }
        progressIndicator.close();
      });
      item.root.addEventListener('mousedown', async (ev) => {
        if (ev.button != 0 || isSettingDescription)
          return;

        isSettingDescription = true;
        const description = queryObject.description
          .replaceAll('src="/', `src="${await getBaseURL()}`)
          .replaceAll('href="/', `href="${await getBaseURL()}`);
        $('.alation-description').empty().html(description);
        isSettingDescription = false;
      });
    });

    return treeRootNode;
  }

  if (objectType === 'table') {
    (objects as types.table[]).forEach((tableObject) => {
      const name =
        (tableObject.title || tableObject.name).trim() || `Unnamed table id ${tableObject.id}`;
      const item = treeRootNode!.item(name, tableObject);
      $(item.root).children().first().before(iconElement);
      item.root.addEventListener('dblclick', async () => {
        const progressIndicator = DG.TaskBarProgressIndicator.create('Opening table...');
        try {
          await connectToDb(tableObject.ds_id, async (c: DG.DataConnection) => getTable(c, tableObject));
        } catch {
          grok.shell.error('Couldn\'t retrieve table');
        }
        progressIndicator.close();
      });
      item.root.addEventListener('mousedown', async (ev) => {
        if (ev.button != 0 || isSettingDescription)
          return;

        isSettingDescription = true;
        const description = tableObject.description
          .replaceAll('src="/', `src="${await getBaseURL()}`)
          .replaceAll('href="/', `href="${await getBaseURL()}`);
        $('.alation-description').empty().html(description);
        isSettingDescription = false;
      });
    });

    return treeRootNode;
  }

  objects.forEach(async (dataSourceObject) => {
    const currentId = dataSourceObject.id;
    const name =
      (dataSourceObject.title || (dataSourceObject as (types.schema | types.table)).name).trim() ||
      `Unnamed ${objectType} id: ${currentId}`;
    const group = treeRootNode!.group(name, dataSourceObject, false);
    isFirstTimeMap[objectType] ??= {};
    isFirstTimeMap[objectType][currentId] = true;

    $(group.root).children().first().children().first().after(iconElement);

    group.root.addEventListener('mousedown', async (ev) => {
      if (ev.button != 0 || isSettingDescription)
        return;
      isSettingDescription = true;
      const pi = DG.TaskBarProgressIndicator.create('Loading child entities...');
      const description = dataSourceObject.description
        .replaceAll('src="/', `src="${await getBaseURL()}`)
        .replaceAll('href="/', `href="${await getBaseURL()}`);
      $('.alation-description').empty().html(description);
      if (isFirstTimeMap[objectType][currentId]) {
        isFirstTimeMap[objectType][currentId] = false;
        await getChildren(objectType, currentId, group);
      }
      pi.close();
      isSettingDescription = false;
    });
  });

  if (objectType == 'schema') {
    const queryGroup = treeRootNode.group('Queries', null, false);
    const currentId = (objects[0] as types.schema).ds_id;

    isFirstTimeMap[objectType] ??= {};
    isFirstTimeMap[objectType][`query ${currentId}`] = true;

    $(queryGroup.root).children().first().children().first().after(`<i class="grok-icon svg-icon svg-data-query"></i>`);

    queryGroup.root.addEventListener('mousedown', async (ev) => {
      if (ev.button != 0 || isSettingDescription)
        return;

      isSettingDescription = true;
      const pi = DG.TaskBarProgressIndicator.create('Loading child entities...');

      if (isFirstTimeMap[objectType][`query ${currentId}`]) {
        isFirstTimeMap[objectType][`query ${currentId}`] = false;
        await getChildren('query', currentId, queryGroup);
      }
      pi.close();
      isSettingDescription = false;
    });
  }

  return treeRootNode;
}

const isFirstTimeMap: {[key: string]: {[key: string]: boolean}} = {};

async function getChildren(
  objectType: types.specialType, currentId: number, group: DG.TreeViewGroup): Promise<DG.TreeViewGroup> {
  let dataList: types.baseEntity[];
  let nextObjectType: types.specialType;

  switch (objectType) {
  case 'data-source':
    dataList = await alationApi.getSchemas(currentId);
    nextObjectType = 'schema';
    break;
  case 'schema':
    dataList = await alationApi.getTables(currentId);
    nextObjectType = 'table';
    break;
  case 'query':
    dataList = await alationApi.getQueries(currentId);
    nextObjectType = objectType;
    break;
  default:
    throw new Error(`Unknown datasource type '${objectType}'`);
  }
  return createTree(dataList, nextObjectType, group);
}

export async function runQuery(conn: DG.DataConnection, queryObject: types.query): Promise<void> {
  try {
    let query = conn.query(queryObject.title, queryObject.content);
    query = await grok.dapi.queries.save(query);

    const df = await query.executeTable();
    df.name = queryObject.title;
    grok.shell.addTableView(df);
  } catch (e) {
    grok.shell.info('Could not run query. See console');
    console.error(e);
  }
}

export async function connectToDb(
  dataSrouceId: number, handler: (conn: DG.DataConnection) => Promise<void>): Promise<void> {
  const dataSource = await alationApi.getDataSourceById(dataSrouceId);
  const connName = dataSource.title || dataSource.qualified_name || dataSource.dbname;
  const thisUser = await grok.dapi.users.current();
  const dsConnections = (await grok.dapi.connections.filter(connName).list()).filter((c) => c.author.id == thisUser.id);

  for (const conn of dsConnections) {
    if (await grok.dapi.permissions.check(conn, "Edit"))
      return handler(conn);
  }

  connectToDbDialog(dataSource, handler);
}

function connectToDbDialog(dataSource: types.dataSource, func: (conn: DG.DataConnection) => Promise<void>) {
  let dbType: string | null = null;
  for (const dsType of constants.DATA_SOURCE_TYPES) {
    if (dataSource.dbtype === dsType.toLowerCase()) {
      dbType = dsType;
      break;
    }
  }

  if (dbType === null)
    throw new Error(`DBTypeError: Unsupported DB type '${dataSource.dbtype}'`);

  const helpText = ui.inlineText(['The database credentials are stored in the secure ',
    ui.link('Datagrok Credentials Management Service', 'https://datagrok.ai/help/govern/security#credentials'),
    ' and the connection is created that can be used as any other ',
    ui.link('Data Connection', 'https://datagrok.ai/help/access/data-connection'), ' in Datagrok.']);
  const helpHost = ui.div(helpText, 'alation-help-host');
  $(helpHost).width(500);
  const usernameField = ui.stringInput('Login', '');
  const passwordField = ui.stringInput('Password', '');
  $(passwordField.root as HTMLInputElement).children('.ui-input-editor').attr('type', 'password');

  const connectionName = dataSource.title || dataSource.qualified_name || dataSource.dbname;
  const dialog = ui.dialog(`Connect to ${connectionName}`)
    .add(ui.divV([helpHost, usernameField, passwordField]))
    .onOK(async () => {
      const progress = DG.TaskBarProgressIndicator.create(`Connecting to ${connectionName}...`);
      try {
        const dcParams = {
          dataSource: dbType!,
          server: `${dataSource.host}:${dataSource.port}`,
          db: dataSource.dbname,
          login: usernameField.stringValue,
          password: passwordField.value,
        };
        let dsConnection = DG.DataConnection.create(connectionName, dcParams);
        dsConnection = await grok.dapi.connections.save(dsConnection);
        await func(dsConnection);
      } catch (e) {
        grok.shell.info('Could not connect to the data source. See console');
        console.error(e);
      }
      progress.close();
    })
    .show();
  return dialog;
}

export async function getTable(dsConnection: DG.DataConnection, tableObject: types.table): Promise<void> {
  try {
    const query = DG.TableQuery.create(dsConnection);
    query.table = `${tableObject.schema_name}.${tableObject.name}`;
    const columns = await alationApi.getColumns(tableObject.id);
    query.fields = columns.map((v) => v.name);
    const df = await query.executeTable();
    df.name = tableObject.name || `Unnamed table id ${tableObject.id}`;
    const tableView = grok.shell.addTableView(df);
    tableView.name = tableObject.title ?? df.name;
  } catch (e) {
    grok.shell.info('Could not retreive table. See console');
    console.error(e);
  }
}
