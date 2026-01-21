/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';
import {DBExplorerConfig} from '@datagrok-libraries/db-explorer/src/types';
import {DBExplorer} from '@datagrok-libraries/db-explorer/src/db-explorer';
import {DBExplorerEditor} from '@datagrok-libraries/db-explorer/src/editor';
import {DB_EXPLORER_OBJ_HANDLER_TYPE, DBExplorerObjectHandler} from '@datagrok-libraries/db-explorer/src/object-handlers';


/** Sets up editing of the db explorers config(s) */
export async function setupGlobalDBExplorer() {
  let savedConfigs: DBExplorerConfig[] | null = null;
  const deferredTreeNodes: DG.TreeViewGroup[] = [];
  const processedNodes = new Set<DG.TreeViewGroup>();

  const processConnectionNode = (node: DG.TreeViewGroup, highlightOnAdd: Boolean = false) => {
    if (processedNodes.has(node))
      return;
    processedNodes.add(node);
    const connection = DG.toJs(node.value) as DG.DataConnection;
    const matchingConfigs = savedConfigs?.filter((c) =>
      c.nqName === connection.nqName && c.dataSourceName === connection.dataSource);
    matchingConfigs?.forEach(async (config) => {
      const hasPermission = await grok.dapi.permissions.check(connection, 'Edit');
      if (!hasPermission)
        return; // remove in future according to what is decided
      const itemElement = ui.divH([ui.iconFA('fingerprint', ()=>{}), ui.divText(`Identifiers - ${config.schemaName}`)],
        {style: {gap: '6px', alignItems: 'center'}});
      ui.tooltip.bind(itemElement, `Configure identifiers for ${config.schemaName} schema...`);
      const tnItem = node.item(itemElement, new DBExplorerConfigWrapper(config));
      tnItem.onSelected.subscribe(async () => {
        let editor: DBExplorerEditor | null = null;
        // no connection section needs to be rendered here
        const view = DG.View.create();
        view.name = `${connection.name}.${config.schemaName} Identifiers`;
        grok.shell.addPreview(view);
        ui.setUpdateIndicator(view.root, true, 'Loading editor...');
        editor = new DBExplorerEditor(() => onSaveAction(editor!), false);
        const editorUI = await editor.getUI();
        view.root.appendChild(editorUI);
        ui.setUpdateIndicator(view.root, true, 'Initializing configuration...');
        await editor.setConfig(config);
        ui.setUpdateIndicator(view.root, false);
      });
      tnItem.root.style.order = '-1'; // make sure it is on top
      if (highlightOnAdd)
        ui.hints.addHintIndicator(tnItem.root, false, 4000);
    });
  };

  grok.shell.browsePanel.mainTree.onChildNodeExpanding.subscribe((node) => { // this event fires once per node
    // if the node in question is a dataConnection node, we add our stuff.
    if (DG.toJs(node?.value) instanceof DG.DataConnection) {
      if (savedConfigs == null)
        deferredTreeNodes.push(node);
      else
        processConnectionNode(node);
    }
  });

  // on context menu for processing db-explorer nodes
  grok.events.onContextMenu.subscribe((args) => {
    if (!(DG.toJs(args?.args?.item) instanceof DG.TreeViewNode) ||
    !(DG.toJs(args.args.item.value) instanceof DBExplorerConfigWrapper) ||
      !args?.args?.menu)
      return;
    const node = DG.toJs(args.args.item) as DG.TreeViewNode<DBExplorerConfigWrapper>;
    const config = node.value.config;
    const menu = args.args.menu as DG.Menu;
    menu.item('Edit...', async () => {
      (DG.toJs(args.args.item) as DG.TreeViewNode).root.click(); // trigger selection, good enaugh
    });
    menu.item('Remove...', async () => {
      ui.dialog('Remove Identifier Configuration')
        .add(ui.divText(`Are you sure you want to remove ${config.connectionName}.${config.schemaName} identifier configuration?`))
        .onOK(async () => {
          const oldConfigs = await getSavedConfigs();
          const newConfigs = oldConfigs.filter((c) =>
            !(c.nqName === config.nqName &&
                    c.dataSourceName === config.dataSourceName && c.schemaName === config.schemaName));
          try {
            await _package.files.writeAsText('db-explorer/configs.json', JSON.stringify(newConfigs, null, 2));
            grok.shell.info('Identifier configuration removed. Refresh the page to see the changes.');
            savedConfigs = newConfigs;
            node.remove();
          } catch (e) {
            _package.logger.error('Failed to save db-explorer configs:' + e);
            grok.shell.error('Failed to remove identifier configuration:');
          }
        }).show();
    }, 4, {description: 'Remove this identifier configuration from the system'});
  });


  // on context menu for processing connection nodes
  grok.events.onContextMenu.subscribe((args) => {
    if (!(DG.toJs(args?.args?.item) instanceof DG.TreeViewGroup) ||
    !(DG.toJs(args.args.item.value) instanceof DG.DataConnection) ||
      !args?.args?.menu ||
      (DG.toJs(args.args.item.value) as DG.DataConnection).dataSource === DG.DataSourceType.Files ||
      DG.DataSourceType.fileDataSources.includes((DG.toJs(args.args.item.value) as DG.DataConnection).dataSource))
      return;
    const connection = DG.toJs(args.args.item.value) as DG.DataConnection;
    const menu = args.args.menu as DG.Menu;
    menu.item('Configure Identifiers...', async () => {
      // ask for schema selection first in a dialog
      const schemas = await grok.dapi.connections.getSchemas(connection);
      if (schemas.length === 0) {
        grok.shell.error('No schemas found for this connection.');
        return;
      }
      const defaultSchema = schemas.includes('public') ? 'public' : schemas[0];

      const schemaChoice = ui.input.choice('Schema', {value: defaultSchema, items: schemas,
        tooltipText: 'Select primary schema for identifiers', nullable: false} );

      ui.dialog('Select primary schema for Identifiers Configuration')
        .add(ui.form([schemaChoice]))
        .onOK(async () => {
          const view = DG.View.create();
          view.name = `${connection.name}.${schemaChoice.value!} Identifiers`;
          ui.setUpdateIndicator(view.root, true, 'Loading editor...');

          grok.shell.addView(view);
          let editor: DBExplorerEditor | null = null;
          editor = new DBExplorerEditor(() => {
            onSaveAction(editor!);
          }, false); // no connection section needs to be rendered here

          const uiEl = await editor.getUI();
          ui.setUpdateIndicator(view.root, true, 'Initializing configuration...');
          // find if there is an existing config for this connection-schema
          const existingConfig = savedConfigs?.find((c) =>
            c.nqName === connection.nqName && c.dataSourceName === connection.dataSource &&
            c.schemaName === schemaChoice.value!);
          if (existingConfig) {
            await editor.setConfig(existingConfig);
            grok.shell.info('Identifiers configuration for this connection and schema already exists.\n Importing...');
          } else {
            await editor.setConfig({
              connectionName: connection.name,
              nqName: connection.nqName,
              dataSourceName: connection.dataSource,
              schemaName: schemaChoice.value ?? undefined,
            }, false);
          } // no warning because schema is deliberately left empty for user to fill in
          ui.setUpdateIndicator(view.root, false);
          view.root.appendChild(uiEl);
        }).show();
    }, 4, {description: 'Configure identifier settings for this database connection'});
  });


  async function getSavedConfigs(): Promise<DBExplorerConfig[]> {
    let configs: DBExplorerConfig[] = [];
    try {
      if (!await _package.files.exists('db-explorer/configs.json'))
        return configs;
      const saved = await _package.files.readAsText('db-explorer/configs.json');
      if (saved)
        configs = JSON.parse(saved);
    } catch (e) {
      _package.logger.error('Failed to read saved db-explorer configs:' + e);
    }
    return configs;
  }


  // we still might be late, so we need to process the nodes
  setTimeout(() => {
    const databasesNode = grok.shell.browsePanel.mainTree.children.find((n) => n.text === 'Databases');
    if (databasesNode && databasesNode instanceof DG.TreeViewGroup && databasesNode.expanded) {
      databasesNode.children.forEach((providersNode) => {
        // next level is providers
        if (providersNode instanceof DG.TreeViewGroup && providersNode.expanded) {
          providersNode.children.forEach((connectionNode) => {
            // this level is datasources
            if (connectionNode instanceof DG.TreeViewGroup && connectionNode.expanded &&
                DG.toJs(connectionNode.value) instanceof DG.DataConnection)
              processConnectionNode(connectionNode);
          });
        }
      });
    }
  }, 1000);


  async function onSaveAction(editor: DBExplorerEditor) {
    // validation is already happening in the editor before calling onSave
    const oldConfigs = await getSavedConfigs();
    const newConfig = editor.getConfig();
    const existingConfiIndex = oldConfigs.findIndex((c) =>
      c.nqName === newConfig.nqName &&
          c.dataSourceName === newConfig.dataSourceName && c.connectionName === newConfig.connectionName &&
          c.schemaName === newConfig.schemaName);
    if (existingConfiIndex >= 0) // one config per connection
      oldConfigs[existingConfiIndex] = newConfig;
    else
      oldConfigs.push(newConfig);
    try {
      // warn the user if there is already a db-explorer for this connection-schema.
      // only one per connection-schema is allowed
      if (existingConfiIndex === -1) {
        await _package.files.writeAsText('db-explorer/configs.json', JSON.stringify(oldConfigs, null, 2));
        grok.shell.info('Identifier configuration saved. Refresh the page to see the changes.');
        savedConfigs = oldConfigs;
        // if we are here, then it is a new config, so we need to add it to the tree
        const connectionNode = Array.from(processedNodes)
          .find((t) => (DG.toJs(t.value) as DG.DataConnection).nqName === newConfig.nqName &&
          (DG.toJs(t.value) as DG.DataConnection).dataSource === newConfig.dataSourceName);
        if (connectionNode) {
          // remove it from the set and re-process
          processedNodes.delete(connectionNode);
          processConnectionNode(connectionNode, true);
        }
      } else {
        ui.dialog('Overwrite existing configuration?')
          .add(ui.divText('A DB Explorer configuration for this connection and schema already exists. \n' +
                  'Do you want to overwrite it?'))
          .onOK(async () => {
            await _package.files.writeAsText('db-explorer/configs.json', JSON.stringify(oldConfigs, null, 2));
            grok.shell.info('Identifier configuration saved. Refresh the page to see the changes.');
            savedConfigs = oldConfigs;
          })
          .show();
      }
    } catch (e) {
      _package.logger.error('Failed to save db-explorer configs:' + e);
      grok.shell.error('Failed to save identifier configurations:');
    }
  }

  savedConfigs = await getSavedConfigs();
  // after this, we process any deferred nodes
  deferredTreeNodes.forEach((node) => processConnectionNode(node));
  savedConfigs.forEach((config) => {
    DBExplorer.initFromConfig(config);
  });
}

class DBExplorerConfigWrapper {
  constructor(public config: DBExplorerConfig) {}
}

/** when querying, table produced can have information in it that says where the table come from.
 * for these purposes, we will setup a generic handler
 */
export async function setupDBQueryCellHandler() {
  const explorePanelRoot = ui.div([], {style: {minHeight: '30px', width: '100%'}});
  let processTimer: any = null;
  const conInfoCache = new Map<string, Promise<{name: string, nqName: string, friendlyName: string} | null>>();
  const getConnectionInfo = async (connectionId: string) => {
    if (conInfoCache.has(connectionId))
      return conInfoCache.get(connectionId)!;
    conInfoCache.set(connectionId, (async () => {
      try {
        const connection = await grok.dapi.connections.find(connectionId);
        if (connection)
          return {name: connection.name, nqName: connection.nqName, friendlyName: connection.friendlyName ?? connection.name};
      } catch (e) {
        console.error('Failed to get connection info for id ' + connectionId, e);
      }
      return null;
    })());
    return conInfoCache.get(connectionId)!;
  };

  const isValidCurrentObject = (cell: DG.Cell) => {
    return (cell?.dart && (cell.rowIndex ?? -1) >= 0 && cell.column?.dart &&
      cell.dataFrame?.dart) && (grok.shell.o == cell.column || // cell equivalency does not work
        (grok.shell.o instanceof DG.SemanticValue && grok.shell.o.cell?.dart && grok.shell.o.cell.column === cell.column && grok.shell.o.cell.rowIndex === cell.rowIndex));
  };

  function processCell(cell: DG.Cell) {
    ui.empty(explorePanelRoot);
    explorePanelRoot.appendChild(ui.divText('No data to explore.'));
    ui.setUpdateIndicator(explorePanelRoot, true);
    if (processTimer != null)
      clearTimeout(processTimer);
    processTimer = setTimeout(async () => {
      if (!isValidCurrentObject(cell)) { // no need to process anything
        ui.setUpdateIndicator(explorePanelRoot, false);
        return;
      }
      // collect necessary information from the cell
      const cId = cell.dataFrame.tags.get(DG.Tags.DataConnectionId);
      const schemaName = cell.column.tags.get(DG.Tags.DbSchema);
      const tableName = cell.column.tags.get(DG.Tags.DbTable);
      const columnName = cell.column.tags.get(DG.Tags.DbColumn);
      if (!cId || !schemaName || !tableName || !columnName) {
        ui.setUpdateIndicator(explorePanelRoot, false);
        return;
      }
      const conInfo = await getConnectionInfo(cId);
      if (!conInfo || !conInfo?.nqName || !conInfo?.name) {
        ui.setUpdateIndicator(explorePanelRoot, false);
        return;
      }
      const handler = DG.ObjectHandler.list().find((h) => { // for entity method does not work, because of instanceof checks.
        // need to do it manually// underneeth, its doing the same with listing all handlers
        return h.type === DB_EXPLORER_OBJ_HANDLER_TYPE && (h as DBExplorerObjectHandler).connectionNqName === conInfo.nqName! && (h as DBExplorerObjectHandler).schemaName === schemaName;
      });
      if (handler) {
        const acc = (handler as DBExplorerObjectHandler).renderPropertiesFromDfRow(cell.row, schemaName, tableName);
        ui.empty(explorePanelRoot);
        explorePanelRoot.appendChild(acc);
      } else {
        // if there is no handler, at this point we can just register it
        console.log('Registering generic handler for ', conInfo.nqName, schemaName);
        const exp = DBExplorer.initFromConfig({
          connectionName: conInfo.name!,
          nqName: conInfo.nqName!,
          schemaName: schemaName,
          joinOptions: [],
          entryPoints: {},
        });
        const newHandler = exp.genericValueHandler; // nex time it will be found
        const acc = newHandler.renderPropertiesFromDfRow(cell.row, schemaName, tableName);
        ui.empty(explorePanelRoot);
        explorePanelRoot.appendChild(acc);
      }
      ui.setUpdateIndicator(explorePanelRoot, false);
    }, 300);
  }

  grok.events.onAccordionConstructed.subscribe((acc) => {
    if (!(acc.context instanceof DG.Column) && !(acc.context instanceof DG.SemanticValue))
      return;
    const getCol = () => {
      if (acc.context instanceof DG.Column)
        return acc.context as DG.Column;
      const sv = acc.context as DG.SemanticValue;
      if (sv.cell && sv.cell.dart && sv.cell.column)
        return sv.cell.column;
      return null;
    };

    const col = getCol();
    if (!col)
      return;
    const df = col?.dart ? col.dataFrame : null;
    if (!df || !df.tags.get(DG.Tags.DataConnectionId) || !col.tags.get(DG.Tags.DbTable) ||
        !col.tags.get(DG.Tags.DbColumn))
      return;
    ui.empty(explorePanelRoot);
    explorePanelRoot.appendChild(ui.divText('Select a cell to explore its value...'));
    const pane = acc.addPane('Explore', () => {
      const children = [getEnrichDiv(col)];
      if (col.tags.get(DG.Tags.DbSchema))
        children.push(explorePanelRoot);
      return ui.divV([children], {style: {width: '100%'}});
    });
    getConnectionInfo(df.tags.get(DG.Tags.DataConnectionId)!).then((conInfo) => {
      if (conInfo) // this way is better than making subscription async...
        pane.name = conInfo.friendlyName;
    });
    // if there is a current cell and its in the given column, process it
    if (col.tags.get(DG.Tags.DbSchema) && col.dataFrame.currentCell && col.dataFrame.currentCell.column === col)
      processCell(col.dataFrame.currentCell!);
  });

  grok.events.onCurrentCellChanged.subscribe((cell?: DG.Cell) => {
    if (!cell || !isValidCurrentObject(cell))
      return;
    const col = cell.column;
    const connectionID = col.dataFrame.tags.get(DG.Tags.DataConnectionId)!;
    const schemaName = col.tags.get(DG.Tags.DbSchema)!;
    const tableName = col.tags.get(DG.Tags.DbTable)!;
    const columnName = col.tags.get(DG.Tags.DbColumn)!;
    if (!connectionID || !schemaName || !tableName || !columnName)
      return;
    processCell(cell);
  });
}

function getEnrichDiv(col: DG.Column): HTMLElement {
  const enrichAcc = ui.accordion('Enrich Column');
  enrichAcc.addPane('Enrich', () => {
    if ([DG.COLUMN_TYPE.BOOL, DG.COLUMN_TYPE.DATE_TIME].includes(col.type as DG.COLUMN_TYPE)) {
      return ui.info(`Cannot use column of ${col.type} type for enrichment. 
Supported types are ${[DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BIG_INT, DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT].join(', ')}.`);
    }

    const df = col.dart ? col.dataFrame : null;
    if (!df || df.rowCount === 0)
      return ui.info('Data frame is empty or not available.');

    const schemaName = col.tags.get(DG.Tags.DbSchema);
    const connId: string = df.tags.get(DG.Tags.DataConnectionNqName) ?? df.tags.get(DG.Tags.DataConnectionId);
    const dbColName = col.tags.get(DG.Tags.DbColumn);
    const tableName = col.tags.get(DG.Tags.DbTable);

    if (!connId || !dbColName || !tableName)
      return ui.info('Column is not linked to database table.');

    return ui.wait(async () => {
      const conn = connId.includes(':')
          ? (await grok.dapi.connections.filter(`nqName = "${connId}"`).list()).find((_) => true)
          : await grok.dapi.connections.find(connId);
      if (!conn)
        return ui.info('Failed to find connection for this column.');
      const tables: DG.TableInfo[] = await grok.dapi.connections.getSchema(conn, schemaName, tableName);
      if (tables.length === 0)
        return ui.info('Could not find a main table used for SQL query. Please, try to specify full "<schema>.<table>" name in SQL query.');
      if (tables.length > 1)
        return ui.info('Ambiguous table name — specify full "<schema>.<table>" name in SQL query.');

      const mainTable: DG.TableInfo = tables[0];
      const enrichmentsRoot = getEnrichmentsDiv(conn, mainTable.tags.get(DG.Tags.TableSchema), mainTable.friendlyName, dbColName, df);
      const addEnrichBtn = document.createElement('button');
      addEnrichBtn.append(ui.icons.add(() => {}), ui.span(['Add enrichment']));
      addEnrichBtn.classList.add('power-pack-enrich-add');
      addEnrichBtn.onclick = () => showEnrichDialog(mainTable, df, dbColName, () => {
        enrichmentsRoot.replaceWith(getEnrichmentsDiv(conn, mainTable.tags.get(DG.Tags.TableSchema), mainTable.friendlyName, dbColName, df));
      });
      return ui.div([
        enrichmentsRoot,
        addEnrichBtn
      ]);
    });
  });
  return ui.divV([enrichAcc.root]); // wraping in div for consistent styling
}

function getEnrichmentsDiv(conn: DG.DataConnection, schema: string, table: string, column: string, df: DG.DataFrame) {
  return ui.wait(async () => {
    const names = await getEnrichConfigsNames(conn.nqName, schema, table, column);
    const empty = ui.div([
      ui.divText('No enrichments yet.', 'power-pack-enrich-empty-title'),
      ui.divText('Add an enrichment for the column to extend this table.', 'power-pack-enrich-empty-sub')
    ], 'power-pack-enrich-empty');
    if (names.length === 0)
      return empty;
    let root: HTMLDivElement | null = null;
    root = ui.divV(names.map((n) => {
      let parent: HTMLDivElement | undefined = undefined;
      const deleteLink = ui.iconFA('times', () => {
        deleteEnrichment(conn.nqName, schema, table, column, n).then((b) => {
          if (b)
            parent?.remove();
          if (root?.children?.length === 0)
            root.replaceWith(empty);
        });
      }, 'Delete');

      const editLink = ui.iconFA('pencil', async () => {
        const enrichment = await readEnrichConfig(conn.nqName, schema, table, column, n);
        if (!enrichment) {
          grok.shell.error('Something went wrong when opening enrichment. Please try again.');
          return;
        }
        const tables: DG.TableInfo[] = await grok.dapi.connections.getSchema(conn, schema, table);
        if (tables.length === 0 || tables.length > 1) {
          grok.shell.error('Could not find key table or table name is ambiguous. Please check the connection and database schema and try again.');
          return;
        }
        const tq = DG.TableQuery.fromTable(tables[0]).build();
        tq.fields = enrichment.fields;
        tq.joins = enrichment.joins;
        showEnrichDialog(tables[0], df, column, () => {
          root?.replaceWith(getEnrichmentsDiv(conn, schema, table, column, df));
        }, n, tq);
      }, 'Edit');

      const runLink = ui.link(n, () => {
        const enrichFunc = DG.Func.find({package: 'PowerPack', name: 'runEnrichment'})[0];
        const progress = DG.TaskBarProgressIndicator.create('Enriching...');
        enrichFunc.prepare({'conn': conn, 'schema': schema, 'table': table, 'column': column, 'name': n, 'df': df})
          .call(false, progress, {report: true, processed: false}).finally(() => progress.close());
      }, 'Apply enrichment.');

      parent = ui.divH([runLink, ui.divH([editLink, deleteLink], 'power-pack-enrichment-actions')], 'power-pack-enrichment-row');
      return parent;
    }));
    return root;
  });
}

function showEnrichDialog(mainTable: DG.TableInfo, df: DG.DataFrame, dbColName: string, onSave: Function, enrichName?: string, tableQuery?: DG.TableQuery) {
  if (!tableQuery) {
    tableQuery = DG.TableQuery.fromTable(mainTable).build();
    tableQuery.fields = [dbColName];
  }
  const pivotView: DG.VisualDbQueryEditor = DG.VisualDbQueryEditor.fromQuery(tableQuery);
  pivotView.showAddToWorkspaceBtn = false;

  const nameInput: DG.InputBase = ui.input.string('Name', {nullable: false, tooltipText: 'Provide name for the enrichment.', value: enrichName});
  const root = ui.wait(async () => {
    await pivotView.isInit();
    // await pivotView.setSingleColumnMode(dbColName);
    ui.setDisplay(pivotView.groupByTag.root.parentElement!, false);
    ui.setDisplay(pivotView.orderTag.root.parentElement!, false);
    ui.setDisplay(pivotView.whereTag.root.parentElement!, false);
    ui.setDisplay(pivotView.aggregateTag.root.parentElement!, false);
    ui.setDisplay(pivotView.pivotTag.root.parentElement!, false);
    ui.setDisplay(pivotView.havingTag.root.parentElement!, false);
    ui.setDisplay(pivotView.root.querySelector('.grok-pivot-grid')!, false);

    pivotView.root.style.flexGrow = 'initial';
    pivotView.root.style.minHeight = '50px';
    pivotView.mainTag.root.children[pivotView.mainTag.root.children.length - 1].classList.add('d4-tag-disabled');
    ui.tooltip.bind(pivotView.mainTag.root, `Enrich supports joining using only a single column. Selected column: “${dbColName}”.`);
    pivotView.grid.root.style.width = '100%';
    return ui.divV([
      pivotView.root,
      pivotView.grid.root
    ]);
  });

  const isValidInput = ui.input.string('');
  isValidInput.nullable = false;
  isValidInput.root.hidden = true;

  const dialog = ui.dialog({title: `Enrich ${dbColName}`})
    .add(nameInput)
    .add(isValidInput)
    .add(root)
    .onOK(async () => {
      pivotView.refreshQuery();
      const enrichment = convertTableQueryToEnrichment(pivotView.query, nameInput.value,
        mainTable.tags.get(DG.Tags.TableSchema), mainTable.friendlyName, dbColName);
      await saveEnrichment(enrichment, pivotView.query.connection.nqName);
      if (enrichName && enrichment.name != enrichName) {
        await deleteEnrichment(pivotView.query.connection.nqName, enrichment.keySchema,
          enrichment.keyTable, enrichment.keyColumn, enrichName);
      }
      onSave();
    });


  dialog.addButton('ENRICH', async () => {
    const progress = DG.TaskBarProgressIndicator.create('Enriching...');
    try {
      pivotView.refreshQuery();
      await executeEnrichQuery(pivotView.query, df, dbColName);
    } finally {
      progress.close();
    }
  });

  dialog.root.style.height = '600px';
  dialog.root.style.width = '800px';
  const saveButton = dialog.getButton('OK');
  saveButton.querySelector('span')!.textContent = 'SAVE';
  const enrichButton = dialog.getButton('ENRICH');
  const isNotValidPivot = () => pivotView.query.joins.length === 0;
  pivotView.onChanged.subscribe((_) => {
    const isNotValid = isNotValidPivot();
    ui.setClass(enrichButton, 'disabled', isNotValid);
    isValidInput.value = isNotValid ? '' : 'not empty';
  });

  const isPivotNotValid = isNotValidPivot();
  ui.setClass(enrichButton, 'disabled', isPivotNotValid);
  isValidInput.value = isPivotNotValid ? '' : 'non empty';
  dialog.show();
}

async function getEnrichConfigsNames(connectionNqName: string, schema: string, table: string, column: string): Promise<string[]> {
  try {
    const files = await grok.dapi.files
      .list(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}/${schema}/${table}/${column}`);
    return files.map((f) => f.name.substring(0, f.name.lastIndexOf('.')));
  } catch (_: any) {
    return [];
  }
}

async function deleteEnrichment(connectionNqName: string, schema: string, table: string, column: string, name: string): Promise<boolean> {
  try {
    await grok.dapi.files
      .delete(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}/${schema}/${table}/${column}/${name}.json`);
    return true;
  } catch (_: any) {
    new DG.Balloon().error('Something went wrong while deleting enrichment. Please, try again.');
    return false;
  }
}

async function readEnrichConfig(connectionNqName: string, schema: string, table: string, column: string, name: string): Promise<Enrichment | null> {
  try {
    const data = await grok.dapi.files
      .readAsText(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}/${schema}/${table}/${column}/${name}.json`);
    return JSON.parse(data);
  } catch (_: any) {
    return null;
  }
}

async function saveEnrichment(e: Enrichment, connectionNqName: string): Promise<void> {
  try {
    await grok.dapi.files
      .writeAsText(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}/${e.keySchema}/${e.keyTable}/${e.keyColumn}/${e.name}.json`,
        JSON.stringify(e));
  } catch (_: any) {
    new DG.Balloon().error('Could not save enrichment config. Please try again.');
  }
}

function nqNameToPath(connectionNqName: string) {
  return connectionNqName.replace(':', '_');
}

function convertTableQueryToEnrichment(query: DG.TableQuery, name: string, keySchema: string, keyTable: string, keyColumn: string): Enrichment {
  return {name: name, keySchema: keySchema, keyTable: keyTable, keyColumn: keyColumn, fields: query.fields, joins: query.joins};
}

function convertEnrichmentToQuery(e: Enrichment, conn: DG.DataConnection): DG.TableQuery {
  const query = conn.tableQuery();
  query.table = `${e.keySchema}.${e.keyTable}`;
  query.fields = e.fields;
  query.joins = e.joins;
  return query;
}

export async function runEnrichmentFromConfig(conn: DG.DataConnection, schema: string, table: string, column: string, name: string, df: DG.DataFrame): Promise<void> {
  const config = await readEnrichConfig(conn.nqName, schema, table, column, name);
  if (!config) {
    new DG.Balloon().error('Could not find enrichment. Please try again or create a new one.');
    return;
  }
  const query = convertEnrichmentToQuery(config, conn);
  await executeEnrichQuery(query, df, column);
}

async function executeEnrichQuery(query: DG.TableQuery, df: DG.DataFrame, keyCol: string): Promise<void> {
  try {
    query.limit = undefined;
    const keyColValues = df.getCol(keyCol).toList();
    let res: DG.DataFrame | undefined;
    for (let i = 0; i < keyColValues.length; i += 200) {
      const inValues = keyColValues.slice(i, i + 200);
      query.where = [{field: `${query.table}.${keyCol}`, pattern: `in (${inValues.join(',')})`, dataType: df.getCol(keyCol).type}];
      const queryCall = query.prepare();
      const run = await queryCall.call(false, undefined,
        {processed: true, report: false});
      const queryResult = run.getOutputParamValue();
      if (!queryResult)
        break;
      if (!res)
        res = queryResult;
      else
        res.append(queryResult, true);
    }
    if (!res) {
      grok.shell.error('Something went wrong and query result returned null.');
      return;
    }

    await grok.data.detectSemanticTypes(res);
    const joinTable = DG.Func.byName('JoinTables');
    const previousName = df.name;
    joinTable.applySync({
      'table1': df,
      'table2': res,
      'keys1': [keyCol],
      'keys2': [keyCol],
      'values1': df.columns.names(),
      'values2': res.columns.names().filter((n) => n !== keyCol),
      'joinType': 'left',
      'inPlace': true
    });
    df.name = previousName;
  } catch (e: any) {
    grok.shell.error(`Failed to enrich:\n ${e}`);
  }
}

interface Enrichment {
  name: string;
  keySchema: string;
  keyTable: string,
  keyColumn: string;
  fields: string[];
  joins: DG.TableJoin[];
}
