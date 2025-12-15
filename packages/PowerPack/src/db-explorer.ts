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
        editor = new DBExplorerEditor(() => editor ? onSaveAction(editor) : Promise.resolve(), false);
        const editorUI = await editor.getUI();
        editor.setConfig(config);
        const view = DG.View.create();
        view.name = `${connection.name}.${config.schemaName} Identifiers`;
        view.root.appendChild(editorUI);
        grok.shell.addPreview(view);
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
      !args?.args?.menu)
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
        tooltipText: 'Select schema for identifiers', nullable: false} );

      ui.dialog('Select Schema for Identifiers Configuration')
        .add(ui.form([schemaChoice]))
        .onOK(async () => {
          let editor: DBExplorerEditor | null = null;
          editor = new DBExplorerEditor(() => {
            editor ? onSaveAction(editor) : Promise.resolve();
          }, false); // no connection section needs to be rendered here

          const uiEl = await editor.getUI();
          // find if there is an existing config for this connection-schema
          const existingConfig = savedConfigs?.find((c) =>
            c.nqName === connection.nqName && c.dataSourceName === connection.dataSource &&
            c.schemaName === schemaChoice.value!);
          if (existingConfig) {
            editor.setConfig(existingConfig);
            grok.shell.info('Identifiers configuration for this connection and schema already exists.\n Importing...');
          } else {
            editor.setConfig({
              connectionName: connection.name,
              nqName: connection.nqName,
              dataSourceName: connection.dataSource,
              schemaName: schemaChoice.value ?? undefined,
            }, false);
          } // no warning because schema is deliberately left empty for user to fill in

          const view = DG.View.create();
          view.name = `${connection.name}.${schemaChoice.value!} Identifiers`;
          view.root.appendChild(uiEl);

          grok.shell.addView(view);
        }).show();
    }, 4, {description: 'Configure identifier settings for this database connection'});
  });


  async function getSavedConfigs(): Promise<DBExplorerConfig[]> {
    let configs: DBExplorerConfig[] = [];
    try {
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
  };

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
  const conIdToNqName = new Map<string, {name?: string, nqName?: string} | null>();

  function processCell(cell: DG.Cell) {
    ui.empty(explorePanelRoot);
    explorePanelRoot.appendChild(ui.divText('No data to explore.'));
    ui.setUpdateIndicator(explorePanelRoot, true);
    if (processTimer != null)
      clearTimeout(processTimer);
    processTimer = setTimeout(async () => {
      if (!cell?.dart || (cell.rowIndex ?? -1) < 0 || !cell.column?.dart ||
        !cell.dataFrame?.dart || grok.shell.o != cell.column
      ) { // no need to process anything
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
      if (!conIdToNqName.has(cId)) {
        const connection = await grok.dapi.connections.find(cId);
        const perm = connection ? await grok.dapi.permissions.check(connection, 'View') : false;
        conIdToNqName.set(cId, connection && perm ? {name: connection.name, nqName: connection.nqName} : null);
      }
      const nqNameInfo = conIdToNqName.get(cId);
      if (!nqNameInfo || !nqNameInfo?.nqName || !nqNameInfo?.name) {
        ui.setUpdateIndicator(explorePanelRoot, false);
        return;
      }
      const handler = DG.ObjectHandler.list().find((h) => { // for entity method does not work, because of instanceof checks.
        // need to do it manually// underneeth, its doing the same with listing all handlers
        return h.type === DB_EXPLORER_OBJ_HANDLER_TYPE && (h as DBExplorerObjectHandler).connectionNqName === nqNameInfo.nqName! && (h as DBExplorerObjectHandler).schemaName === schemaName;
      });
      if (handler) {
        const acc = (handler as DBExplorerObjectHandler).renderPropertiesFromDfRow(cell.row, tableName);
        ui.empty(explorePanelRoot);
        explorePanelRoot.appendChild(acc);
      } else {
        // if there is no handler, at this point we can just register it
        console.log('Registering generic handler for ', nqNameInfo.nqName, schemaName);
        const exp = DBExplorer.initFromConfig({
          connectionName: nqNameInfo.name!,
          nqName: nqNameInfo.nqName!,
          schemaName: schemaName,
          joinOptions: [],
          entryPoints: {},
        });
        const newHandler = exp.genericValueHandler; // nex time it will be found
        const acc = newHandler.renderPropertiesFromDfRow(cell.row, tableName);
        ui.empty(explorePanelRoot);
        explorePanelRoot.appendChild(acc);
      }
      ui.setUpdateIndicator(explorePanelRoot, false);
    }, 300);
  }


  grok.events.onAccordionConstructed.subscribe((acc) => {
    if (!(acc.context instanceof DG.Column))
      return;
    const col = acc.context as DG.Column;
    const df = col.dart ? col.dataFrame : null;
    if (!df || !df.tags.get(DG.Tags.DataConnectionId) || !col.tags.get(DG.Tags.DbSchema) ||
            !col.tags.get(DG.Tags.DbTable) || !col.tags.get(DG.Tags.DbColumn))
      return;
    ui.empty(explorePanelRoot);
    explorePanelRoot.appendChild(ui.divText('Select a cell to explore its value...'));
    acc.addPane('Explore', () => explorePanelRoot);
    // if there is a current cell and its in the given column, process it
    console.log('Accordion constructed event received in db-explorer setup.');
    if (col.dataFrame.currentCell && col.dataFrame.currentCell.column === col)
      processCell(col.dataFrame.currentCell!);
  });

  grok.events.onCurrentCellChanged.subscribe((cell?: DG.Cell) => {
    console.log('Current cell changed event received in db-explorer setup.');
    if (!cell?.dart || !cell.column?.dart || !cell.dataFrame?.dart || grok.shell.o != cell.column)
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
