/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';
import {DBExplorerConfig} from '@datagrok-libraries/db-explorer/src/types';
import {DBExplorer} from '@datagrok-libraries/db-explorer/src/db-explorer';
import {DBExplorerEditor} from '@datagrok-libraries/db-explorer/src/editor';

export const WIDGETS_STORAGE = 'widgets';

export interface UserWidgetSettings {
  factoryName?: string;
  caption?: string;
  ignored?: boolean;
}

export interface UserWidgetsSettings {
  [index: string]: UserWidgetSettings;
}

export let settings: UserWidgetsSettings;

export function getSettings(): UserWidgetsSettings {
  if (!settings) {
    const savedSettings: {[key: string]: any} = grok.userSettings.get(WIDGETS_STORAGE) ?? {};
    for (const key of Object.keys(savedSettings))
      savedSettings[key] = JSON.parse(savedSettings[key]);
    settings = savedSettings;
  }
  return settings;
}

export function saveSettings(): void {
  const s: {[key: string]: any} = {};
  for (const key of Object.keys(settings))
    s[key] = JSON.stringify(settings[key]);
  grok.userSettings.addAll(WIDGETS_STORAGE, s);
}


function initWidgetHost(host: HTMLDivElement, w: DG.Widget, title?: string) {
  function remove(): void {
    host.remove();
    if (w.factory?.name) {
      const widgetSettings = settings[w.factory.name] ?? (settings[w.factory.name] = { });
      widgetSettings.ignored = true;
      saveSettings();
    }
  }

  if (w.props.hasProperty('order'))
    host.style.order = w.props.order;

  const header = host.querySelector('.d4-dialog-header')! as HTMLElement;
  if (!title) {
    header.classList.add('d4-dialog-header-hidden');
    w.root.appendChild(ui.icons.close(remove, 'Remove'));
  } else
    header.appendChild(ui.icons.close(remove, 'Remove'));

  if (w.root.classList.contains('widget-narrow'))
    host.classList.add('widget-narrow');
  if (w.root.classList.contains('widget-wide'))
    host.classList.add('widget-wide');

  host.querySelector('.power-pack-widget-content')!.appendChild(w.root);
  ui.tools.setHoverVisibility(host, Array.from(host.querySelectorAll('i')));
  if (w.factory?.name) {
    const widgetSettings = settings[w.factory.name] ?? (settings[w.factory.name] = { });
    if (widgetSettings.ignored === undefined || widgetSettings.ignored === null)
      widgetSettings.ignored = false;
    saveSettings();
  }
}

function createWidgetHost(title: string): HTMLDivElement {
  const header = ui.div([ui.divText(title, 'd4-dialog-title')], 'd4-dialog-header');
  const host = ui.box(null, 'power-pack-widget-host');
  host.appendChild(header);
  host.appendChild(ui.box(null, 'power-pack-widget-content'));
  return host;
}

export function widgetHostFromFunc(f: DG.Func) {
  const title = f.options['showName'] === 'false' ? '' : f.friendlyName;
  const host: HTMLDivElement = createWidgetHost(title);
  const contentDiv: HTMLElement = (host.querySelector('.power-pack-widget-content')!) as HTMLElement;

  f.apply().then(function(w: DG.Widget) {
    if (w) {
      w.factory = f;
      initWidgetHost(host, w, title);
    } else
      host.remove();
  })
    .catch((e) => {
      host.style.display = 'none';
      host.remove();
      console.error(`Error creating widget ${f.name}`, e);
    })
    .finally(() => ui.setUpdateIndicator(contentDiv, false, ''));

  setTimeout(() => {
    if (contentDiv!.children.length == 0)
      ui.setUpdateIndicator(contentDiv, true, '');
  }, 1000);

  if (f.options['order'] !== null)
    host.style.order = f.options['order'];

  return host;
}


export function widgetHost(w: DG.Widget/*, widgetHeader?: HTMLDivElement*/): HTMLElement {
  const title = w.props.hasProperty('caption') ? w.props.caption ?? '' : '';
  const host = createWidgetHost(title);
  initWidgetHost(host, w, title);
  //widgetHeader ??= ui.div();

  return host;
}

/** Sets up editing of the db explorers config(s) */
export async function setupGlobalDBExplorer() {
  // first thing we do is to listen to the browse tree node expanding, to add stuff to connections
  let savedConfigs: DBExplorerConfig[] | null = null;
  const deferredTreeNodes: DG.TreeViewGroup[] = [];
  const processedNodes = new Set<DG.TreeViewGroup>();

  const processConnectionNode = (node: DG.TreeViewGroup) => {
    if (processedNodes.has(node))
      return;
    processedNodes.add(node);
    const connection = DG.toJs(node.value) as DG.DataConnection;
    const matchingConfigs = savedConfigs?.filter((c) =>
      c.nqName === connection.nqName && c.dataSourceName === connection.dataSource);
    matchingConfigs?.forEach((config) => {
      const itemElement = ui.divH([ui.iconFA('fingerprint', ()=>{}), ui.divText(`Identifiers - ${config.schemaName}`)],
        {style: {gap: '6px', alignItems: 'center'}});
      const tnItem = node.item(itemElement, config);
      tnItem.onSelected.subscribe(async () => {
        let editor: DBExplorerEditor | null = null;
        editor = new DBExplorerEditor(() => editor ? onSaveAction(editor) : Promise.resolve());
        const editorUI = await editor.getUI();
        editor.setConfig(config);
        const view = DG.View.create();
        view.name = `${connection.name} - ${config.schemaName} DB Explorer Editor`;
        view.root.appendChild(editorUI);
        grok.shell.addPreview(view);
      });
      tnItem.root.style.order = '-1'; // make sure it is on top
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
          processConnectionNode(connectionNode);
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

  // set up the entry points for connections.

  grok.events.onContextMenu.subscribe((args) => {
    if (!(DG.toJs(args?.args?.item) instanceof DG.TreeViewGroup) ||
    !(DG.toJs(args.args.item.value) instanceof DG.DataConnection) ||
      !args?.args?.menu)
      return;
    const connection = DG.toJs(args.args.item.value) as DG.DataConnection;
    const menu = args.args.menu;
    menu.item('Configure Identifiers', async () => {
      let editor: DBExplorerEditor | null = null;
      editor = new DBExplorerEditor(() => {
        editor ? onSaveAction(editor) : Promise.resolve();
      });

      const uiEl = await editor.getUI();
      editor.setConfig({
        connectionName: connection.name,
        nqName: connection.nqName,
        dataSourceName: connection.dataSource,
      }, false); // no warning because schema is deliberately left empty for user to fill in

      const view = DG.View.create();
      view.name = connection.name + ' DB Explorer Editor';
      view.root.appendChild(uiEl);

      grok.shell.addView(view);
    });
  });
}
