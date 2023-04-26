import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import '../../css/demo.css';

type Direction = {
  category: string,
  name: string,
  description: string,
  func: DG.Func
};


export class DemoView extends DG.ViewBase {
  dockPanel: DG.DockNode = new DG.DockNode(undefined);
  tree: DG.TreeViewGroup = ui.tree();
  searchInput: DG.InputBase = ui.searchInput('', '');

  constructor() {
    super();
    this._initDockPanel();
    this._initContent();
  }

  static findDemoFunc(demoPath: string) {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  async startDemoFunc(func: DG.Func, viewPath: string) {
    this._closeAll();

    ui.setUpdateIndicator(grok.shell.tv.root, true);
    grok.shell.windows.showHelp = true;

    await func.apply();
    ui.setUpdateIndicator(grok.shell.tv.root, false);

    grok.shell.v.path.includes('/apps/Tutorials/Demo') ?
      grok.shell.v.path = grok.shell.v.basePath = `/${viewPath}` :
      grok.shell.v.path = grok.shell.v.basePath = `/apps/Tutorials/Demo/${viewPath}`;
  }


  private _closeAll() {
    grok.shell.closeAll();
    this._closeDemoScript();
  }

  private _closeDemoScript() {
    const scriptDockNode = Array.from(grok.shell.dockManager.rootNode.children)[1];
    if (scriptDockNode?.container.containerElement.classList.contains('tutorials-demo-script-container')) {
      grok.shell.dockManager.close(scriptDockNode);
    }
  }

  private _closeDockPanel() {
    const panelDockNode = Array.from(grok.shell.dockManager.rootNode.children)[0];
    if (panelDockNode?.container.containerElement.classList.contains('tutorials-demo-container')) {
      grok.shell.dockManager.close(panelDockNode);
    }
  }

  private _initContent() {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;

    this.root.innerHTML = '';
    
    this.name = 'Demo app';

    const tree = ui.tree();
    tree.root.classList.add('demo-app-group-view');
    
    const tempGroups:String[] = [];

    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
      const path = pathOption.split('|').map((s) => s.trim());
      tempGroups.push(path[0]);
    }

    const groups: String[] = [...new Set(tempGroups)];

    for (let i=0; i<groups.length; i++){
      const name = groups[i] as string;
      tree.group(name, null, true).root.lastChild?.appendChild(this.groupRoot(name));
    }
    this.root.append(ui.div([tree.root], 'grok-gallery-grid'));
  }

  groupRoot (groupName: string) {
    const root = ui.div([], 'demo-app-group-view grok-gallery-grid');

    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      if (f.options[DG.FUNC_OPTIONS.DEMO_PATH].includes(groupName)) {
        const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
        const path = pathOption.split('|').map((s) => s.trim());
        const demo = path[path.length - 1];

        const imgPath = `${_package.webRoot}images/demoapp/${f.name}.jpg`;
        const img = ui.div('', 'ui-image');

        fetch(imgPath)
          .then(res => {
            if (res.ok)
              img.style.backgroundImage = `url(${imgPath})`
            else
              img.style.backgroundImage = `url(${_package.webRoot}images/demoapp/emptyImg.jpg)`
          })
          .catch()

        let item = ui.card(ui.divV([
          img,
          ui.div([demo], 'tutorials-card-title'),
          ui.div([f.description], 'tutorials-card-description')
        ], 'demo-app-card'));

        item.onclick = () => {
          let node = this.tree.items.find(node => node.text == demo)?.root;
          node?.click();
        };

        if (f.description != '')
          ui.tooltip.bind(item, f.description)

        root.append(item);
      }
    }
    return root
  }
  
  nodeView(viewName: string) {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;

    grok.shell.closeAll();

    const view = grok.shell.newView(viewName);
    view.basePath = '/apps/Tutorials/Demo';
    view.path = `/${viewName}`;

    const root = ui.div([], 'grok-gallery-grid');
    const tree = ui.tree();
    let treeNode = tree.group(viewName, null, true);

    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      if (f.options[DG.FUNC_OPTIONS.DEMO_PATH].includes(viewName)) {
        const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
        const path = pathOption.split('|').map((s) => s.trim());
        const demo = path[path.length - 1];

        const imgPath = `${_package.webRoot}images/demoapp/${f.name}.jpg`;
        const img = ui.div('', 'ui-image');

        fetch(imgPath)
          .then(res => {
            if (res.ok)
              img.style.backgroundImage = `url(${imgPath})`
            else
              img.style.backgroundImage = `url(${_package.webRoot}images/demoapp/emptyImg.jpg)`
          })
          .catch()

        let item = ui.card(ui.divV([
          img,
          ui.div([demo], 'tutorials-card-title'),
          ui.div([f.description], 'tutorials-card-description')
        ], 'demo-app-card'));

        item.onclick = () => {
          let node = this.tree.items.find(node => node.text == demo)?.root;
          node?.click();
        };

        if (f.description)
          ui.tooltip.bind(item, f.description)

        root.append(item);
      }
    }
    treeNode.root.lastChild?.appendChild(root);
    tree.root.classList.add('demo-app-group-view');
    grok.shell.v.root.append(ui.div([tree.root], 'grok-gallery-grid'));
  }

  private _initDockPanel() {
    if (this._isDockPanelInit()) {
      this._closeDockPanel();
    }

    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
      const path = pathOption.split('|').map((s) => s.trim());
      const folder = this.tree.getOrCreateGroup(path.slice(0, path.length - 1).join(' | '));
      const item = folder.item(path[path.length - 1]);

      if (folder.text === 'Viewers')
        folder.root.style.order = '10'
      else  
        folder.root.style.order = '1';

      item.root.onmouseover = (event) => {
        const packageMessage = `Part of the ${f.package.name} package`;
        ui.tooltip.show(f.description ? ui.divV([f.description, ui.element('br'), packageMessage]) : ui.div(packageMessage),
          event.clientX, event.clientY);
      };

      item.root.onmouseout = (_) => {
        ui.tooltip.hide();
      };
    }

    this.searchInput.onChanged(() => {
      const dom = this.tree.root.getElementsByClassName('d4-tree-view-node');

      for (let i = 0; i < dom.length; i++) {
        const item = dom[i] as HTMLElement;
        if (item.innerText.toLowerCase().includes(this.searchInput.value.toLowerCase()))
          item.classList.remove('hidden');
        else
          item.classList.add('hidden');
      }
    });

    this.searchInput.input.onkeyup = (event) => {
      if (event.key === 'Escape')
        this.searchInput.fireChanged();
    };

    const closeIcon = this.searchInput.root.getElementsByClassName('ui-input-icon-right')[0] as HTMLElement;
    closeIcon.onclick = () => {
      this.searchInput.value = '';
      this.searchInput.fireChanged();
    };

    let homeNode = this.tree.group('Home');
    homeNode.root.classList.add('demo-app-tree-home-node');
    homeNode.root.getElementsByClassName('d4-tree-view-node')[0]?.prepend(ui.iconFA('home'));
    homeNode.root.getElementsByClassName('d4-tree-view-tri')[0].remove();

    DG.debounce(this.tree.onSelectedNodeChanged, 300).subscribe(async (value) => {
      if (DemoScript.currentObject) {
        DemoScript.currentObject.cancelScript();
        const scriptDockNode = Array.from(grok.shell.dockManager.rootNode.children)[1];
        if (scriptDockNode.container.containerElement.classList.contains('tutorials-demo-script-container')) {
            grok.shell.dockManager.close(scriptDockNode);
            grok.shell.closeAll();
            grok.shell.addView(new DemoView());
        }
      }

      if (value.root.classList.contains('d4-tree-view-item')) {
        const categoryName = value.root.parentElement?.parentElement
          ?.getElementsByClassName('d4-tree-view-group-label')[0].innerHTML;
        const viewerName = value.text;
        const demoFunc = DemoView.findDemoFunc(`${categoryName} | ${viewerName}`);
        const demoPath = `${categoryName}/${viewerName}`;
        await this.startDemoFunc(demoFunc, demoPath);
        this.tree.root.focus();
      } else if (value.root.classList.contains('demo-app-tree-home-node')) { 
        this._initContent();
        grok.shell.windows.showToolbox = false;
        grok.shell.windows.showHelp = false;
        grok.shell.windows.showProperties = false;
        grok.shell.closeAll();
        const view = grok.shell.addView(this);
        view.basePath = '/apps/Tutorials/Demo';
        view.path = `/`;
      } else {
        this.tree.root.focus();
        this.nodeView(value.text);
      }
    });

    this.dockPanel = grok.shell.dockManager.dock(ui.panel([
      this.searchInput.root,
      this.tree.root,
    ]), 'left', null, 'Categories');
    this.dockPanel.container.containerElement.classList.add('tutorials-demo-container');

    this.tree.root.classList.add('demo-app-tree-group');

    this._initWindowOptions();

    // TODO: on click on viewer demo set viewer help url in property panel (func helpUrl)
    // TODO: implement search in demo - search on meta.keywords, name, description

    // TODO: if there empty space - add viewer/filter/etc.

    // TODO: add to script demo class grok.shell.windows.showPropertyPanel = true and showHelp = false
    // TODO: add GIS
    // TODO: add breadcrumbs instead of name
  }

  private _isDockPanelInit(): boolean {
    const panelDockNode = Array.from(grok.shell.dockManager.rootNode.children)[0];
    return panelDockNode?.container.containerElement.classList.contains('tutorials-demo-container');
  }

  private _initWindowOptions() {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showRibbon = true;
    grok.shell.windows.showHelp = true;
    grok.shell.windows.showProperties = false;

    grok.shell.windows.help.syncCurrentObject = false;
  }
}
