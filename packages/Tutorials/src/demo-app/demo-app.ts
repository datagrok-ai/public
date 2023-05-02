import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {sortFunctionsByDemoPath} from './utils';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import '../../css/demo.css';

const DIRECTIONS: string[] = [
  'Cheminformatics',
  'Bioinformatics',
  'Viewers',
];

// const CHEMINFORMATICS_CATEGORIES: string[] = [];
// const BIOINFORMATICS_CATEGORIES: string[] = [];
// const VIEWERS_CATEGORIES: string[] = [
//   'Data flow and hierarchy',
//   'Data separation',
//   'General',
//   'Geographical',
//   'Input and edit',
//   'Statistical',
//   'Time and date'
// ];

// const SUBCATEGORIES = {
//   Cheminformatics: CHEMINFORMATICS_CATEGORIES,
//   Bioinformatics: BIOINFORMATICS_CATEGORIES,
//   Viewers: VIEWERS_CATEGORIES,
// };

type DemoFunc = {
  name: string;
  func: DG.Func,
  category: string;
  path: string;
  keywords: string;
  imagePath: string;
};


export class DemoView extends DG.ViewBase {
  dockPanel: DG.DockNode = new DG.DockNode(undefined);
  tree: DG.TreeViewGroup = ui.tree();
  searchInput: DG.InputBase = ui.searchInput('', '');
  funcs: DemoFunc[] = [];

  constructor() {
    super();
    this._initFunctions();
    this._initDockPanel();
    this._initContent();
  }


  static findDemoFunc(demoPath: string): DG.Func {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  async startDemoFunc(func: DG.Func, viewPath: string): Promise<void> {
    const path = viewPath.split('|').map((s) => s.trim()).join('/');

    this._closeAll();

    ui.setUpdateIndicator(grok.shell.tv.root, true);
    grok.shell.windows.showHelp = true;

    await func.apply();
    ui.setUpdateIndicator(grok.shell.tv.root, false);

    grok.shell.v.path.includes('/apps/Tutorials/Demo') ?
      grok.shell.v.path = grok.shell.v.basePath = `/${path}` :
      grok.shell.v.path = grok.shell.v.basePath = `/apps/Tutorials/Demo/${path}`;
  }


  private _closeAll(): void {
    grok.shell.closeAll();
    this._closeDemoScript();
  }

  private _closeDemoScript(): void {
    const scriptDockNode = Array.from(grok.shell.dockManager.rootNode.children)[1];
    if (scriptDockNode?.container.containerElement.classList.contains('tutorials-demo-script-container'))
      grok.shell.dockManager.close(scriptDockNode);
  }

  private _closeDockPanel(): void {
    const panelDockNode = Array.from(grok.shell.dockManager.rootNode.children)[0];
    if (panelDockNode?.container.containerElement.classList.contains('tutorials-demo-container'))
      grok.shell.dockManager.close(panelDockNode);
  }

  private _initFunctions(): void {
    const funcs = DG.Func.find({meta: {'demoPath': null}}).sort(sortFunctionsByDemoPath);
    for (let i = 0; i < DIRECTIONS.length; ++i) {
      const directionFuncs = funcs.filter((func) => {
        return (func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DIRECTIONS[i]);
      });

      for (let j = 0; j < directionFuncs.length; ++j) {
        let imgPath = `${_package.webRoot}images/demoapp/${directionFuncs[j].name}.jpg`;
        fetch(imgPath)
          .then(res => {
            imgPath = res.ok ? imgPath : `${_package.webRoot}images/demoapp/emptyImg.jpg`;
          })
          .catch();

        const path = directionFuncs[j].options[DG.FUNC_OPTIONS.DEMO_PATH] as string;
        const pathArray = path.split('|').map((s) => s.trim());

        this.funcs[this.funcs.length] = {
          name: pathArray[pathArray.length - 1],
          func: directionFuncs[j],
          category: DIRECTIONS[i],
          path: path,
          keywords: (directionFuncs[j].options['keywords'] !== undefined) ? directionFuncs[j].options['keywords'] as string : '',
          imagePath: imgPath,
        };
      }
    }
  }

  private _initContent(): void {
    this._initWindowOptions();

    this.root.innerHTML = '';
    this.name = 'Demo app';

    const tree = ui.tree();
    tree.root.classList.add('demo-app-group-view');
    
    for (let i = 0; i < DIRECTIONS.length; i++) {
      const name = DIRECTIONS[i];
      tree.group(name, null, true).root.lastChild?.appendChild(this._groupRoot(name));
    }

    this.root.append(ui.div([tree.root], 'grok-gallery-grid'));
  }

  private _createViewRootElement(viewOrGroupName: string): HTMLDivElement {
    const root = ui.div([]);

    const directionFuncs = this.funcs.filter((func) => {
      return (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(viewOrGroupName);
    });

    for (let i = 0; i < directionFuncs.length; i++) {
      const img = ui.div('', 'ui-image');
      img.style.backgroundImage = `url(${directionFuncs[i].imagePath}`;

      let item = ui.card(ui.divV([
        img,
        ui.div([directionFuncs[i].name], 'tutorials-card-title'),
        ui.div([directionFuncs[i].func.description], 'tutorials-card-description')
      ], 'demo-app-card'));

      item.onclick = () => {
        const node = this.tree.items.find(node => node.text === directionFuncs[i].name)?.root;
        node?.click();
      };

      const packageMessage = `Part of the ${directionFuncs[i].func.package.name === 'Tutorials' ?
        'platform core' : `${directionFuncs[i].func.package.name} package`}`;
      ui.tooltip.bind(item, () => directionFuncs[i].func.description ?
        ui.divV([directionFuncs[i].func.description, ui.element('br'), packageMessage]) : ui.div(packageMessage));

      root.append(item);
    }

    return root;
  }

  private _groupRoot(groupName: string): HTMLDivElement {
    const root = this._createViewRootElement(groupName);
    root.classList.add('demo-app-group-view');
    root.classList.add('grok-gallery-grid');

    return root;
  }
  
  // TODO: in view show tree too
  // TODO: make config for easy showing (Cheminformatics -> Bio -> Viewers -> ...)
  // TODO: fix search, also on input show them in view, make it work on description and keywords

  // TODO: in DemoScript - loading circle to end when delay ends (make interactive icon) to disappear like in balloon

  // TODO: try to load breadcrumbs instead of view name (onViewLoaded)

  // TODO: pause on exceptions in browser and vscode, check PowerGrid problem

  // TODO: demos: FileManager: show files in folder with demo, show molecules table
  // TODO: demos: Table linking: make custom view with 2 grids and link them with the proper API (filter in one table will set uo the second table)
  // TODO: demos: Grid customizations (in PowerGrid): have to add some sparklines, also add frozen columns (check in PowerGrid)

  nodeView(viewName: string, path: string): void {
    this._initWindowOptions();
    this._closeAll();

    const view = grok.shell.newView(viewName);
    view.basePath = '/apps/Tutorials/Demo';
    view.path = `/${path}`;

    const root = this._createViewRootElement(viewName);
    root.classList.add('grok-gallery-grid');

    const tree = ui.tree();
    const treeGroup = tree.group(viewName, null, true);

    treeGroup.root.lastChild?.appendChild(root);
    tree.root.classList.add('demo-app-group-view');
    view.root.append(ui.div([tree.root], 'grok-gallery-grid'));
  }

  private _createHomeNode(): void {
    const homeNode = this.tree.group('Home');
    homeNode.root.classList.add('demo-app-tree-home-node');
    (homeNode.root.firstElementChild as HTMLElement).dataset.name = 'Home';
    homeNode.root.getElementsByClassName('d4-tree-view-node')[0]?.prepend(ui.iconFA('home'));
    homeNode.root.getElementsByClassName('d4-tree-view-tri')[0].remove();
  }

  private _initDockPanel(): void {
    if (this._isDockPanelInit())
      this._closeDockPanel();

    this._createHomeNode();

    for (let i = 0; i < DIRECTIONS.length; ++i) {
      const directionFuncs = this.funcs.filter((func) => {
        return (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DIRECTIONS[i]);
      });

      for (let j = 0; j < directionFuncs.length; ++j) {
        const path = directionFuncs[j].path.split('|').map((s) => s.trim());

        if (path.length > 2) {
          let groupPath = path[0];
          let treePath = this.tree.getOrCreateGroup(path[0], {path: groupPath});
          (treePath.root.firstElementChild as HTMLElement).dataset.name = path[0];
          for (let i = 1; i < path.length - 1; i++) {
            groupPath += `/${path[i]}`;
            treePath = treePath.getOrCreateGroup(path[i], {path: groupPath});
            (treePath.root.firstElementChild as HTMLElement).dataset.name = path[i];
          }

          const item = treePath.item(directionFuncs[j].name, {path: directionFuncs[j].path});
          item.root.onmouseover = (event) => {
            const packageMessage = `Part of the ${directionFuncs[j].func.package.name === 'Tutorials' ?
              'platform core' : `${directionFuncs[j].func.package.name} package`}`;
            ui.tooltip.show(directionFuncs[j].func.description ?
              ui.divV([directionFuncs[j].func.description, ui.element('br'), packageMessage]) :
              ui.div(packageMessage), event.clientX, event.clientY);
          };

          item.root.onmouseout = (_) => {
            ui.tooltip.hide();
          };
        } else {
          const folder = this.tree.getOrCreateGroup(directionFuncs[j].category, {path: path[0]});
          (folder.root.firstElementChild as HTMLElement).dataset.name = directionFuncs[j].category;
          const item = folder.item(directionFuncs[j].name, {path: directionFuncs[j].path});

          item.root.onmouseover = (event) => {
            const packageMessage = `Part of the ${directionFuncs[j].func.package.name === 'Tutorials' ?
              'platform core' : `${directionFuncs[j].func.package.name} package`}`;
            ui.tooltip.show(directionFuncs[j].func.description ?
              ui.divV([directionFuncs[j].func.description, ui.element('br'), packageMessage]) :
              ui.div(packageMessage), event.clientX, event.clientY);
          };

          item.root.onmouseout = (_) => {
            ui.tooltip.hide();
          };
        }
      }
    }

    this.searchInput.onChanged(() => {
      this._searchItem();
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

    DG.debounce(this.tree.onSelectedNodeChanged, 300).subscribe(async (value) => {
      if (DemoScript.currentObject) {
        DemoScript.currentObject.cancelScript();
        this._closeAll();
        grok.shell.addView(new DemoView());
      }

      if (value.root.classList.contains('d4-tree-view-item')) {
        const demoFunc = DemoView.findDemoFunc(value.value.path);        
        await this.startDemoFunc(demoFunc, value.value.path);
        this.tree.root.focus();
      } else if (value.root.classList.contains('demo-app-tree-home-node')) { 
        this._initContent();
        this._closeAll();
        const view = grok.shell.addView(this);
        view.basePath = '/apps/Tutorials/Demo';
      } else {
        this.tree.root.focus();
        this.nodeView(value.text, value.value.path);
      }
    });

    this.dockPanel = grok.shell.dockManager.dock(ui.panel([
      this.searchInput.root,
      this.tree.root,
    ]), DG.DOCK_TYPE.LEFT, null, 'Categories');
    this.dockPanel.container.containerElement.classList.add('tutorials-demo-container');

    this.tree.root.classList.add('demo-app-tree-group');
  }

  private _isDockPanelInit(): boolean {
    const panelDockNode = Array.from(grok.shell.dockManager.rootNode.children)[0];
    return panelDockNode?.container.containerElement.classList.contains('tutorials-demo-container');
  }

  private _searchItem(): void {
    const foundFuncs = this.funcs.filter((func) => {
      return func.name.toLowerCase().includes(this.searchInput.value.toLowerCase()) ||
        func.func.description.toLowerCase().includes(this.searchInput.value.toLowerCase()) ||
        func.keywords.toLowerCase().includes(this.searchInput.value.toLowerCase())
    });
    
    const dom = this.tree.root.getElementsByClassName('d4-tree-view-node');

    for (let i = 0; i < dom.length; i++) {
      const item = dom[i] as HTMLElement;
      const foundFunc = foundFuncs.find((func) => func.name.toLowerCase() === item.innerText.toLowerCase());
      if (foundFunc) {
        const foundFuncPath = foundFunc.path.split('|').map((s) => s.trim());
        item.classList.remove('hidden');
        if (item.classList.contains('d4-tree-view-item')) {
          for (let i = foundFuncPath.length - 2; i >= 0; i--) {
            const currentCategory = this.tree.root.querySelector(`[data-name="${foundFuncPath[i]}"]`);
            currentCategory?.classList.remove('hidden');
          }
        }
      }
      else if (item.innerText.toLowerCase().includes(this.searchInput.value.toLowerCase())) {
        item.classList.remove('hidden');
        // if (!DIRECTIONS.includes(this.searchInput.value.toLowerCase())) {

        // }
      }
      else
        item.classList.add('hidden');
    }
  }

  private _initWindowOptions(): void {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showRibbon = true;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;

    grok.shell.windows.help.syncCurrentObject = false;
  }
}
