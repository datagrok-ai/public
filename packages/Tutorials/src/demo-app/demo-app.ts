import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {sortFunctionsByHierarchy, getParentCategoryName} from './utils';
import {DEMO_APP_HIERARCHY} from './const';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import '../../css/demo.css';

export type DemoFunc = {
  name: string;
  func: DG.Func,
  category: string;
  path: string;
  keywords: string;
  imagePath: string;
};

const resultContainer = ui.div([], 'hidden');
const dockRoot = ui.panel([]);
let treeNodeY: number = 0;

export class DemoView extends DG.ViewBase {
  dockPanel: DG.DockNode = new DG.DockNode(undefined);
  tree: DG.TreeViewGroup = ui.tree();
  searchInput: DG.InputBase = ui.input.search('', {value: ''});
  funcs: DemoFunc[] = [];
  subCategories: string[] = [];
  isInBrowse: boolean = false;
  browseView?: DG.BrowseView;

  constructor(isInBrowse: boolean = false, browseView?: DG.BrowseView) {
    super();
    this.isInBrowse = isInBrowse;
    this.browseView = browseView;
    this._initFunctions();
    if (!isInBrowse)
      this._initDockPanel();
    this._initContent();
  }

  static findDemoFunc(demoPath: string): DG.Func {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  public async startDemoFunc(func: DG.Func, viewPath: string): Promise<void> {
    const path = viewPath.split('|').map((s) => s.trim()).join('/');

    if (!this.isInBrowse)
      this._closeAll();

    // TODO: change to getElementById('elementContent') when it is fixed
    const updateIndicatorRoot = this.isInBrowse ?
      document.querySelector('.layout-dockarea .view-tabs .dock-container.dock-container-fill > .tab-host > .tab-content > .d4-dock-container > .grok-view.grok-view-browse > .d4-root > .splitter-container-row > .dock-container.dock-container-fill > .tab-host > .tab-content')! as HTMLElement :
      document.querySelector('.layout-dockarea .view-tabs .dock-container.dock-container-fill > .tab-host > .tab-content')! as HTMLElement;

    if (func.options['isDemoScript'] == 'True') {
      ui.setUpdateIndicator(updateIndicatorRoot, true);
      const pathElements = viewPath.split('|').map((s) => s.trim());
      if (this.isInBrowse && this.browseView) {
        const v = DG.View.create('demo-app-script-view');
        v.name = pathElements[pathElements.length - 1];
        v.appendAll([ui.panel([
          ui.h1(pathElements[pathElements.length - 1]),
          ui.divText(func.description),
          ui.bigButton('Start', async () => {
            grok.shell.isInDemo = true;
            try {
              await func.apply();
            } catch (e) {
              grok.shell.isInDemo = false;
              console.error(e);
            } finally {
              this.tree.root.focus();
            }
          })
        ])]);
        this.browseView.preview = v;
      }
      else {
        grok.shell.newView(pathElements[pathElements.length - 1], [ ui.panel([
            ui.h1(pathElements[pathElements.length - 1]),
            ui.divText(func.description),
            ui.bigButton('Start', async () => { await func.apply() })
          ], 'demo-app-script-view')
        ]);
      }
      ui.setUpdateIndicator(updateIndicatorRoot, false);
    } else {
      ui.setUpdateIndicator(updateIndicatorRoot, true);
      if (!this.isInBrowse)
        await func.apply();
      else {
        grok.shell.isInDemo = true;
        try {
          await func.apply();
        }
        finally {
          grok.shell.isInDemo = false;
        }
        if (this.browseView?.preview instanceof DG.TableView)
          await grok.data.detectSemanticTypes(this.browseView.preview.dataFrame);
      }
      ui.setUpdateIndicator(updateIndicatorRoot, false);
    }
    grok.shell.v.path = grok.shell.v.basePath = '';
    if (grok.shell.v.basePath.includes('/apps/Tutorials/Demo')) {
      grok.shell.v.path = `/${path}`;
    }
    else {
      grok.shell.v.basePath = '/apps/Tutorials/Demo';
      grok.shell.v.path = `/${path}`;
    }
    // temporary fix, change after all demo scripts are using meta.isDemoScript
    if ((func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes('Visualization'))    
      grok.shell.v.path = grok.shell.v.basePath = `/apps/Tutorials/Demo/${path}`
    
    this._setBreadcrumbsInViewName(viewPath.split('|').map((s) => s.trim()), grok.shell.v);
  }


  private _setBreadcrumbsInViewName(viewPath: string[], view: DG.ViewBase): void {
    const breadcrumbs = ui.breadcrumbs(viewPath);

    breadcrumbs.onPathClick.subscribe((value) => {
      const currentFunc = this.funcs.filter((func) => {
        return (func.name === value[value.length - 1]);
      });
      if (currentFunc.length !== 0)
        return;
      const v = this.nodeView(value[value.length - 1], value.join('/'));
      this.isInBrowse && this.browseView ? this.browseView.preview = v : grok.shell.addView(v);
    });

    const viewNameRoot = view.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
    if (viewNameRoot) {
      viewNameRoot.textContent = '';
      viewNameRoot.appendChild(breadcrumbs.root);
    }
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


  public _initFunctions(): DemoFunc[] {
    const funcs = DG.Func.find({meta: {'demoPath': null}}).sort(sortFunctionsByHierarchy);

    for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
      const directionFuncs = funcs.filter((func) => {
        return (func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DEMO_APP_HIERARCHY.children[i].name);
      });
      let tempArr: string[] = [];

      for (let j = 0; j < directionFuncs.length; ++j) {
        let imgPath = `${_package.webRoot}images/demoapp/${directionFuncs[j].name}.jpg`;

        const path = directionFuncs[j].options[DG.FUNC_OPTIONS.DEMO_PATH] as string;
        const pathArray = path.split('|').map((s) => s.trim());
        
        if (pathArray.length > 2) {
          tempArr.push(pathArray[1]);
        }  

        this.funcs[this.funcs.length] = {
          name: pathArray[pathArray.length - 1],
          func: directionFuncs[j],
          category: DEMO_APP_HIERARCHY.children[i].name,
          path: path,
          keywords: (directionFuncs[j].options['keywords'] !== undefined) ? directionFuncs[j].options['keywords'] as string : '',
          imagePath: imgPath,
        };
      }
      this.subCategories = [...new Set(tempArr)];
    }

    return this.funcs;
  }

  public _initContent(): void {
    this._initWindowOptions();

    this.root.innerHTML = '';
    this.root.append(resultContainer);
    this.name = 'Demo app';

    const tree = ui.tree();
    tree.root.classList.add('demo-app-group-view');
    
    for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
      const name = DEMO_APP_HIERARCHY.children[i].name;
      //tree.group(name, null, true).root.lastChild?.appendChild(this._groupRoot(name));
      const root = this._createViewRootElement(name);
      root.classList.add('grok-gallery-grid');

      const treeGroup = tree.group(name, null, true);
      
      let tempArr: string[] = [];
      let groupRoot = ui.div([], 'grok-gallery-grid');

      const collection = root.querySelectorAll('.d4-item-card')

      for (let i = 0; i < collection.length; i++) {
        if (collection[i].hasAttribute('data-sub-category')) {
          tempArr.push(collection[i].getAttribute('data-sub-category') as string)
        }
        else {
          groupRoot.append(collection[i] as HTMLElement)
        }

      }

      treeGroup.root.lastChild?.appendChild(groupRoot);
      let subCategories = [...new Set(tempArr)];
      
      if (subCategories.length > 1) {
        for (let i = 0; i < subCategories.length; i++){
          const subGroupRoot = ui.div([], 'grok-gallery-grid');
          const subTreeGroup = treeGroup.group(String(subCategories[i]), null, true);
          for (let j = 0; j < collection.length; j++) {
            if (collection[j].getAttribute('data-sub-category') === subCategories[i]) {
              subGroupRoot.append(collection[j])
            }
          }
          subTreeGroup.root.lastChild?.appendChild(subGroupRoot);
        }
      } else {
        for (let i = 0; i < collection.length; i++) {
          groupRoot.append(collection[i])
        }
      }
    }

    this.root.append(ui.div([tree.root], 'grok-gallery-grid'));
  }

  private _createViewRootElement(viewOrGroupName: string): HTMLDivElement {
    const root = ui.div([]);

    const directionFuncs = this.funcs.filter((func) => {
      return (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(viewOrGroupName);
    });
    
    for (let i = 0; i < directionFuncs.length; i++) {

      const path = directionFuncs[i].path.split('|').map((s) => s.trim());
      //const img = ui.div('', 'ui-image');
      
      const img = ui.div([ui.wait(async () => {
        let root = ui.div('','img');
        root.className = 'ui-image';
        await fetch(`${directionFuncs[i].imagePath}`)
          .then(response => {
            if (response.ok) {
              return Promise.resolve(response.url)
            } else if(response.status === 404) {
              return Promise.reject(`${_package.webRoot}images/demoapp/emptyImg.jpg`)
            }
          })
          .then((data) => root.style.backgroundImage = `url(${data})`)
          .catch((data) => root.style.backgroundImage = `url(${data})`);
        return root;
        })
      ]);

      let item = ui.card(ui.divV([
        img,
        ui.div([directionFuncs[i].name], 'tutorials-card-title'),
        ui.div([directionFuncs[i].func.description], 'tutorials-card-description')
      ], 'demo-app-card'));
      
      if ( path.length > 2 ){
        item.setAttribute('data-category', path[0]);
        item.setAttribute('data-sub-category', path[1]);
      } else {
        item.setAttribute('data-category', path[0]);
      }

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

  // TODO: add demoScript node to class

  nodeView(viewName: string, path: string): DG.View {
    this._initWindowOptions();
    if (this.isInBrowse)
      this._closeDemoScript();
    else
      this._closeAll();

    const view = DG.View.create();
    view.name = viewName;
    view.append(resultContainer);
    view.basePath = '/apps/Tutorials/Demo';
    view.path = `/${path}`;

    const root = this._createViewRootElement(viewName);
    root.classList.add('grok-gallery-grid');

    const tree = ui.tree();
    const treeGroup = tree.group(viewName, null, true);
    
    let tempArr: string[] = [];
    let groupRoot = ui.div([], 'grok-gallery-grid');

    const collection = root.querySelectorAll('.d4-item-card')

    for (let i = 0; i < collection.length; i++) {
      if (collection[i].hasAttribute('data-sub-category')) {
        tempArr.push(collection[i].getAttribute('data-sub-category') as string)
      }
      else {
        groupRoot.append(collection[i] as HTMLElement)
      }

    }

    treeGroup.root.lastChild?.appendChild(groupRoot);
    let subCategories = [...new Set(tempArr)];

    if (subCategories.length > 1) {
      for (let i = 0; i < subCategories.length; i++){
        const subGroupRoot = ui.div([], 'grok-gallery-grid');
        const subTreeGroup = treeGroup.group(String(subCategories[i]), null, true);
        for (let j = 0; j < collection.length; j++) {
          if (collection[j].getAttribute('data-sub-category') === subCategories[i]) {
            subGroupRoot.append(collection[j])
          }
        }
        subTreeGroup.root.lastChild?.appendChild(subGroupRoot);
      }
    } else {
      for (let i = 0; i < collection.length; i++) {
        groupRoot.append(collection[i])
      }
    }

    tree.root.classList.add('demo-app-group-view');
    view.root.append(ui.div([tree.root], 'grok-gallery-grid'));

    this._setBreadcrumbsInViewName(path.split('/').map((s) => s.trim()), view);

    this.tree.root.focus();

    return view;
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

    dockRoot.innerHTML = '';

    this._createHomeNode();

    this._initTree();

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

    dockRoot.append(this.searchInput.root);
    dockRoot.append(this.tree.root);
    this.dockPanel = grok.shell.dockManager.dock(dockRoot, DG.DOCK_TYPE.LEFT, null, 'Categories', 0.2);
    this.dockPanel.container.containerElement.classList.add('tutorials-demo-container');
  }

  public _initTree(tree?: DG.TreeViewGroup): void {
    for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
      const directionFuncs = this.funcs.filter((func) => {
        return (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DEMO_APP_HIERARCHY.children[i].name);
      });

      for (let j = 0; j < directionFuncs.length; ++j) {
        const path = directionFuncs[j].path.split('|').map((s) => s.trim());

        if (path.length > 2) {
          let groupPath = path[0];
          let treePath = (tree ? tree : this.tree).getOrCreateGroup(path[0], {path: groupPath}, false);
          (treePath.root.firstElementChild as HTMLElement).dataset.name = path[0];
          for (let i = 1; i < path.length - 1; i++) {
            groupPath += `/${path[i]}`;
            treePath = treePath.getOrCreateGroup(path[i], {path: groupPath}, false);
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
          const folder = (tree ? tree : this.tree).getOrCreateGroup(directionFuncs[j].category, {path: path[0]}, false);
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
    (tree ? tree : this.tree).root.classList.add('demo-app-tree-group');

    DG.debounce((tree ? tree : this.tree).rootNode.onSelectedNodeChanged, 300).subscribe(async (value) => {
      if (!this.tree.root.contains(value.root))
        return;

      const panelRoot = this.isInBrowse ? this.tree.rootNode.root.parentElement! : dockRoot.parentElement!;
      treeNodeY = panelRoot.scrollTop!;

      if (DemoScript.currentObject) {
        DemoScript.currentObject.cancelScript();
        if (this.isInBrowse && this.browseView)
          this.browseView.preview = new DemoView() as unknown as DG.View;
        else {
          this._closeAll();
          grok.shell.addView(new DemoView());
        }
        panelRoot.scrollTo(0, treeNodeY);
      }

      if (value.root.classList.contains('d4-tree-view-item')) {
        const demoFunc = DemoView.findDemoFunc(value.value.path);
        await this.startDemoFunc(demoFunc, value.value.path);
        (tree ? tree : this.tree).root.focus();
      } else if (value.root.classList.contains('demo-app-tree-home-node')) {
        this._initContent();
        this._closeAll();
        const view = grok.shell.addView(this);
        view.basePath = '/apps/Tutorials/Demo';
      } else {
        this.tree.root.focus();
        const v = this.nodeView(value.text, value.value.path);
        this.isInBrowse && this.browseView ? this.browseView.preview = v : grok.shell.addView(v);
      }

      panelRoot.scrollTo(0, treeNodeY);
    });

    if (this.isInBrowse && tree)
      this.tree = tree;
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
        let parentCategoryName = getParentCategoryName(item.innerText);
        while (parentCategoryName !== '') {
          const parentCategory = this.tree.root.querySelector(`[data-name="${parentCategoryName}"]`);
          parentCategory?.classList.remove('hidden');

          parentCategoryName = getParentCategoryName(parentCategoryName);
        }
      }
      else
        item.classList.add('hidden');
    }

    for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
      for (let j = 0; j < DEMO_APP_HIERARCHY.children[i].children.length; j++){
        if (grok.shell.v.path === `/apps/Tutorials/Demo/${DEMO_APP_HIERARCHY.children[i].name}` || 
        grok.shell.v.path === `/apps/Tutorials/Demo` ||
        grok.shell.v.path === `/apps/Tutorials/Demo/${DEMO_APP_HIERARCHY.children[i].name}/${DEMO_APP_HIERARCHY.children[i].children[j].name}` 
        ){
          grok.shell.v.root.lastElementChild?.classList.add('hidden');
          
          const root = ui.div([]);
          for (let i = 0; i < foundFuncs.length; i++) {
            const img = ui.div('', 'ui-image');
            img.style.backgroundImage = `url(${foundFuncs[i].imagePath}`;

            let item = ui.card(ui.divV([
              img,
              ui.div([foundFuncs[i].name], 'tutorials-card-title'),
              ui.div([foundFuncs[i].func.description], 'tutorials-card-description')
            ], 'demo-app-card'));

            item.onclick = () => {
              const node = this.tree.items.find(node => node.text === foundFuncs[i].name)?.root;
              node?.click();
            };

            const packageMessage = `Part of the ${foundFuncs[i].func.package.name === 'Tutorials' ?
              'platform core' : `${foundFuncs[i].func.package.name} package`}`;
            ui.tooltip.bind(item, () => foundFuncs[i].func.description ?
              ui.divV([foundFuncs[i].func.description, ui.element('br'), packageMessage]) : ui.div(packageMessage));

            root.append(item);
          }

          root.classList.add('grok-gallery-grid');

          const tree = ui.tree();
          const treeGroup = tree.group(`${foundFuncs.length} results found`, null, true);

          treeGroup.root.lastChild?.appendChild(root);
          tree.root.classList.add('demo-app-group-view');
          resultContainer.innerHTML = '';
          resultContainer.classList.remove('hidden');
          resultContainer.append(ui.div([tree.root], 'grok-gallery-grid'));
        }
      }
         
    }

    if (this.searchInput.value === ''){
        resultContainer.innerHTML = '';
        resultContainer.classList.add('hidden');
        grok.shell.v.root.lastElementChild?.classList.remove('hidden');
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

function imageExists(image_url: string){
  var http = new XMLHttpRequest();

  http.open('HEAD', image_url, false);
  http.send();

  return http.status != 404;
}