import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {sortFunctionsByHierarchy} from './utils';
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
let treeNodeY: number = 0;

export class DemoView extends DG.ViewBase {
  funcs: DemoFunc[] = [];
  subCategories: string[] = [];
  // browseView: DG.BrowseView = grok.shell.view('Browse') as DG.BrowseView;
  tree: DG.TreeViewGroup;
  DEMO_APP_PATH: string = 'apps/Tutorials/Demo';

  constructor(initVisual: boolean = true) {
    super();
    // this.browseView.showTree = true;
    // this.tree = this.browseView.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Demo');
    this.tree = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Demo');
    this._initFunctions();
    if (this.tree.items.length === 0)
      this._initTree();
    if (initVisual)
      this._initContent();
  }

  static findDemoFunc(demoPath: string): DG.Func {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  public async startDemoFunc(func: DG.Func, viewPath: string): Promise<void> {
    const path = viewPath.split('|').map((s) => s.trim()).join('/');
    const updateIndicatorRoot = (Array.from(document.querySelectorAll('#elementContent')) as HTMLElement[])
      .find((el) => el.classList.contains('d4-dock-container'))!;

    if (func.options['isDemoScript'] == 'True') {
      ui.setUpdateIndicator(updateIndicatorRoot, true);
      const pathElements = viewPath.split('|').map((s) => s.trim());
      const v = DG.View.create('demo-app-script-view');
      v.name = pathElements[pathElements.length - 1];
      v.appendAll([ui.panel([
        ui.h1(pathElements[pathElements.length - 1]),
        ui.divText(func.description),
        ui.bigButton('Start', async () => {
          try {
            await func.apply();
          } catch (e) {
            console.error(e);
          } finally {
            // this.tree.root.focus();
            this.tree.rootNode.root.focus();
          }
        })
      ])]);
      grok.shell.addView(v);
      // this.browseView.preview = v;
      ui.setUpdateIndicator(updateIndicatorRoot, false);
    } else {
      ui.setUpdateIndicator(updateIndicatorRoot, true);
      try {
        // const sub = grok.events.onViewAdded.subscribe((view) => {
        //   this.browseView.preview = view;
        //   grok.shell.v = this.browseView;
        // });
        await func.apply();
        // sub.unsubscribe();
      }
      finally { }
      if (grok.shell.tv instanceof DG.TableView)
        await grok.data.detectSemanticTypes(grok.shell.tv.dataFrame);
      ui.setUpdateIndicator(updateIndicatorRoot, false);
    }

    grok.shell.v.path = `${this.DEMO_APP_PATH}/${path.replaceAll(' ', '-')}`;
    this._setBreadcrumbsInViewName(viewPath.split('|').map((s) => s.trim()));
  }


  private _setBreadcrumbsInViewName(viewPath: string[]): void {
    const path = viewPath.includes('Home') ? viewPath : ['Home', ...viewPath];
    const breadcrumbs = ui.breadcrumbs(path);

    breadcrumbs.onPathClick.subscribe(async (value) => {
      const currentFunc = this.funcs.filter((func) => {
        return (func.name === value[value.length - 1]);
      });
      if (currentFunc.length !== 0)
        return;
      if (value.length === 1 && value[0] === 'Home') {
        // await this.browseView.setHomeView();
        return;
      }
      this.nodeView(value[value.length - 1], (value[0] === 'Home' ? value.slice(1) : value).join('/'));
    });

    // const viewNameRoot = this.browseView.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
    // if (viewNameRoot) {
    //   viewNameRoot.textContent = '';
    //   viewNameRoot.appendChild(breadcrumbs.root);
    // }
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


  public _initFunctions(): DemoFunc[] {
    const funcs = DG.Func.find({meta: {'demoPath': null}}).sort(sortFunctionsByHierarchy);

    for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
      const directionFuncs = funcs.filter((func) => {
        return (func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DEMO_APP_HIERARCHY.children[i].name);
      });
      let tempArr: string[] = [];

      for (let j = 0; j < directionFuncs.length; ++j) {
        let imgPath = `${_package.webRoot}images/demoapp/${directionFuncs[j].name}.png`;

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
      const directionFuncs = this.funcs.filter((func) => (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(name));
      //tree.group(name, null, true).root.lastChild?.appendChild(this._groupRoot(name));
      const root = this._createViewRootElement(name);
      root.classList.add('grok-gallery-grid');

      const treeGroup = tree.group(name, null, true);

      let tempArr: string[] = [];
      let groupRoot = ui.div([], 'grok-gallery-grid');

      const collection = root.querySelectorAll('.d4-item-card')

      for (let i = 0; i < collection.length; i++) {
        if (collection[i].hasAttribute('data-sub-category'))
          tempArr.push(collection[i].getAttribute('data-sub-category') as string)
        else
          groupRoot.append(collection[i] as HTMLElement)
      }

      let subCategories = [...new Set(tempArr)];
      const subCategoryFuncs = directionFuncs.filter((func) => {
        for (let i = 0; i < subCategories.length; i++)
          if (func.path.includes(subCategories[i]))
            return true;
        return false;
      });

      if (subCategories.length > 1 || (subCategories.length > 0 && subCategoryFuncs.length !== directionFuncs.length)) {
        for (let i = 0; i < subCategories.length; i++){
          const subGroupRoot = ui.div([], 'grok-gallery-grid');
          const subTreeGroup = treeGroup.group(String(subCategories[i]), null, true);
          for (let j = 0; j < collection.length; j++)
            if (collection[j].getAttribute('data-sub-category') === subCategories[i])
              subGroupRoot.append(collection[j])
          subTreeGroup.root.lastChild?.appendChild(subGroupRoot);
        }
      } else
        for (let i = 0; i < collection.length; i++)
          groupRoot.append(collection[i])

      treeGroup.root.lastChild?.appendChild(groupRoot);
    }

    const searchInput = this._createSearchInput(this.funcs, tree);
    this.root.append(ui.div([searchInput.root, tree.root], 'grok-gallery-grid'));
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
              return Promise.reject(`${_package.webRoot}images/demoapp/emptyImg.png`)
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
        const node = this.tree.items.find(node => node.text.replaceAll('amp;', '') === directionFuncs[i].name)?.root;
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

  nodeView(viewName: string, path: string): void {
    this._initWindowOptions();
    this._closeDemoScript();

    const view = DG.View.create();
    view.name = viewName;
    view.append(resultContainer);
    // this.browseView.path = `${this.DEMO_APP_PATH}/${path.replaceAll(' ', '-')}`;
    grok.shell.v.path = `${this.DEMO_APP_PATH}/${path.replaceAll(' ', '-')}`;

    const directionFuncs = this.funcs.filter((func) => (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(viewName));
    const root = this._createViewRootElement(viewName);
    root.classList.add('grok-gallery-grid');

    const tree = ui.tree();
    const treeGroup = tree.group(viewName, null, true);

    let tempArr: string[] = [];
    let groupRoot = ui.div([], 'grok-gallery-grid');

    const collection = root.querySelectorAll('.d4-item-card')

    for (let i = 0; i < collection.length; i++) {
      if (collection[i].hasAttribute('data-sub-category'))
        tempArr.push(collection[i].getAttribute('data-sub-category') as string)
      else
        groupRoot.append(collection[i] as HTMLElement)
    }

    let subCategories = [...new Set(tempArr)];
    const subCategoryFuncs = directionFuncs.filter((func) => {
      for (let i = 0; i < subCategories.length; i++)
        if (func.path.includes(subCategories[i]))
          return true;
      return false;
    });

    if (subCategories.length > 1 || (subCategories.length > 0 && subCategoryFuncs.length !== directionFuncs.length)) {
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
    } else
      for (let i = 0; i < collection.length; i++)
        groupRoot.append(collection[i])

    treeGroup.root.lastChild?.appendChild(groupRoot);

    tree.root.classList.add('demo-app-group-view');
    const searchInput = this._createSearchInput(directionFuncs, tree);
    view.root.append(ui.div([searchInput.root, tree.root], 'grok-gallery-grid'));
    this._setBreadcrumbsInViewName(path.split('/').map((s) => s.trim()));
    // this.tree.root.focus();
    this.tree.rootNode.root.focus();
    // this.browseView.preview = view;
  }

  private _createSearchInput(directionFuncs: DemoFunc[], tree: DG.TreeViewGroup): DG.InputBase<string> {
    const searchInput = ui.input.search('', {
      onValueChanged: (value) => {
        const foundFuncs = directionFuncs.filter((func) => {
          return func.name.toLowerCase().includes(value.toLowerCase()) ||
              func.func.description.toLowerCase().includes(value.toLowerCase()) ||
              func.keywords.toLowerCase().includes(value.toLowerCase())
        });
        const cards = tree.root.querySelectorAll('.d4-item-card.ui-div');
        for (let i = 0; i < cards.length; i++) {
          const cardTitle = cards[i].querySelector('.tutorials-card-title.ui-div') as HTMLElement;
          if (foundFuncs.some((func) => func.name.toLowerCase() === cardTitle.innerText.toLowerCase()))
            cards[i].classList.remove('hidden');
          else
            cards[i].classList.add('hidden');
        }
      },
      elementOptions: {classes: 'demo-app-search-input'}
    });

    searchInput.input.onkeyup = (event) => {
      if (event.key === 'Escape')
        searchInput.fireChanged();
    }
    const closeIcon = searchInput.root.getElementsByClassName('ui-input-icon-right')[0] as HTMLElement;
    closeIcon.onclick = () => {
      searchInput.value = '';
      searchInput.fireChanged();
    };

    return searchInput;
  }

  private _initTree(): void {
    for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
      const directionFuncs = this.funcs.filter((func) => {
        return (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DEMO_APP_HIERARCHY.children[i].name);
      });

      for (let j = 0; j < directionFuncs.length; ++j) {
        const path = directionFuncs[j].path.split('|').map((s) => s.trim());

        if (path.length > 2) {
          let groupPath = path[0];
          let treePath = this.tree.getOrCreateGroup(path[0], {path: groupPath}, false);
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
          const folder = this.tree.getOrCreateGroup(directionFuncs[j].category, {path: path[0]}, false);
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
    this.tree.root.classList.add('demo-app-tree-group');

    DG.debounce(this.tree.rootNode.onSelectedNodeChanged, 300).subscribe(async (value) => {
      if (!value || !this.tree.root.contains(value.root) || value.text === 'Demo')
        return;

      const panelRoot = this.tree.rootNode.root.parentElement!;
      treeNodeY = panelRoot.scrollTop!;

      if (DemoScript.currentObject) {
        DemoScript.currentObject.cancelScript();
        // this.browseView.preview = new DemoView() as unknown as DG.View;
        grok.shell.v = new DemoView() as unknown as DG.View;
        panelRoot.scrollTo(0, treeNodeY);
      }

      if (value.root.classList.contains('d4-tree-view-item')) {
        const demoFunc = DemoView.findDemoFunc(value.value.path);
        await this.startDemoFunc(demoFunc, value.value.path);
        // this.tree.root.focus();
        this.tree.rootNode.root.focus();
      } else {
        // this.tree.root.focus();
        this.tree.rootNode.root.focus();
        this.nodeView(value.text, value.value.path);
      }

      panelRoot.scrollTo(0, treeNodeY);
    });
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
