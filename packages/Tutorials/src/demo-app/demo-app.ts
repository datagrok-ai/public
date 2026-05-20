import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {sortFunctionsByHierarchy} from './utils';
import {DEMO_APP_HIERARCHY} from './const';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import '../../css/demo.css';

/** Structural shape of what the demo app needs from a `func` (real DG.Func or project shim). */
export interface DemoEntryFunc {
  name: string;
  description: string;
  options: {[key: string]: any};
  package: {name: string};
  apply(): Promise<any>;
}

export type DemoFunc = {
  name: string;
  func: DemoEntryFunc,
  category: string;
  path: string;
  keywords: string;
  imagePath: string;
};

const resultContainer = ui.div([], 'hidden');
let treeNodeY: number = 0;

export class DemoView extends DG.ViewBase {
  funcs: DemoFunc[] = [];
  tree: DG.TreeViewGroup;
  currentView: DG.ViewBase | null = null;
  private focusSub: { unsubscribe(): void } | null = null;
  private _focusGuardCleanup: (() => void) | null = null;
  DEMO_APP_PATH: string = 'apps/Tutorials/Demo';

  private _addedProjectIds = new Set<string>();
  public projectsReady: Promise<void>;

  constructor(initVisual: boolean = true) {
    super();
    this.tree = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Demo');
    this._initFunctions();
    if (this.tree.items.length === 0)
      this._initTree();
    if (initVisual)
      this._initContent();
    this.projectsReady = this._initProjects();
  }

  private async _initProjects(): Promise<void> {
    const found = await grok.dapi.projects.filter('entityMetaParams.name = "demoPath"').list({pageSize: 50});
    if (found.length === 0)
      return;
    const projects = await Promise.all(found.map((p) => grok.dapi.projects.find(p.id)));
    let added = 0;
    for (const p of projects) {
      if (this._addedProjectIds.has(p.id))
        continue;
      const path = p.meta?.demoPath as string | undefined;
      if (!path)
        continue;
      const dir = DEMO_APP_HIERARCHY.children.find((d) => path.includes(d.name));
      if (!dir)
        continue;
      const parts = path.split('|').map((s) => s.trim());
      const entry: DemoFunc = {
        name: parts[parts.length - 1],
        func: this._projectAsFunc(p, path),
        category: dir.name, path, keywords: '',
        imagePath: p.pictureUrl,
      };
      this.funcs.push(entry);
      this._addProjectToTree(entry);
      this._addedProjectIds.add(p.id);
      added++;
    }
    if (added > 0) {
      this.funcs.sort((a, b) => sortFunctionsByHierarchy(a.func as DG.Func, b.func as DG.Func));
      if (this.root)
        this._initContent();
    }
  }

  private _projectAsFunc(p: DG.Project, demoPath: string): DemoEntryFunc {
    return {
      name: p.name,
      description: p.description ?? '',
      options: {[DG.FUNC_OPTIONS.DEMO_PATH]: demoPath},
      package: {name: 'Project'},
      apply: () => p.open(),
    };
  }

  private _addProjectToTree(entry: DemoFunc): void {
    const path = entry.path.split('|').map((s) => s.trim());
    let group: DG.TreeViewGroup = this.tree;
    let parentName: string;
    if (path.length > 2) {
      let gp = path[0];
      group = this.tree.getOrCreateGroup(path[0], {path: gp}, false);
      for (let i = 1; i < path.length - 1; i++) {
        gp += `/${path[i]}`;
        group = group.getOrCreateGroup(path[i], {path: gp}, false);
      }
      parentName = path[path.length - 2];
    }
    else {
      group = this.tree.getOrCreateGroup(entry.category, {path: path[0]}, false);
      parentName = entry.category;
    }
    if (group.items.some((n) => n.text === entry.name))
      return;
    const node = group.item(entry.name, {path: entry.path});
    node.root.onmouseover = (event) => {
      ui.tooltip.show(entry.func.description ?
        ui.divV([entry.func.description, ui.element('br'), 'Datagrok project']) :
        ui.div('Datagrok project'), event.clientX, event.clientY);
    };
    node.root.onmouseout = (_) => ui.tooltip.hide();
    const targetIndex = this._hierarchyInsertIndex(group, parentName, entry.name);
    if (targetIndex !== null && targetIndex !== group.children.length - 1) {
      node.remove();
      group.addNode(node, targetIndex);
    }
  }

  /** Index where `name` should sit in `group` based on DEMO_APP_HIERARCHY order, or null if not orderable. */
  private _hierarchyInsertIndex(group: DG.TreeViewGroup, parentName: string, name: string): number | null {
    const children = this._findHierarchyChildren(parentName);
    if (children === null)
      return null;
    const targetPos = children.findIndex((c) => c.name === name);
    if (targetPos < 0)
      return null;
    const order = new Map(children.map((c, i) => [c.name, i] as [string, number]));
    const siblings = group.children.filter((n) => n.text !== name);
    let index = 0;
    for (const sib of siblings) {
      const sibPos = order.get(sib.text) ?? Infinity;
      if (sibPos < targetPos)
        index++;
      else
        break;
    }
    return index;
  }

  private _findHierarchyChildren(name: string, list: {name: string; children?: any[]}[] = DEMO_APP_HIERARCHY.children): {name: string}[] | null {
    for (const c of list) {
      if (c.name === name && c.children)
        return c.children;
      if (c.children) {
        const r = this._findHierarchyChildren(name, c.children);
        if (r)
          return r;
      }
    }
    return null;
  }

  static findDemoFunc(demoPath: string): DG.Func {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  public async startDemoFunc(func: DemoEntryFunc, viewPath: string): Promise<void> {
    const splitViewPath = viewPath.split('|');
    const path = splitViewPath.map((s) => s.trim()).join('/');
    const dockNodes = Array.from(grok.shell.dockManager.rootNode.children)
      .filter((dockNode) => dockNode.container.containerElement.getElementsByClassName('tab-content').length > 0);
    if (dockNodes.length === 0)
      return;
    const dockNode = dockNodes.filter((node) => node.container.containerElement.classList.contains('document-manager') ||
      node.container.containerElement.getElementsByClassName('document-manager').length > 0)[0];
    const activeRoot = grok.shell.v?.root ?? this.currentView?.root;
    let activeTabContent: HTMLElement | null = null;
    if (activeRoot) {
      let el: HTMLElement | null = activeRoot.parentElement;
      while (el && !el.classList.contains('tab-content'))
        el = el.parentElement;
      activeTabContent = el;
    }
    const tabContent = (activeTabContent ?? dockNode.container.containerElement.getElementsByClassName('tab-content')[0]) as HTMLElement;
    const updateIndicatorRoot = (tabContent.parentElement ?? tabContent) as HTMLElement;

    this.currentView = null;

    const viewsBefore = new Set<DG.View>(Array.from(grok.shell.views));

    if (func.options['isDemoScript'] == 'True') {
      updateIndicatorRoot.classList.add('demo-app-loading');
      ui.setUpdateIndicator(updateIndicatorRoot, true);
      try {
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
              this.tree.rootNode.root.focus();
            }
          })
        ])]);
        grok.shell.addView(v);
        this._closePrevDemoViews();
        this.currentView = v;
        v.path = `${this.DEMO_APP_PATH}/${path.replaceAll(' ', '-')}`;
      } finally {
        this._tagNewDemoViews(viewsBefore, func.name);
        ui.setUpdateIndicator(updateIndicatorRoot, false);
        updateIndicatorRoot.classList.remove('demo-app-loading');
      }
    } else {
      updateIndicatorRoot.classList.add('demo-app-loading');
      ui.setUpdateIndicator(updateIndicatorRoot, true);
      const prevAutoShowToolbox = grok.shell.windows.autoShowToolbox;
      grok.shell.windows.autoShowToolbox = false;
      const viewBeforeApply = grok.shell.v;
      try {
        await func.apply();
        this._closePrevDemoViews();
        this.currentView = grok.shell.v;
        if (grok.shell.tv instanceof DG.TableView)
          await grok.data.detectSemanticTypes(grok.shell.tv.dataFrame);
      } finally {
        this._tagNewDemoViews(viewsBefore, func.name);
        grok.shell.windows.autoShowToolbox = prevAutoShowToolbox;
        ui.setUpdateIndicator(updateIndicatorRoot, false);
        updateIndicatorRoot.classList.remove('demo-app-loading');
      }
      this.tree.rootNode.root.focus();
      this._guardTreeFocus();
      if (grok.shell.v !== viewBeforeApply) {
        grok.shell.v.name = splitViewPath[splitViewPath.length - 1].trim();
        grok.shell.v.path = `${this.DEMO_APP_PATH}/${path.replaceAll(' ', '-')}`;
        this._setBreadcrumbsInViewName(viewPath.split('|').map((s) => s.trim()));
      }
    }
  }

  private _closePrevDemoViews(): void {
    const toClose = Array.from(grok.shell.views).filter((v) => v.temp?.['demoApp']);
    for (const v of toClose) {
      try {
        v.close();
      } catch (_) {}
    }
  }

  private _tagNewDemoViews(viewsBefore: Set<DG.View>, funcName: string): void {
    for (const v of grok.shell.views) {
      if (!viewsBefore.has(v))
        v.temp['demoApp'] = funcName;
    }
  }

  private _setBreadcrumbsInViewName(viewPath: string[], view?: DG.View): void {
    const usedView = view ? view === grok.shell.v ? view : null : grok.shell.v;
    const path = ['Home', 'Demo', ...viewPath.filter((v) => v !== 'Home' && v !== 'Demo')];
    const breadcrumbs = ui.breadcrumbs(path);

    breadcrumbs.onPathClick.subscribe(async (value) => {
      const actualItem = value[value.length - 1];
      if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
        return;
      this.tree.currentItem = actualItem === 'Demo' ? this.tree : this.tree.items.find((item) => item.text === actualItem)!;
    });

    if (usedView) {
      if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') { // integrate it to the actual breadcrumbs element
        const homeIcon = ui.iconFA('home', () => {
          grok.shell.v.close();
          grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME);
        }, 'Home');
        homeIcon.classList.add('demo-breadcrumbs-home-element');
        breadcrumbs.root.firstElementChild!.replaceWith(homeIcon);
      }
      const viewNameRoot = usedView.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
      if (viewNameRoot) {
        viewNameRoot.textContent = '';
        viewNameRoot.appendChild(breadcrumbs.root);
      }
    }
  }

  private _closeDemoScript(): void {
    const scriptDockNodes = Array.from(Array.from(grok.shell.dockManager.rootNode.children)[0].children)
      .filter((dockNode) => dockNode.container.containerElement.classList.contains('tutorials-demo-script-container'));
    scriptDockNodes.forEach((dockNode) => grok.shell.dockManager.close(dockNode));
  }

  public _initFunctions(): DemoFunc[] {
    const funcs = DG.Func.find({meta: {'demoPath': null}}).sort(sortFunctionsByHierarchy);

    for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
      const directionFuncs = funcs.filter((func) => {
        return (func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DEMO_APP_HIERARCHY.children[i].name);
      });

      for (let j = 0; j < directionFuncs.length; ++j) {
        let imgPath = `${_package.webRoot}images/demoapp/${directionFuncs[j].name}.png`;

        const path = directionFuncs[j].options[DG.FUNC_OPTIONS.DEMO_PATH] as string;
        const pathArray = path.split('|').map((s) => s.trim());

        this.funcs[this.funcs.length] = {
          name: pathArray[pathArray.length - 1],
          func: directionFuncs[j],
          category: DEMO_APP_HIERARCHY.children[i].name,
          path: path,
          keywords: (directionFuncs[j].options['keywords'] !== undefined) ? directionFuncs[j].options['keywords'] as string : '',
          imagePath: imgPath,
        };
      }
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
    this.root.append(ui.div([searchInput.root, tree.root], 'grok-gallery-grid grok-gallery-grid-view-demo-app'));
    this.subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
      if (grok.shell.v?.name === 'Demo app')
        this._setBreadcrumbsInViewName([]);
    }));
  }

  private _createViewRootElement(viewOrGroupName: string): HTMLDivElement {
    const root = ui.div([]);

    const directionFuncs = this.funcs.filter((func) => {
      return (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(viewOrGroupName);
    });

    for (let i = 0; i < directionFuncs.length; i++) {

      const path = directionFuncs[i].path.split('|').map((s) => s.trim());

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

      const packageMessage = directionFuncs[i].func.package.name === 'Project' ? 'Datagrok project' :
        `Part of the ${directionFuncs[i].func.package.name === 'Tutorials' ?
          'platform core' : `${directionFuncs[i].func.package.name} package`}`;
      ui.tooltip.bind(item, () => directionFuncs[i].func.description ?
        ui.divV([directionFuncs[i].func.description, ui.element('br'), packageMessage]) : ui.div(packageMessage));

      root.append(item);
    }

    return root;
  }

  nodeView(viewName: string, path: string): void {
    this._closeDemoScript();
    this.currentView = null;

    const view = DG.View.create();
    view.name = viewName;
    view.append(resultContainer);
    view.path = `${this.DEMO_APP_PATH}/${path.replaceAll(' ', '-')}`;

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
    view.root.append(ui.div([searchInput.root, tree.root], 'grok-gallery-grid grok-gallery-grid-view-demo-app'));
    grok.shell.addView(view);
    this._closePrevDemoViews();
    view.temp['demoApp'] = `_categoryView:${viewName}`;
    this.currentView = view;
    this._setBreadcrumbsInViewName(path.split('/').map((s) => s.trim()));
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
          item.root.onmouseout = (_) => { ui.tooltip.hide(); };
        }
      }
    }
    this.tree.root.classList.add('demo-app-tree-group');

    DG.debounce(this.tree.rootNode.onSelectedNodeChanged, 300).subscribe(async (value) => {
      if (!value || !this.tree.root.contains(value.root) || value.text === 'Demo')
        return;

      this._closeDemoScript();
      const panelRoot = this.tree.rootNode.root.parentElement!;
      treeNodeY = panelRoot.scrollTop!;

      if (DemoScript.currentObject) {
        DemoScript.currentObject.cancelScript();
        grok.shell.v = new DemoView() as unknown as DG.View;
        panelRoot.scrollTo(0, treeNodeY);
      }

      this.focusSub?.unsubscribe();
      if (value.root.classList.contains('d4-tree-view-item')) {
        this.focusSub = grok.events.onCurrentViewChanged.subscribe(() => {
          this.focusSub?.unsubscribe();
          this.focusSub = null;
          this._initWindowOptions();
          this.tree.rootNode.root.focus();
        });
        const demoFunc = this.funcs.find((f) => f.path === value.value.path)?.func ?? DemoView.findDemoFunc(value.value.path);
        await this.startDemoFunc(demoFunc, value.value.path);
      } else {
        this.focusSub = grok.events.onCurrentViewChanged.subscribe(() => {
          this.focusSub?.unsubscribe();
          this.focusSub = null;
          this._initWindowOptions();
          this.tree.rootNode.root.focus();
          this._guardTreeFocus();
        });
        this.nodeView(value.text, value.value.path);
      }
      this.close();

      panelRoot.scrollTo(0, treeNodeY);
    });
  }

  /** Redirects focus back to the browse tree if anything steals it (e.g. grid auto-focus
   *  after data/renderers load). Cleared on the first mousedown or after a timeout,
   *  whichever comes first. */
  private _guardTreeFocus(): void {
    this._focusGuardCleanup?.();
    const onFocusIn = (e: FocusEvent) => {
      if (!this.tree.rootNode.root.contains(e.target as Node))
        this.tree.rootNode.root.focus();
    };
    const clear = () => {
      document.removeEventListener('focusin', onFocusIn, true);
      document.removeEventListener('mousedown', clear, true);
      clearTimeout(timer);
      this._focusGuardCleanup = null;
    };
    const timer = setTimeout(clear, 1500);
    this._focusGuardCleanup = clear;
    document.addEventListener('focusin', onFocusIn, true);
    document.addEventListener('mousedown', clear, true);
  }

  private _initWindowOptions(): void {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showRibbon = true;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;
    grok.shell.windows.help.syncCurrentObject = false;
  }
}
