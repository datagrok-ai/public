import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';

import '../../css/demo.css';


type Direction = {
  name: string,
  icon: string
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
    grok.shell.closeAll();
    ui.setUpdateIndicator(grok.shell.tv.root, true);
    grok.shell.windows.showHelp = true;

    await func.apply();
    ui.setUpdateIndicator(grok.shell.tv.root, false);

    grok.shell.v.path.includes('/apps/Tutorials/Demo') ?
      grok.shell.v.path = grok.shell.v.basePath = `/${viewPath}` :
      grok.shell.v.path = grok.shell.v.basePath = `/apps/Tutorials/Demo/${viewPath}`;
  }


  private _initContent() {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;

    this.name = 'Demo app';

    const title = ui.divText('Datagrok Platform Demo', 'demo-app-view-title');
    const description = ui.divText(`Explore Datagrok functionality features.
      Select a the category or choose the demo from the list.`, 'demo-app-view-subtitle');
    const root = ui.div([], 'demo-app-view grok-gallery-grid');

    const groups: Direction[] = [];

    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
      const path = pathOption.split('|').map((s) => s.trim());
      const categoryName = path[0];

      const item: Direction = {
        name: categoryName,
        icon: f.package.getIconUrl()
      };

      if (!groups.some((elem) => elem.name === categoryName))
        groups[groups.length] = item;
    }

    for (let i = 0; i < groups.length; i++) {
      const item = ui.card(ui.divV([
        ui.image(groups[i].icon, 80, 80),
        ui.div([groups[i].name], 'tutorials-card-title')
      ]));

      item.onclick = () => {
        let node = this.tree.items.find(node => node.text == groups[i].name)?.root;
        node?.click();
      };
      root.append(item);
    }

    this.root.append(ui.divV([
      title,
      description,
      root
    ]));
  }

  nodeView(viewName: string) {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = false;

    grok.shell.closeAll();

    const view = grok.shell.newView(viewName);
    view.basePath = '/apps/Tutorials/Demo';
    view.path = `/${viewName}`;

    const root = ui.div([], 'demo-app-group-view grok-gallery-grid');

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

        if (f.description != '')
          ui.tooltip.bind(item, f.description)

        root.append(item);
      }

      grok.shell.v.root.append(root);
    }
  }

  private _initDockPanel() {
    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
      const path = pathOption.split('|').map((s) => s.trim());
      const folder = this.tree.getOrCreateGroup(path.slice(0, path.length - 1).join(' | '));
      const item = folder.item(path[path.length - 1]);

      item.root.onmouseover = (event) => {
        const packageMessage = `Part of the ${f.package.name} package`;
        ui.tooltip.show(f.description ? ui.divV([f.description, packageMessage]) : ui.div(packageMessage),
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

    DG.debounce(this.tree.onSelectedNodeChanged, 300).subscribe(async (value) => {
      if (value.root.classList.contains('d4-tree-view-item')) {
        const categoryName = value.root.parentElement?.parentElement
          ?.getElementsByClassName('d4-tree-view-group-label')[0].innerHTML;
        const viewerName = value.text;
        const demoFunc = DemoView.findDemoFunc(`${categoryName} | ${viewerName}`);
        const demoPath = `${categoryName}/${viewerName}`;
        await this.startDemoFunc(demoFunc, demoPath);
        this.tree.root.focus();
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

    // TODO: if loading ended in 0.1s, then no div, if not - then div - DG.debounce, merge etc.
    // TODO: on click on viewer demo set viewer help url in property panel (func helpUrl)
    // TODO: implement search in demo - search on meta.keywords, name, description
    // TODO: add all the platform viewers to demo (make demo functions in Tutorials)

    // TODO: if there empty space - add viewer/filter/etc.
    // TODO: write API for step control and example, steps are written in context panel - first priority

    // TODO: add to script demo class grok.shell.windows.showPropertyPanel = true and showHelp = false
    // TODO: add GIS
  }

  private _initWindowOptions() {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showRibbon = true;
    grok.shell.windows.showHelp = true;
    grok.shell.windows.showProperties = false;
  }
}
