import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import '../../css/demo.css';
//import { _package } from '../package-test';
import { _package } from '../package';

export class DemoView extends DG.ViewBase {
  dockPanel: DG.DockNode = new DG.DockNode(undefined);
  tree: DG.TreeViewGroup = ui.tree();
  search: DG.InputBase = ui.searchInput('', '');
  
  constructor() {
    super();
    this._initDockPanel();
    this._initContent();
    this.tree.root.classList.add('demo-app-tree-group');
  }

  static findDemoFunc(demoPath: string) {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  async startDemoFunc(func: DG.Func, viewPath: string) {
    grok.shell.closeAll();
    const loadingScreen = ui.div('Loading...', 'loading');
    grok.shell.tv.root.appendChild(loadingScreen);

    await func.apply();
    loadingScreen.remove();
    grok.shell.v.path = grok.shell.v.basePath = `/apps/Tutorials/Demo/${viewPath}`;
  }

  private _initContent() {
    this.name = 'Demo app';
    this.root.appendChild(ui.divText('Select a demo from the toolbox on the left', 'demo-text'));
  }

  private nodeView(viewName: string, viewPath: string) {
    
    if (viewName === 'Viewers') {
      let root = ui.div([], 'grok-gallery-grid');
      grok.shell.closeAll();
      grok.shell.newView(viewName);
      grok.shell.v.path = grok.shell.v.basePath = `/apps/Tutorials/Demo/${viewPath}`;

      for (const f of DG.Func.find({meta: {'demoPath': null}})) {
        if (f.options[DG.FUNC_OPTIONS.DEMO_PATH].includes('Viewers')){
          const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
          const path = pathOption.split('|').map((s) => s.trim());
          const viewer = path[path.length - 1];
          const demoPath = `Viewers/${viewer}`; 

          let image = ui.image(`${_package.webRoot}images/viewers/${f.friendlyName}Img.jpg`, 0, 0);
          if (viewer === 'Globe')
            image.style.backgroundColor = 'black';

          let item = ui.card(ui.divV([
            image,
            ui.div([viewer],'tutorials-card-title'),
            ui.div([f.description], 'tutorials-card-description')
          ], 'demo-app-card'));
          item.addEventListener('click', async () => {
            let node = this.tree.items.find(node => node.text == viewer)?.root;
            node?.click();
            //const demoPath = `Viewers/${viewer}`;
            //await this.startDemoFunc(f, demoPath);
          });
          root.append(item);
        }
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
        if (f.description)
          ui.tooltip.show(f.description, event.clientX, event.clientY);
      };

      item.root.onmouseout = (_) => {
        ui.tooltip.hide();
      };
      
    }

    this.search.onChanged(() => {
      const dom = this.tree.root.getElementsByClassName('d4-tree-view-node');
      
      for (let i = 0; i < dom.length; i++) {
        const item = dom[i] as HTMLElement;
        if (item.innerText.toLowerCase().includes(this.search.value.toLowerCase())){
          item.classList.remove('hidden');
        }
        else
          item.classList.add('hidden');
      }
    });

    this.search.input.onkeyup = (event) => {
      if (event.key === 'Escape')
        this.search.fireChanged();
    };

    const closeIcon = this.search.root.getElementsByClassName('ui-input-icon-right')[0] as HTMLElement;
    closeIcon.onclick = () => {
      this.search.value = '';
      this.search.fireChanged();
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
        this.nodeView(value.text, value.text);
      }
    });

    this.dockPanel = grok.shell.dockManager.dock(ui.panel([
      this.search.root,
      this.tree.root,
    ]), 'left', null, 'Categories');
    this.dockPanel.container.containerElement.classList.add('tutorials-demo-container');

    this._initWindowOptions();

    // grok.events.onCurrentViewChanged.subscribe((view) => this.tree.root.focus());

    // TODO: if loading ended in 0.1s, then no div, if not - then div - DG.debounce, merge etc.
    // TODO: add starting demo app viewer on just up/down arrows
    // TODO: on click on viewer demo set viewer help url in property panel (func helpUrl)
    // TODO: implement search in demo - search on meta.keywords, name, description
    // TODO: add all the platform viewers to demo (make demo functions in Tutorials)

    // TODO: main viewer: 3/4
    // TODO: additional: grid - histogram - barchart (all by 33%)
    // TODO: for standard the same
    // TODO: if there empty space - add viewer/filter/etc.
    // TODO: write API for step control and example, steps are written in context panel - first priority
    // TODO: add DG.debounce
  }

  private _initWindowOptions() {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showRibbon = true;
    grok.shell.windows.showHelp = true;
    grok.shell.windows.showProperties = false;
  }
}
