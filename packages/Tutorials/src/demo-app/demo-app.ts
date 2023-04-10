import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import '../../css/demo.css';


export class DemoView extends DG.ViewBase {
  dockPanel: DG.DockNode = new DG.DockNode(undefined);
  tree: DG.TreeViewGroup = ui.tree();

  constructor() {
    super();

    this._initDockPanel();
    this._initContent();
  }

  static findDemoFunc(demoPath: string) {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  startDemoFunc(func: DG.Func, viewPath: string) {
    grok.shell.closeAll();
    const loadingScreen = ui.div('Loading...', 'loading');
    grok.shell.tv.root.appendChild(loadingScreen);

    func.apply().then((_) => {
      loadingScreen.remove();
      grok.shell.tv.path = grok.shell.tv.basePath = `/apps/Tutorials/Demo/${viewPath}`;
    });
  }

  private _initContent() {
    this.name = 'Demo app';
    this.root.appendChild(ui.divText('Select a demo from the toolbox on the left', 'demo-text'));
  }

  private _initDockPanel() {
    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
      const path = pathOption.split('|').map((s) => s.trim());
      const folder = this.tree.getOrCreateGroup(path.slice(0, path.length - 1).join(' | '));
      const item = folder.item(path[path.length - 1]);

      item.root.onmousedown = (_) => {
        this.startDemoFunc(f, `${path[0]}/${path[1]}`);
      };

      item.root.onmouseover = (event) => {
        if (f.description)
          ui.tooltip.show(f.description, event.clientX, event.clientY);
      };

      item.root.onmouseout = (_) => {
        ui.tooltip.hide();
      };
    }

    this.tree.onNodeEnter.subscribe((value) => {
      if (value.root.classList.contains('d4-tree-view-item')) {
        const categoryName = value.root.parentElement?.parentElement
          ?.getElementsByClassName('d4-tree-view-group-label')[0].innerHTML;
        const viewerName = value.text;

        const demoFunc = DemoView.findDemoFunc(`${categoryName} | ${viewerName}`);
        const demoPath = `${categoryName}/${viewerName}`;
        this.startDemoFunc(demoFunc, demoPath);
        // TODO: add focus return to dock panel
      }
    });

    this.dockPanel = grok.shell.dockManager.dock(ui.div(
      [ui.searchInput('', ''), this.tree]), 'left', null, 'Categories');
    this.dockPanel.container.containerElement.classList.add('tutorials-demo-container');

    this._initWindowOptions();

    // TODO: if loading ended in 0.1s, then no div, if not - then div - DG.debounce, merge etc.
    // TODO: add starting demo app viewer on just up/down arrows
    // TODO: on click on viewer demo set viewer help url in property panel (func helpUrl)
    // TODO: implement search in demo - search on meta.keywords, name, description
    // TODO: add all the platform viewers to demo (make demo functions in Tutorials)
    // TODO: add curves to demo
  }

  private _initWindowOptions() {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showRibbon = true;
    grok.shell.windows.showHelp = true;
  }
}
