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

  startDemoFunc(func: DG.Func) {
    grok.shell.closeAll();
    const loadingScreen = ui.div('Loading...', 'loading');
    grok.shell.tv.root.appendChild(loadingScreen);

    func.apply().then((_) => {loadingScreen.remove();});
  }

  _initContent() {
    this.root.appendChild(ui.divText('Select a demo from the toolbox on the left', 'demo-text'));
  }

  _initDockPanel() {
    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
      const path = pathOption.split('|').map((s) => s.trim());
      const folder = this.tree.getOrCreateGroup(path.slice(0, path.length - 1).join(' | '));
      const item = folder.item(path[path.length - 1]);

      item.root.onmousedown = (_) => {
        this.startDemoFunc(f);
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

        this.startDemoFunc(DemoView.findDemoFunc(`${categoryName} | ${viewerName}`));
        // TODO: add focus return to dock panel
      }
    });

    this.dockPanel = grok.shell.dockManager.dock(ui.div(
      [ui.searchInput('', ''), this.tree]), 'left', null, 'Categories');
    this.dockPanel.container.containerElement.classList.add('tutorials-demo-container');

    // TODO: make div with loading at center of viewer
    // TODO: if loading ended in 0.1s, then no div, if not - then div - DG.debounce, merge etc.
    // TODO: also fix routing things
    // TODO: add starting demo app viewer on just up/down arrows
    // TODO: add switch start on arrow up/down
    // TODO: turn off the toolbox and the menu at top of Datagrok
    // TODO: on click on viewer demo set viewer help url in property panel (func helpUrl)
    // TODO: fix search in demo - search on meta.keywords, name, description
  }
}
