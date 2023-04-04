import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

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
    func.apply().then((_) => { });
  }

  _initContent() {
    this.root.appendChild(ui.divText('Select a demo from the toolbox on the right'));
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

    this.tree.root.onkeyup = (event) => {
      const currentElement = document.getElementsByClassName('d4-tree-view-node-selected')[0];
      if (event.code === 'Enter' && currentElement.classList.contains('d4-tree-view-item')) {
        const categoryName = currentElement.parentElement?.parentElement
          ?.getElementsByClassName('d4-tree-view-group-label')[0].innerHTML;
        const viewerName = currentElement.getElementsByClassName('d4-tree-view-item-label')[0].innerHTML;

        this.startDemoFunc(DemoView.findDemoFunc(`${categoryName} | ${viewerName}`));
        // TODO: add focus return to dock panel
      }
    };

    this.dockPanel = grok.shell.dockManager.dock(ui.div(
      [ui.searchInput('', ''), this.tree]), 'left', null, 'Categories');
    this.dockPanel.container.containerElement.style.maxWidth = '250px';

    // TODO: make div with loading at center of viewer
    // TODO: if loading ended in 0.1s, then no div, if not - then div - DG.debounce, merge etc.
    // TODO: also fix routing things
  }
}
