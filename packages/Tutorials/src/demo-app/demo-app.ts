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
        grok.shell.closeAll();
        f.apply().then((_) => { });
      };
    }
    this.dockPanel = grok.shell.dockManager.dock(ui.div(
      [ui.searchInput('', ''), this.tree]), 'left', null, 'Categories');
    this.dockPanel.container.containerElement.style.maxWidth = '250px';
  }
}
