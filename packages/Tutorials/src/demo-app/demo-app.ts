// import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// import {filter} from "rxjs/operators";

export class DemoView extends DG.ViewBase {
  tree: DG.TreeViewGroup = ui.tree();

  constructor() {
    super();

    this._initToolbox();
    this._initContent();
  }

  static findDemoFunc(demoPath: string) {
    return DG.Func.find({meta: {'demoPath': demoPath}})[0];
  }

  _initContent() {
    this.root.appendChild(ui.divText('Select a demo from the toolbox on the right'));
  }

  _initToolbox() {
    for (const f of DG.Func.find({meta: {'demoPath': null}})) {
      const pathOption = <string>f.options[DG.FUNC_OPTIONS.DEMO_PATH];
      const path = pathOption.split('|').map((s) => s.trim());
      const folder = this.tree.getOrCreateGroup(path.slice(0, path.length - 1).join(' | '));
      const item = folder.item(path[path.length - 1]);
      item.root.onmousedown = (_) => {
        // grok.shell.closeAll();
        f.apply().then((_) => { });
      };
    }

    this.toolbox = this.tree.root;
  }
}
