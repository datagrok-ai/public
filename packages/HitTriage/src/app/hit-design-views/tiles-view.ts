import {HitDesignBaseView} from './base-view';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {TileCategoriesColName} from '../consts';

export class HitDesignTilesView extends HitDesignBaseView {
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Progress Tracker';
  }

  render(): void {
    ui.empty(this.root);

    //const tv = DG.TableView.create(this.app.dataFrame!, false);
    const v = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, this.app.dataFrame!,
      {lanesColumnName: TileCategoriesColName, autoGenerate: true});
    //this.root.appendChild(tv.root);
    this.root.appendChild(v.root);
    v.root.style.height = '100%';
    setTimeout(() => {
      v.view?._onAdded();
      const ribbons = grok.shell.v?.getRibbonPanels();
      if (ribbons) {
        const submitButton = ui.div(ui.bigButton('Submit', () => {
          const dialogContent = this.app._submitView?.render();
          if (dialogContent)
            ui.dialog('Submit').add(dialogContent).show();
        }));
        ribbons.push([submitButton]);
        grok.shell.v.setRibbonPanels(ribbons);
      }
    }, 5);

    //const tv = DG.TableView.create(this.app.dataFrame!, false);
  }

  onActivated(): void {
    this.render();
  }
}
