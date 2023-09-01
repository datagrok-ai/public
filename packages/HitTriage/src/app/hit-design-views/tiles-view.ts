import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {TileCategoriesColName} from '../consts';
import {HitBaseView} from '../base-view';
import {HitDesignTemplate} from '../types';

export class HitDesignTilesView extends HitBaseView<HitDesignTemplate, HitDesignApp> {
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
        const hasSubmit = ribbons.reduce((prev, cur) => {
          return prev || cur.reduce((p, c) =>
            p || c.classList.contains('hit-design-submit-button') ||
              Array.from(c.children).some((child) =>child.classList.contains('hit-design-submit-button')), false);
        }, false);
        if (hasSubmit)
          return;
        const submitButton = ui.div(ui.bigButton('Submit', () => {
          const dialogContent = this.app._submitView?.render();
          if (dialogContent)
            ui.dialog('Submit').add(dialogContent).show().getButton('CANCEL').textContent = 'Ok';
        }), {classes: 'hit-design-submit-button'});
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
