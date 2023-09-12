import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {TileCategoriesColName} from '../consts';
import {HitBaseView} from '../base-view';
import {HitDesignTemplate} from '../types';
import {checkRibbonsHaveSubmit} from '../utils';

export class HitDesignTilesView extends HitBaseView<HitDesignTemplate, HitDesignApp> {
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Progress Tracker';
  }

  async render(): Promise<void> {
    return new Promise<void>((resolve) =>{
      ui.empty(this.root);

      //const tv = DG.TableView.create(this.app.dataFrame!, false);
      const v = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, this.app.dataFrame!,
        {lanesColumnName: TileCategoriesColName, autoGenerate: true});
      //this.root.appendChild(tv.root);
      this.root.appendChild(v.root);
      v.root.style.height = '100%';
      v.root.style.width = '100%';
      setTimeout(() => {
        v.view?._onAdded();
        const ribbons = grok.shell.v?.getRibbonPanels();

        if (ribbons) {
          const hasSubmit = checkRibbonsHaveSubmit(ribbons);
          if (hasSubmit) {
            resolve();
            return;
          }
          const submitButton = ui.bigButton('Submit', () => {
            const dialogContent = this.app._submitView?.render();
            if (dialogContent) {
              const dlg = ui.dialog('Submit');
              dlg.add(dialogContent);
              dlg.addButton('Save', ()=>{this.app.saveCampaign(); dlg.close();});
              dlg.addButton('Submit', ()=>{this.app._submitView?.submit(); dlg.close();});
              dlg.show();
            }
          });
          submitButton.classList.add('hit-design-submit-button');
          ribbons.push([submitButton]);
          grok.shell.v.setRibbonPanels(ribbons);
        }
        resolve();
      }, 5);
    });

    //const tv = DG.TableView.create(this.app.dataFrame!, false);
  }

  onActivated(): void {
    this.render();
  }
}
