import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {TileCategoriesColName} from '../consts';
import {HitBaseView} from '../base-view';
import {HitDesignTemplate} from '../types';
import {checkRibbonsHaveSubmit} from '../utils';
import './utils.css';

export class HitDesignTilesView extends HitBaseView<HitDesignTemplate, HitDesignApp> {
  private v: DG.Viewer | null = null;
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Progress Tracker';
  }

  async render(): Promise<void> {
    return new Promise<void>((resolve) =>{
      ui.empty(this.root);

      // const filteredDf = DG.DataFrame.fromColumns(
      //   this.app.dataFrame!.columns.toList().filter((col) => !col.name.startsWith('~')));
      // const _hiidenTv = DG.TableView.create(filteredDf, true);
      //const tv = DG.TableView.create(this.app.dataFrame!, false);
       this.app.dataFrame!.columns.names().forEach((name) => {
         if (name.startsWith('~'))
          this.app.dataFrame!.columns.remove(name);
       });
       this.v = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, this.app.dataFrame!,
         {lanesColumnName: TileCategoriesColName, lanes: this.app.template?.stages ?? []});
       //this.root.appendChild(tv.root);
       this.root.appendChild(this.v.root);
       this.v.root.style.height = '100%';
       this.v.root.style.width = '100%';
       setTimeout(() => {
         this.v?.view?._onAdded();
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

  public destroy(): void {
    this.v?.detach();
  }
}
