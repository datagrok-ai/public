import * as DG from 'datagrok-api/dg';
import {AppName, HitDesignTemplate, HitTriageTemplate, PeptiHitTemplate} from './types';
import {HitDesignApp} from './hit-design-app';
import {HitTriageApp} from './hit-triage-app';
import {_package} from '../package';

export class HitBaseView<Ttemplate extends HitDesignTemplate | HitTriageTemplate | PeptiHitTemplate,
    Tapp extends HitDesignApp | HitTriageApp> extends DG.ViewBase {
  app: Tapp;
  protected deletedCampaigns: string[] = [];
  constructor(app: Tapp) {
    super();

    this.app = app;
  }

  get template(): Ttemplate {return this.app.template! as Ttemplate;}

  /** Override to initialize the view based on the session. */
  onActivated(): void {
  }

  /** Gets called when any of the previous views change. */
  onReset() {
  }

  async process(): Promise<any> {
  }

  async deleteCampaign(appName: AppName, campaignId: string): Promise<void> {
    await _package.files.delete(`${appName}/campaigns/${campaignId}`);
  }
}
