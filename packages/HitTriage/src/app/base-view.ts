import * as DG from 'datagrok-api/dg';
import {AppName, HitDesignTemplate, HitTriageTemplate, PeptiHitTemplate} from './types';
import type {HitDesignApp} from './hit-design-app';
import type {HitTriageApp} from './hit-triage-app';
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

  /** Info views implement this to rebuild their content. */
  init?(presetTemplate?: any): Promise<void>;

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
    delete (_package.campaignsCache[appName] as {[name: string]: unknown} | undefined)?.[campaignId];
  }

  /** The full campaign-deletion flow the UI and the AI share: delete the files,
   * hide the campaign from the cached list, and rebuild the view. */
  async deleteCampaignAndRefresh(appName: AppName, campaignId: string): Promise<void> {
    await this.deleteCampaign(appName, campaignId);
    this.deletedCampaigns.push(campaignId);
    await this.init?.();
  }
}
