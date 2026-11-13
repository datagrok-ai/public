import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Scope} from '@datagrok-libraries/u2';
import type {Control} from '@datagrok-libraries/u2';
import {host} from '@datagrok-libraries/u2/dg';
import {AppName, HitDesignTemplate, HitTriageTemplate, PeptiHitTemplate} from './types';
import type {HitDesignApp} from './hit-design-app';
import type {HitTriageApp} from './hit-triage-app';
import {CampaignIdKey} from './consts';
import {_package} from '../package';
import * as api from '../package-api';
import {HitHome, HitHomeOptions, hitHome} from './shared/hit-home';
import {AnyCampaign, CampaignCatalogOptions, CampaignFamily} from './shared/campaign-catalog';

export class HitBaseView<Ttemplate extends HitDesignTemplate | HitTriageTemplate | PeptiHitTemplate,
    Tapp extends HitDesignApp | HitTriageApp> extends DG.ViewBase {
  app: Tapp;
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
}

/** A campaign the landing cannot open for a reason the user is told (its template or its table is gone): an expected
 * state of the data, not a failure of the app, so it is not logged as an error. */
export class CampaignUnavailableError extends Error {
  constructor(reason: string) {
    super(reason);
    this.name = 'CampaignUnavailableError';
  }
}

/** An app's landing ("Hit Design", "Hit Triage", …): its HitHome, and the campaign a link names opened over it. */
export abstract class HitLandingView<Ttemplate extends HitDesignTemplate | HitTriageTemplate | PeptiHitTemplate,
    Tapp extends HitDesignApp | HitTriageApp> extends HitBaseView<Ttemplate, Tapp> {
  protected home?: Control & HitHome;
  /** Deleted in this page: a share listing may show their folders a while longer. */
  protected deletedCampaigns: string[] = [];

  /** [readme]: the app's help page in the package (README_HD.md), kept as the landing's `helpUrl` — the app gives it
   * to the view the shell shows, whose help the platform shows (a help shown by hand loses to the panel's own
   * page for the startup placeholder). */
  constructor(app: Tapp, readme: string) {
    super(app);
    this.helpUrl = _package.webRoot + readme;
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(this.helpUrl);
  }

  protected abstract homeOptions(preset?: Ttemplate): HitHomeOptions<Ttemplate>;

  /** Opens the campaign with this code; answers nothing when there is none the user may view. Throws when it cannot
   * be opened. */
  protected abstract checkCampaign(campaignId: string): Promise<unknown>;

  /** Builds the landing again: the table reads the page's campaign list, the picker starts on `preset`. */
  async init(preset?: Ttemplate): Promise<void> {
    this.home?.dispose();
    // a tree of its own, not the caller's: a campaign view's crumb asks for the landing from a scope that is closing
    const owner = Scope.root('hit-home');
    this.home = owner.build(() => hitHome<Ttemplate>(this.homeOptions(preset)));
    this.home.own(() => owner.dispose());
    this.root.replaceChildren(host(this.home));
  }

  /** The campaign's files, then what the database keeps of it (its last save, its lock, its V-iD links);
   * rows the database fails to drop are only logged — the campaign is gone. */
  async deleteCampaign(appName: AppName, campaignId: string): Promise<void> {
    await _package.files.delete(`${appName}/campaigns/${campaignId}`);
    delete (_package.campaignsCache[appName] as {[name: string]: unknown} | undefined)?.[campaignId];
    this.deletedCampaigns.push(campaignId);
    await api.queries.deleteCampaignRecords(appName, campaignId).catch((e) => _package.logger.error(e));
  }

  /** The campaign deletion the UI and the assistant share: the files go, and the table drops the row. */
  async deleteCampaignAndRefresh(appName: AppName, campaignId: string): Promise<void> {
    await this.deleteCampaign(appName, campaignId);
    await this.home?.catalog.reload();
  }

  /** A campaign's name link, the assistant's open: a campaign that cannot open says why. */
  async setCampaign(campaignId: string): Promise<void> {
    await this.open(() => this.checkCampaign(campaignId), campaignId);
  }

  /** The campaign of the link the app was opened with (`?campaignId=`), opened over the landing; one that cannot open
   * says why and leaves the landing in place. The router passes the link's id to the app's call: the address may have
   * changed since. */
  protected async openLinked(): Promise<void> {
    const linked = this.app.parentCall.getParamValue(CampaignIdKey) ??
      new URLSearchParams(location.search).get(CampaignIdKey);
    if (linked)
      await this.open(() => this.checkCampaign(linked), linked);
  }

  /** This app's campaigns table: whose settings it shares, what its balloons call a campaign. */
  protected catalogOptions(family: CampaignFamily, label: (c: AnyCampaign) => string,
    onLoaded?: (campaigns: readonly AnyCampaign[]) => void): CampaignCatalogOptions {
    return {app: this.app.appName, family, deleted: this.deletedCampaigns, label, onLoaded,
      open: (code) => this.setCampaign(code), remove: (code) => this.deleteCampaign(this.app.appName, code)};
  }

  private async open(check: () => Promise<unknown>, campaignId: string): Promise<void> {
    const previous = this.app.campaign;
    const table = this.app.dataFrame;
    try {
      await check();
    }
    catch (e) {
      // the campaign still open saves as itself, not as the one that failed to open; one the app took whole stays
      if (this.app.dataFrame === table)
        this.app.campaign = previous;
      if (!(e instanceof CampaignUnavailableError))
        _package.logger.error(e);
      grok.shell.error(`Campaign ${campaignId} could not be opened: ${e instanceof Error ? e.message : e}`);
      this.app.setBaseUrl();
    }
  }
}
