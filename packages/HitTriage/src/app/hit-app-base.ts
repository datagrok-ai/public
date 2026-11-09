/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {AppName, ComputeFunctions, HitDesignCampaign, HitTriageCampaign} from './types';
import {CampaignJsonName, HitTriageComputeFunctionTag} from './consts';
import {discoverComputeFunctions} from '@datagrok-libraries/statistics/src/compute-functions/discovery';
import {checkEditPermissions} from './utils';
import * as api from '../package-api';

/** How many codes a first save tries before it gives up: each one taken by another session meanwhile. */
const MaxNewCodeAttempts = 20;

/** The campaign's lock stayed with another session for the whole wait: what "someone else might be editing" means. */
export class CampaignLockHeldError extends Error {}

export abstract class HitAppBase<T> {
  public dataFrame?: DG.DataFrame;
  public template?: T;
  public baseUrl!: string;
  public computeFunctions: ComputeFunctions;
  protected hasEditPermission = false;
  public get canEdit(): boolean {return this.hasEditPermission;}
  private _appName: AppName;
  public get appName(): AppName {return this._appName;}
  /** Who holds the campaign locks this app takes: the user, and this app's instance — a lock is dropped only by the
   * session that took it, not by another tab of the same user. */
  protected readonly lockHolder =
    `${DG.User.current()?.friendlyName || DG.User.current()?.name || 'unknown'} · ${DG.Utils.randomString(6)}`;
  public static molFileExtReaders: { ext: string; handlerFunc: DG.Func }[] = ['sdf', 'mol', 'smi', 'mol2']
    .map((ext) => {
      const handlerFunc =
      DG.Func.find({meta: {role: DG.FUNC_TYPES.FILE_IMPORTER}}).concat(DG.Func.find({tags: ['file-handler']}))
        .find((f) => f?.options?.ext && typeof f.options.ext === 'string' && f.options.ext.split(',').includes(ext));
      return {ext, handlerFunc: handlerFunc as DG.Func};
    })
    .filter((it) => it.handlerFunc);

  // public layouts: Promise<DG.ViewLayout[]>;
  constructor(public parentCall: DG.FuncCall, appN: AppName) {
    this.baseUrl = new URL(window.location.href).origin + `/apps/HitTriage/${appN.replace(' ', '')}`;
    this._appName = appN;
    this.computeFunctions = discoverComputeFunctions(HitTriageComputeFunctionTag);

    // this.layouts = new Promise<DG.ViewLayout[]>(async (resolve) => {
    //   // we need all layouts, as applicable ones might not be enaugh
    //   const layouts = await grok.dapi.layouts.list();
    //   resolve(layouts);
    // });
  }

  protected async setCanEdit(campaign: HitTriageCampaign | HitDesignCampaign) {
    this.hasEditPermission = await this.canEditCampaign(campaign);
  }

  /** Whether the user may edit `campaign`: always a new one (none, or no author yet). */
  protected async canEditCampaign(campaign?: HitTriageCampaign | HitDesignCampaign): Promise<boolean> {
    return !campaign?.authorUserId || !campaign.permissions ||
      await checkEditPermissions(campaign.authorUserId, campaign.permissions);
  }

  /** The table's molecule columns, never a hidden `~` one (Chem's `~<column>.canonicalSmiles` cache, say). */
  static moleculeColumns(df: DG.DataFrame): DG.Column[] {
    return df.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE).filter((col) => !col.name.startsWith('~'));
  }

  static moleculeColumn(df: DG.DataFrame): DG.Column | null {
    return HitAppBase.moleculeColumns(df)[0] ?? null;
  }

  protected getFilterType(colName: string): DG.FILTER_TYPE {
    const col = this.dataFrame!.col(colName);
    if (col?.semType === DG.SEMTYPE.MOLECULE)
      return DG.FILTER_TYPE.SUBSTRUCTURE;
    if (col?.type === DG.COLUMN_TYPE.BOOL)
      return DG.FILTER_TYPE.BOOL_COLUMNS;
    if (col?.type === DG.COLUMN_TYPE.STRING)
      return DG.FILTER_TYPE.CATEGORICAL;
    return DG.FILTER_TYPE.HISTOGRAM;
  }

  public download(df: DG.DataFrame, name: string, onlyFiltered = false): void {
    const element = document.createElement('a');
    const result = DG.DataFrame.fromColumns(df.columns.toList().filter((c) => !c.name.startsWith('~')))
      .toCsv({filteredRowsOnly: onlyFiltered});
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
    element.setAttribute('download', name + '.csv');
    element.click();
  }

  public async getNewCampaignName(folderName: string, templateKey: string) {
    const templateCampaigns = (await _package.files.list(folderName))
      .map((file) => file.name)
      .filter((name) => name.startsWith(`${templateKey}-`) && !!parseInt(name.substring(templateKey.length + 1)));
    if (templateCampaigns.length === 0)
      return templateKey + '-1';
    const postFixes = templateCampaigns
      .map((c) => c.substring(templateKey.length + 1))
      .filter(Boolean).map((c) => parseInt(c, 10))
      .sort((a, b) => a - b);// sort in ascending order
    return templateKey + '-' + ((postFixes[postFixes.length - 1] + 1).toString());
  }

  /** The code a campaign's first save writes under, its lock held: when another session saved a campaign under
   * `campaignId` meanwhile, the next free code of the template's `key` is locked instead — the check is made under
   * the lock, which every writer of a campaign holds, so two sessions never share a code. Throws, holding no lock,
   * when no code could be locked; the caller reports it. */
  protected async lockNewCampaignCode(campaignId: string, key: string): Promise<string> {
    const numberOf = (code: string) => parseInt(code.substring(key.length + 1), 10) || 0;
    for (let attempt = 0; attempt < MaxNewCodeAttempts; attempt++) {
      const lock = await this.aquireAndWaitCampaignLock(campaignId);
      if (!lock.aquired)
        throw new CampaignLockHeldError(`Could not lock campaign ${campaignId} within 30 seconds`);
      let free = false;
      try {
        free = !await _package.files.exists(`${this.appName}/campaigns/${campaignId}/${CampaignJsonName}`);
      } finally {
        if (!free)
          await this.dropCampaignLock(campaignId);
      }
      if (free)
        return campaignId;
      const listed = await this.getNewCampaignName(`${this.appName}/campaigns`, key);
      campaignId = `${key}-${Math.max(numberOf(campaignId) + 1, numberOf(listed))}`;
    }
    throw new Error(`No free campaign code of ${key} after ${MaxNewCodeAttempts} attempts`);
  }

  public setBaseUrl() {
    const title = document.title;
    if (history.replaceState) {
      const obj = {
        Title: title,
        Url: this.baseUrl,
      };
      history.replaceState(obj, obj.Title, obj.Url);
    }
  }

  protected async aquireAndWaitCampaignLock(campaignId: string, maxWaitMs = 30_000): Promise<{aquired: boolean, expires_at?: null | string}> {
    const startTime = Date.now();
    let errorCount = 0;
    let lastError: any = null;
    while (Date.now() - startTime < maxWaitMs && errorCount < 5) {
      try {
        const result = await api.queries.acquireCampaignLock(this.appName, campaignId, this.lockHolder);
        if (result && result.rowCount > 0) {
          const res = {aquired: true, expires_at: result?.col('expires_at')?.getString(0)};
          return res;
        }
        lastError = null;
        await DG.delay(200);
      } catch (err) {
        // console.error('Error acquiring campaign lock:', err);
        lastError = err;
        errorCount++;
      }
    }
    if (lastError)
      throw lastError;
    return {aquired: false, expires_at: null};
  }

  /** Lets this session's lock go without logging a save (nothing was written under it); a failure is only logged. */
  protected async dropCampaignLock(campaignId: string): Promise<void> {
    await api.queries.dropCampaignLock(this.appName, campaignId, this.lockHolder)
      .catch((e) => _package.logger.error(e));
  }

  protected async releaseCampaignLock(campaignId: string): Promise<string | null> {
    return (await api.queries.releaseCampaignLock(this.appName, campaignId))?.col('updated_at')?.getString(0) || null;
  }
}
