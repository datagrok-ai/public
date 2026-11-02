/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {AppHeaderOptions} from '@datagrok-libraries/u2';
import {HitTriageApp} from '../hit-triage-app';
import {_package} from '../../package';
import {HitTriageCampaign, HitTriageTemplate} from '../types';
import {CampaignIdKey, CampaignJsonName, HitTriageDataSourceTag, i18n} from '../consts';
import '../../../css/hit-triage.css';
import {addBreadCrumbsToRibbons, checkViewPermissions, loadTemplate, modifyUrl, popRibbonPannels} from '../utils';
import {INewCampaignResult, newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import {createTemplateAccordeon} from '../accordeons/new-template-accordeon';
import {CampaignUnavailableError, HitLandingView} from '../base-view';
import {HitHomeOptions} from '../shared/hit-home';
import type {TemplateFormContext} from '../shared/template-picker';
import {hitInfoAiDescription} from '../ai/view-functions';

export class InfoView extends HitLandingView<HitTriageTemplate, HitTriageApp> {
  private dataSourceFunctionsMap: {[key: string]: DG.Func | DG.DataQuery} = {};
  constructor(app: HitTriageApp) {
    super(app, 'README_HT.md');
    this.name = 'Hit Triage';
    // TODO: Remove before api release
    (this as any).aiDescription = hitInfoAiDescription(app.appName);
    this.onActivated();
    void this.init().then(() => this.openLinked());
  }

  protected header(): AppHeaderOptions {
    return {
      icon: _package.webRoot + '/images/icons/hit-triage-icon.png',
      learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md',
      items: [
        'Configure your own workflow using the template editor.',
        'Calculate different molecular properties.',
        'Filter molecules using different criteria.',
        'Submit processed dataframe to the function of your choice.',
        'Save campaigns and continue any time from where you left off.',
      ],
    };
  }

  protected getTemplatesFolder(): string {
    return `${this.app.appName}/templates`;
  }

  protected getDataSourceTag(): string {
    return HitTriageDataSourceTag;
  }

  protected getDataSourceFunctionsMap(): {[key: string]: DG.Func | DG.DataQuery} {
    return this.dataSourceFunctionsMap;
  }

  protected homeOptions(preset?: HitTriageTemplate): HitHomeOptions<HitTriageTemplate> {
    return {
      header: this.header(),
      catalog: this.catalogOptions('HT', (c) => c.friendlyName ?? c.name),
      picker: {app: this.app.appName, preset, edit: (template) => void this.createNewTemplate(template),
        form: (template, context) => this.getNewCampaignForm(template, context)},
    };
  }

  protected async checkCampaign(campaignId: string) {
    // check if such campaign exists
    if (!await _package.files.exists(`${this.app.appName}/campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: HitTriageCampaign =
      JSON.parse(await _package.files.readAsText(`${this.app.appName}/campaigns/${campaignId}/${CampaignJsonName}`));
    // in case if the link was opened and user has no permissions to view the campaign
    if (campaign.authorUserId && campaign.permissions &&
        !await checkViewPermissions(campaign.authorUserId, campaign.permissions)
    ) {
      this.app.setBaseUrl();
      return;
    }
    // what the open reads is checked before the app takes the campaign: a failed open leaves the open one in place
    const templatePath = `${this.getTemplatesFolder()}/${campaign.templateName}.json`;
    if (!await _package.files.exists(templatePath))
      throw new CampaignUnavailableError('its template is missing');
    if (campaign.ingest.type === 'File' && !await grok.dapi.files.exists(campaign.ingest.query))
      throw new CampaignUnavailableError(`its table ${campaign.ingest.query} is missing`);
    const template: HitTriageTemplate = await loadTemplate<HitTriageTemplate>(templatePath);
    this.app.campaign = campaign;
    // modify the template with path to the campaign's precalculated table
    await this.app.setTemplate(template, campaign.filters, campaignId, campaign.ingest);

    modifyUrl(CampaignIdKey, campaignId);

    return campaign;
  }

  protected getNewCampaignForm(template: HitTriageTemplate, {owner, column}: TemplateFormContext): Promise<HTMLElement> {
    return newCampaignAccordeon(template, this.getDataSourceFunctionsMap(), this.getDataSourceTag(),
      {owner, column, start: (camp) => this.startCampaign(template, camp)});
  }

  /** A new campaign: nothing of the campaign opened before it carries over. */
  protected async startCampaign(template: HitTriageTemplate, camp: INewCampaignResult) {
    this.app.campaign = undefined;
    this.app.dataFrame = camp.df;
    this.app._fileInputType = camp.type;
    await this.app.setTemplate(template, undefined, undefined, undefined, camp.friendlyName);
    this.app.campaignProps = camp.campaignProps;
    const campaign: HitTriageCampaign | undefined = await this.app.saveCampaign(undefined, false);
    if (template.layoutViewState && campaign)
      campaign.layout = template.layoutViewState;
  }

  protected async createNewTemplate(preset?: HitTriageTemplate) {
    ui.setUpdateIndicator(this.root, true);
    try {
      const newTemplateAccordeon = await createTemplateAccordeon(this.app, this.getDataSourceFunctionsMap(), preset);
      this.showTemplateEditor(newTemplateAccordeon);
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  protected showTemplateEditor(newTemplateAccordeon: {root: HTMLElement, template: Promise<HitTriageTemplate>,
    cancelPromise: Promise<void>}) {
    const newView = new DG.ViewBase();
    const curView = grok.shell.v;
    newView.name = 'New Template';
    newView.root.appendChild(newTemplateAccordeon.root);
    newView.parentCall = this.app.parentCall;
    grok.shell.addView(newView);
    newView.path = new URL(this.app.baseUrl).pathname + '/new-template';
    newView.parentView = curView;
    const {sub} = addBreadCrumbsToRibbons(grok.shell.v, this.app.appName, i18n.createNewTemplate, () => {
      grok.shell.v = curView;
      newView.close();
    });

    newTemplateAccordeon.template.then(async (t) => {
      await this.init(t);
      newView.close();
      popRibbonPannels(newView);
      grok.shell.v = curView;
      sub.unsubscribe();
    });
    newTemplateAccordeon.cancelPromise.then(() => {
      sub.unsubscribe();
      newView.close();
      grok.shell.v = curView;
    });
  }
};
