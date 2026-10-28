/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import type {AppHeaderOptions} from '@datagrok-libraries/u2';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {CampaignJsonName, CampaignTableName, HitDesignCampaignIdKey} from '../consts';
import {HitDesignCampaign, HitDesignTemplate} from '../types';
import {checkViewPermissions, loadTemplate, modifyUrl} from '../utils';
import {CampaignUnavailableError, HitLandingView} from '../base-view';
import {HitHomeOptions} from '../shared/hit-home';
import {newCampaignForm} from '../shared/new-campaign-form';
import {openTemplateEditor} from '../shared/template-editor';
import {hitInfoAiDescription} from '../ai/view-functions';

export class HitDesignInfoView
  <T extends HitDesignTemplate = HitDesignTemplate, K extends HitDesignApp = HitDesignApp>
  extends HitLandingView<T, K> {
  constructor(app: K) {
    super(app, 'README_HD.md');
    this.name = 'Hit Design';
    // TODO: Remove before api release
    (this as any).aiDescription = hitInfoAiDescription(app.appName);
    this.onActivated();
    void this.init().then(() => this.openLinked());
  }

  protected header(): AppHeaderOptions {
    return {
      icon: _package.webRoot + '/images/icons/hit-design-icon.png',
      learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HD.md',
      items: [
        'Configure your own workflow using the template editor',
        'Sketch molecules in the molecular spreadsheet',
        'Annotate and share ideas with the team',
        'Calculate different molecular properties',
        'Save campaigns and continue from where you left off',
        'Submit final selection to the function of your choice',
      ],
    };
  }

  protected homeOptions(preset?: T): HitHomeOptions<T> {
    return {
      header: this.header(),
      catalog: this.catalogOptions('HD', (c) => c.name, (campaigns) =>
        this.app.existingStatuses = Array.from(new Set(campaigns.map((c) => c.status).filter((s) => !!s)))),
      picker: {app: this.app.appName, preset,
        edit: (template) => void openTemplateEditor({app: this.app, preset: template, onCreated: (t) => this.init(t)}),
        form: (template, {column}) => newCampaignForm({template, column,
          nextCode: () => this.app.getNewCampaignName(`${this.app.appName}/campaigns`, template.key),
          start: ({fields, name}) => this.app.startNewCampaign(template, fields, name)})},
    };
  }

  protected async checkCampaign(campId?: string) {
    const url = location.search;
    const urlParams = new URLSearchParams(url);
    if (!urlParams.has(HitDesignCampaignIdKey) && !campId)
      return;
    const campaignId = campId ?? urlParams.get(HitDesignCampaignIdKey);
    // check if such campaign exists
    if (!await _package.files.exists(`${this.app.appName}/campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: HitDesignCampaign =
      JSON.parse(await _package.files.readAsText(`${this.app.appName}/campaigns/${campaignId}/${CampaignJsonName}`));
    // in case if the link was opened and user has no permissions to view the campaign
    if (campaign.authorUserId && campaign.permissions &&
      !await checkViewPermissions(campaign.authorUserId, campaign.permissions)
    ) {
      this.app.setBaseUrl();
      return;
    }
    // what the open reads is checked before the app takes the campaign: a failed open leaves the open one in place
    const templatePath = `${this.app.appName}/templates/${campaign.templateName}.json`;
    if (!campaign.template && !await _package.files.exists(templatePath))
      throw new CampaignUnavailableError('its template is missing');
    const tablePath = campaign.savePath ??
      `System.AppData/HitTriage/${this.app.appName}/campaigns/${campaignId}/${CampaignTableName}`;
    if (!await grok.dapi.files.exists(tablePath))
      throw new CampaignUnavailableError(`its table ${tablePath} is missing`);
    const template: T = (campaign.template as T) ?? await loadTemplate<T>(templatePath);
    await this.app.setTemplate(template, campaignId!, campaign);

    if (campId)
      modifyUrl(HitDesignCampaignIdKey, campId);

    campaign.template = template;
    return campaign;
  }
}
