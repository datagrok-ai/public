import {HitTriageBaseView} from './base-view';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {_package} from '../../package';
import {ICampaign, ICampaignStatus} from '../types';
import {CampaignJsonName, CampaignTableName} from '../consts';
import {saveCampaignDialog} from '../dialogs/save-campaign-dialog';
import {toFormatedDateString} from '../utils';

export class SubmitView extends HitTriageBaseView {
  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Submit';
  }

  render(): void {
    ui.empty(this.root);
    const submitDiv = ui.divH([]);
    if (this.app.template?.submit && this.app.template.submit.fName)
      submitDiv.appendChild(ui.bigButton('SUBMIT', () => this.submit()));

    submitDiv.appendChild(ui.bigButton('Save Campaign', () => this.saveCampaign()));

    const content = ui.divV([
      ui.h1('Summary'),
      ui.div([ui.tableFromMap(this.app.getSummary())]),
      submitDiv,
    ]);
    this.root.appendChild(content);
  }

  onActivated(): void {
    this.render();
  }

  async submit(): Promise<any> {
    const submitParams= this.app.template?.submit;
    if (!submitParams)
      return;
    const submitFn = DG.Func.find({name: submitParams.fName, package: submitParams.package})[0];
    if (!submitFn) {
      grok.shell.error(`Function ${submitParams.fName} not found.`);
      return;
    }
    const filteredDf = DG.DataFrame.fromCsv(this.app.dataFrame!.toCsv({filteredRowsOnly: true}));
    await submitFn.apply({df: filteredDf, molecules: this.app.molColName});
    this.app.campaign && (this.app.campaign.status = 'Submitted');
    this.saveCampaign('Submitted');
    grok.shell.info('Submitted successfully.');
  }

  async saveCampaign(status?: ICampaignStatus): Promise<any> {
    const campaignId = this.app.campaignId!;
    const filters = this.app.filterSettings!;
    const templateName = this.app.template!.name;
    const enrichedDf = this.app.dataFrame!;
    const campaignPrefix = `System:AppData/HitTriage/campaigns/${campaignId}/`;
    const campaignName = campaignId ?? await saveCampaignDialog(campaignId);
    const campaign: ICampaign = {
      name: campaignName,
      templateName,
      filters: filters.slice(1), // TODO: because the first filter is chem substructure which does not work
      ingest: {
        type: 'File',
        query: `${campaignPrefix}${CampaignTableName}`,
        molColName: this.app.molColName!,
      },
      status: status ?? this.app.campaign?.status ?? 'In Progress',
      createDate: this.app.campaign?.createDate ?? toFormatedDateString(new Date()),
      campaignFields: this.app.campaign?.campaignFields ?? this.app.campaignProps,
    };
    await _package.files.writeAsText(`campaigns/${campaignId}/${CampaignJsonName}`, JSON.stringify(campaign));
    await _package.files.writeAsText(`campaigns/${campaignId}/${CampaignTableName}`, enrichedDf.toCsv());
    grok.shell.info('Campaign saved successfully.');
  }
}
