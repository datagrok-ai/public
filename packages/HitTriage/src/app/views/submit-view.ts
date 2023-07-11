import {HitTriageBaseView} from './base-view';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {_package} from '../../package';
import {ICampaign} from '../types';
import {CampaignJsonName, CampaignTableName} from '../consts';

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
    const submitFname = this.app.template!.submit!.fName;
    const submitFn = await grok.functions.find(submitFname);
    if (!submitFn) {
      grok.shell.error(`Function ${submitFname} not found.`);
      return;
    }
    await submitFn.apply({df: this.app.dataFrame, molecules: this.app.template?.ingest.molColName});
    grok.shell.info('Submitted successfully.');
  }

  async saveCampaign(): Promise<any> {
    const campaignId = this.app.campaignId!;
    const filters = this.app.filterSettings!;
    const templateName = this.app.template!.name;
    const enrichedDf = this.app.dataFrame!;
    const campaign: ICampaign = {
      templateName,
      filters: filters.slice(1), // TODO: because the first filter is chem substructure which does not work
    };
    await _package.files.writeAsText(`campaigns/${campaignId}/${CampaignJsonName}`, JSON.stringify(campaign));
    await _package.files.writeAsText(`campaigns/${campaignId}/${CampaignTableName}`, enrichedDf.toCsv());
    grok.shell.info('Campaign saved successfully.');
  }
}
