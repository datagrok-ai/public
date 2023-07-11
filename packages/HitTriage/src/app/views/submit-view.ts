import {HitTriageBaseView} from './base-view';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
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
    const content = ui.divV([
      ui.h1('Summary'),
      ui.div([ui.tableFromMap(this.app.getSummary())]),
      ui.divH([ui.bigButton('SUBMIT', () => this.submit()), ui.bigButton('Save Campaign', () => this.saveCampaign())]),
    ]);
    this.root.appendChild(content);
  }

  onActivated(): void {
    this.render();
  }

  async submit(): Promise<any> {
    //TODO
    // const d = new Date();
    // const time =
    //   `${d.getFullYear()}-${d.getMonth()}-${d.getDate()}_${d.getHours()}-${d.getMinutes()}-${d.getSeconds()}`;
    // const folder = `${time}_${grok.shell.user.login}`;
    // await _package.files.writeAsText(`${folder}/session.json`, JSON.stringify(this.app.template));

    // await _package.files.writeAsText(`${folder}/molecules.csv`, this.app.template.enrichedTable!.toCsv());

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
