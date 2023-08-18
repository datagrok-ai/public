import {HitDesignBaseView} from './base-view';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {HitDesignCampaign, HitTriageCampaignStatus} from '../types';
import {toFormatedDateString} from '../utils';
import {CampaignJsonName, CampaignTableName} from '../consts';

export class HitDesignSubmitView extends HitDesignBaseView {
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Submit';
  }

  render(): HTMLDivElement {
    ui.empty(this.root);
    const submitDiv = ui.divH([]);
    if (this.app.template?.submit && this.app.template.submit.fName)
      submitDiv.appendChild(ui.bigButton('SUBMIT', () => this.submit()));

    submitDiv.appendChild(ui.bigButton('Save Campaign', () => this.saveCampaign()));
    submitDiv.style.gap = '10px';
    const content = ui.divV([
      ui.h1('Summary'),
      ui.div([ui.tableFromMap(this.app.getSummary())]),
      submitDiv,
    ]);
    return content;
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

  async saveCampaign(status?: HitTriageCampaignStatus): Promise<any> {
    const campaignId = this.app.campaignId!;
    const templateName = this.app.template!.name;
    const enrichedDf = this.app.dataFrame!;
    const campaignName = campaignId;
    const campaign: HitDesignCampaign = {
      name: campaignName,
      templateName,
      status: status ?? this.app.campaign?.status ?? 'In Progress',
      createDate: this.app.campaign?.createDate ?? toFormatedDateString(new Date()),
      campaignFields: this.app.campaign?.campaignFields ?? this.app.campaignProps,
    };
    await _package.files.writeAsText(`Hit Design/campaigns/${campaignId}/${CampaignJsonName}`,
      JSON.stringify(campaign));

    const csvDf = DG.DataFrame.fromColumns(
      enrichedDf.columns.toList().filter((col) => !col.name.startsWith('~')),
    ).toCsv();
    await _package.files.writeAsText(`Hit Design/campaigns/${campaignId}/${CampaignTableName}`, csvDf);
    grok.shell.info('Campaign saved successfully.');
  }
}
