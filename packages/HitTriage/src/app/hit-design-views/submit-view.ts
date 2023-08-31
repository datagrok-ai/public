import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {HitDesignTemplate} from '../types';
import {HitBaseView} from '../base-view';

export class HitDesignSubmitView extends HitBaseView<HitDesignTemplate, HitDesignApp> {
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Submit';
  }

  render(): HTMLDivElement {
    ui.empty(this.root);
    const submitDiv = ui.divH([]);
    if (this.app.template?.submit && this.app.template.submit.fName)
      submitDiv.appendChild(ui.bigButton('SUBMIT', () => this.submit()));

    submitDiv.appendChild(ui.bigButton('Save Campaign', () => this.app.saveCampaign()));
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
    this.app.saveCampaign('Submitted');
    grok.shell.info('Submitted successfully.');
  }
}
