import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

export class CustomDataUploader extends DG.Widget {
  constructor(private func: DG.Func) {
    const root = ui.panel([ui.divV([
      ui.label('My custom data uploader')
    ])]);
    super(root);
  }

  async getFunccalls(): Promise<DG.FuncCall[]> {
    // Replace this code by DB connection / API call / etc.
    const funccall = await this.func.prepare({
      'ambTemp': 22,
      'initTemp': 100,
      'desiredTemp': 30,
      'area': 0.06,
      'heatCap': 4200,
      'heatTransferCoeff': 8.3,
      'simTime': 21600,
      }).call();

    return [funccall];
  }
}