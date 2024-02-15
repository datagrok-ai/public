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
    const funccall = await this.func.prepare().call();

    return [funccall];
  }
}