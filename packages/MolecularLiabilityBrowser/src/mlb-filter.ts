import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class MlbFilterViewer extends DG.JsViewer {
  constructor() {
    super();
  }

  get mlbDf(): DG.DataFrame {
    return this.dataFrame;
  }

  public async init(): Promise<void> {

  }

  override async onTableAttached() {
    await this.init();
  }

  async setData(mlbDf: DG.DataFrame): Promise<void> {
    await this.destroyView();
    this.dataFrame = mlbDf;
    await this.buildView();
  }

  async destroyView() {

  }

  async buildView() {

  }
}
