import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {FileSource} from 'datagrok-api/dg';

export class OligoSdFileApp {
  df!: DG.DataFrame;

  constructor() {}

  async init(srcDf?: DG.DataFrame): Promise<void> {

    let dataDf: DG.DataFrame;
    if (srcDf) {
      dataDf = srcDf;
    } else {
      const dfFn: string = 'System:AppData/SequenceTranslator/test input_Nov28_Duplex_dimer.xlsx';
      dataDf = await this.loadData(dfFn);
    }

    this.setData(dataDf);
  }

  async loadData(dfFn: string): Promise<DG.DataFrame> {
    //
    const dataDf: DG.DataFrame = await grok.data.files.openTable(dfFn);
    return dataDf;
  }

  async setData(df: DG.DataFrame): Promise<void> {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this.df = df;

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  private viewed: boolean = false;
  private tView?: DG.TableView;

  async destroyView(): Promise<void> {
    this.tView!.close();
    delete this.tView;
  }

  async buildView(): Promise<void> {
    console.debug('SequenceTranslator: OligoSdFileApp.buildView() ');

    this.tView = grok.shell.addTableView(this.df);
    this.tView.path = this.tView.basePath = 'func/SequenceTranslator.oligoSdFileApp';
  }
}