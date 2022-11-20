import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

export class NglForGridTestApp {
  df: DG.DataFrame;

  constructor() {

  }

  async init(): Promise<void> {
    const dfCsv: string = await _package.files.readAsText('pdb_data.csv');
    this.df = DG.DataFrame.fromCsv(dfCsv);
    const tv = grok.shell.addTableView(this.df);
    tv.path = tv.basePath = '/func/BiostructureViewer.nglForGridTestApp';
  }
}
