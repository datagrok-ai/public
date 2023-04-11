import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {delay} from '@datagrok-libraries/utils/src/test';

const dataFn = 'data/sample_FASTA_DNA.csv';

export async function demoBio01UI(funcPath: string) {
  let view: DG.TableView;
  let df: DG.DataFrame;

  async function step01() {
    grok.shell.info(`Loading DNA notation 'fasta'.`);
    const pi = DG.TaskBarProgressIndicator.create('Loading \'fasta\' data...');
    try {
      df = await _package.files.readCsv(dataFn);
      view = grok.shell.addTableView(df);
      view.path = view.basePath = funcPath;
    } finally {
      pi.close();
    }
  }

  async function step02() {
    grok.shell.info('Similarity search analysis for sequences');
    const pi = DG.TaskBarProgressIndicator.create('Similarity search...');
    try {
      const simViewer = await df.plot.fromType('Sequence Similarity Search') as DG.Viewer;
      view.dockManager.dock(simViewer, DG.DOCK_TYPE.RIGHT, null, 'Similarity search', 0.35);
    } finally {
      pi.close();
    }
  }

  async function step03() {
    grok.shell.info('Diversity search analysis for sequences');
    const pi = DG.TaskBarProgressIndicator.create('Diversity sequences ...');
    try {
      const divViewer = await df.plot.fromType('Sequence Diversity Search') as DG.Viewer;
      view.dockManager.dock(divViewer, DG.DOCK_TYPE.DOWN, null, 'Diversity search', 0.27);
    } finally {
      pi.close();
    }
  }

  async function step04a() {
    grok.shell.info('Current row  3');
    df.currentRowIdx = 3;
  }

  async function step04b() {
    grok.shell.info('Current row  7');
    df.currentRowIdx = 7;
  }

  Promise.resolve()
    .then(async () => { await step01(); })
    .then(async () => { await delay(1600); })
    .then(async () => { await step02(); })
    .then(async () => { await delay(1600); })
    .then(async () => { await step03(); })
    .then(async () => { await delay(1600); })
    .then(async () => { await step04a(); })
    .then(async () => { await delay(1600); })
    .then(async () => { await step04b(); });
}
