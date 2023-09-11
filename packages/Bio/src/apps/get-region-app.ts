import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS, positionSeparator} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {IWebLogoViewer} from '@datagrok-libraries/bio/src/viewers/web-logo';

import {_package} from '../package';

const csv = `seq,value
ATCCGTCGT,0.5
TGTTCGTCA,0.4
ATGGTCGTA,0.7
ATCCGTGCA,0.1`;

const positionNames = ['1', '1A', '1C', '2', '4', '4A', '4B', '5', '6'].join(positionSeparator);

const regions = [
  {name: 'first region', start: '1', end: '2'},
  {name: 'second region', start: '1C', end: '4'},
  {name: 'overlapping second', start: '1C', end: '4A'},
  {name: 'whole sequence', start: '1', end: '6'},
  {name: 'bad start', start: '0', end: '6'},
  {name: 'bad end', start: '1', end: '4C'},
  {name: 'bad start & end', start: '0', end: '4C'},
];

export class GetRegionApp {
  df: DG.DataFrame;
  view: DG.TableView;

  constructor(
    private readonly urlParams: URLSearchParams,
    private readonly funcName: string
  ) {}

  async init(): Promise<void> {
    this.df = DG.DataFrame.fromCsv(csv);
    const seqCol = this.df.getCol('seq');
    seqCol.setTag(TAGS.positionNames, positionNames);
    seqCol.setTag(TAGS.regions, JSON.stringify(regions));

    await this.buildView();
  }

  // -- View --

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.funcName}`;

    const viewer: DG.Viewer & IWebLogoViewer = (await this.view.dataFrame.plot
      .fromType('WebLogo')) as DG.Viewer & IWebLogoViewer;
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, 'WebLogo', 0.35);
  }
}
