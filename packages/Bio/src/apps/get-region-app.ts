import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS, positionSeparator} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {IWebLogoViewer} from '@datagrok-libraries/bio/src/viewers/web-logo';

import {_package} from '../package';

const defaultData: GetRegionAppData = {
  df: DG.DataFrame.fromCsv(`seq,value
ATCCGTCGT,0.5
TGTTCGTCA,0.4
ATGGTCGTA,0.7
ATCCGTGCA,0.1`),
  colName: 'seq',
  positionNames: ['1', '1A', '1C', '2', '4', '4A', '4B', '5', '6'].join(positionSeparator),
  regions: [
    {name: 'first region', start: '1', end: '2'},
    {name: 'second region', start: '1C', end: '4'},
    {name: 'overlapping second', start: '1C', end: '4A'},
    {name: 'whole sequence', start: '1', end: '6'},
    {name: 'bad start', start: '0', end: '6'},
    {name: 'bad end', start: '1', end: '4C'},
    {name: 'bad start & end', start: '0', end: '4C'},
  ]
};

export type GetRegionAppData = {
  df: DG.DataFrame, colName: string,
  positionNames?: string, regions?: { name: string, start: string, end: string }[]
};

export class GetRegionApp {
  view: DG.TableView;
  data!: GetRegionAppData;

  constructor(
    private readonly urlParams: URLSearchParams,
    private readonly funcName: string
  ) {}

  async init(data?: GetRegionAppData): Promise<void> {
    this.data = data ?? defaultData;
    const seqCol = this.data.df.getCol(this.data.colName);
    if (!!this.data.positionNames) seqCol.setTag(TAGS.positionNames, this.data.positionNames);
    if (!!this.data.regions) seqCol.setTag(TAGS.regions, JSON.stringify(this.data.regions));

    await this.buildView();
  }

  // -- View --

  async buildView(): Promise<void> {
    // To allow showing a WebLogoViewer
    await grok.data.detectSemanticTypes(this.data.df);

    this.view = grok.shell.addTableView(this.data.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.funcName}`;

    const viewer: DG.Viewer & IWebLogoViewer = (await this.view.dataFrame.plot
      .fromType('WebLogo', {sequenceColumnName: this.data.colName})) as DG.Viewer & IWebLogoViewer;
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, 'WebLogo', 0.35);
  }
}
