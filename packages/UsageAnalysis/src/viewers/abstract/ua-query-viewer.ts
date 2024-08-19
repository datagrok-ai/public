// import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../../filter';
import {UaViewer} from './ua-viewer';
import ColorHash from 'color-hash';

const colorHash = new ColorHash();

export abstract class UaQueryViewer extends UaViewer {
  queryName?: string;
  createViewer: (t: DG.DataFrame) => DG.Viewer;
  getDataFrame?: () => Promise<DG.DataFrame>;
  processDataFrame?: (t: DG.DataFrame) => DG.DataFrame;
  dataFrame?: Promise<DG.DataFrame>;
  // staticFilter: Object = {};
  filter: Object = {};
  viewer: DG.Viewer | undefined;
  activated = false;

  protected constructor(name: string, options: {queryName?: string, createViewer: (t: DG.DataFrame) => DG.Viewer,
    setStyle?: Function | null, staticFilter?: Object | null, filter?: UaFilter | null,
    getDataFrame?: () => Promise<DG.DataFrame>, processDataFrame?: (t: DG.DataFrame) => DG.DataFrame, activated?: boolean}) {
    super(name, options.setStyle);
    this.queryName = options.queryName;
    this.createViewer = options.createViewer;
    this.getDataFrame = options.getDataFrame;
    this.processDataFrame = options.processDataFrame;
    this.activated = options.activated ?? false;
    // if (staticFilter) this.staticFilter = staticFilter;
    if (options.filter)
      this.filter = options.filter;
    this.root.appendChild(this.loader);
  }

  reloadViewer(staticFilter?: object) {
    // console.log('reloading');
    this.dataFrame = new Promise<DG.DataFrame>((resolve, reject) => {
      this.loader.classList.add('ua-wait');
      const filter = {...this.filter, ...staticFilter};
      // console.log(this.queryName);
      if (this.queryName === undefined) {
        this.getDataFrame!().then(this.postQuery.bind(this)).then(resolve.bind(this));
        return;
      }
      grok.functions.call('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
        if (dataFrame.columns.byName('count') != null)
          dataFrame.columns.byName('count').meta.format = '#';
        const userColumn = dataFrame.columns.byName('user');
        if (userColumn != null) {
          const users: {[key: string]: string} = {};
          userColumn.categories.forEach((u: string) => {users[u] = colorHash.hex(u);});
          userColumn.meta.colors.setCategorical(users);
        }
        this.postQuery(dataFrame);
        return dataFrame;
      }).then(resolve.bind(this));
    });
  }

  postQuery(dataFrame: DG.DataFrame): DG.DataFrame {
    if (this.processDataFrame)
      dataFrame = this.processDataFrame!(dataFrame);
    this.viewer ??= this.createViewer(dataFrame);
    this.viewer.dataFrame = dataFrame;
    this.viewer.setOptions({
      markerMinSize: 10,
      markerMaxSize: 30,
    });
    this.root.appendChild(this.viewer!.root);
    this.viewer!.dataFrame = dataFrame ?? this.viewer!.dataFrame;
    this.viewer!.root.style.display = 'flex';
    this.loader.classList.remove('ua-wait');
    return dataFrame;
  }
}
