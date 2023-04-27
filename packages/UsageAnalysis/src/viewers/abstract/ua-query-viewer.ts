// import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../../filter';
import {UaViewer} from './ua-viewer';
import ColorHash from 'color-hash';

const colorHash = new ColorHash();

export abstract class UaQueryViewer extends UaViewer {
  queryName?: string;
  viewerFunction: (t: DG.DataFrame) => DG.Viewer;
  getDataFrame?: () => Promise<DG.DataFrame>;
  dataFrame?: Promise<DG.DataFrame>;
  staticFilter: Object = {};
  filter: Object = {};
  viewer: DG.Viewer | undefined;

  protected constructor(name: string, options: {queryName?: string, viewerFunction: (t: DG.DataFrame) => DG.Viewer,
    setStyle?: Function | null, staticFilter?: Object | null, filter?: UaFilter | null,
    getDataFrame?: () => Promise<DG.DataFrame>}) {
    super(name, options.setStyle);
    this.queryName = options.queryName;
    this.viewerFunction = options.viewerFunction;
    this.getDataFrame = options.getDataFrame;
    // if (staticFilter) this.staticFilter = staticFilter;
    if (options.filter)
      this.filter = options.filter;
    this.root.appendChild(this.loader);
  }

  reloadViewer(staticFilter?: object) {
    // console.log('reloading');
    this.dataFrame = new Promise<DG.DataFrame>((resolve, reject) => {
      this.loader.style.display = 'block';
      if (this.viewer !== undefined)
        this.viewer.root.style.display = 'none';
      const filter = {...this.filter, ...staticFilter};
      console.log(this.queryName);
      if (this.queryName === undefined) {
        this.getDataFrame!().then(this.postQuery.bind(this)).then(resolve.bind(this));
        return;
      }
      grok.data.query('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
        if (dataFrame.columns.byName('count') != null)
          dataFrame.columns.byName('count').tags['format'] = '#';
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
    this.viewer = this.viewerFunction(dataFrame);
    this.root.appendChild(this.viewer!.root);
    this.viewer!.dataFrame = dataFrame ?? this.viewer!.dataFrame;
    this.viewer!.root.style.display = 'flex';
    this.loader.style.display = 'none';
    return dataFrame;
  }
}
