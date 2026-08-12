import * as ui from 'datagrok-api/ui';
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
  errorDiv: HTMLElement | null = null;

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
    this.dataFrame = new Promise<DG.DataFrame>((resolve) => {
      this.errorDiv?.remove();
      this.errorDiv = null;
      this.root.appendChild(this.loader);
      this.loader.classList.add('ua-wait');
      const filter = {...this.filter, ...staticFilter};
      const load = this.queryName === undefined ? this.getDataFrame!() :
        grok.functions.call('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
          const countColumn = dataFrame.columns.byName('count');
          if (countColumn != null)
            countColumn.meta.format = '#';
          const userColumn = dataFrame.columns.byName('user');
          if (userColumn != null) {
            const users: {[key: string]: string} = {};
            for (const u of userColumn.categories)
              users[u] = colorHash.hex(u);
            userColumn.meta.colors.setCategorical(users);
          }
          return dataFrame;
        });
      load.then(this.postQuery.bind(this)).then(resolve).catch((e: any) => {
        this.loader.remove();
        this.errorDiv = ui.divText(`${this.name}: ${e?.message ?? e}`, 'd4-viewer-error');
        this.root.appendChild(this.errorDiv);
      });
    });
  }

  postQuery(dataFrame: DG.DataFrame): DG.DataFrame {
    if (this.processDataFrame)
      dataFrame = this.processDataFrame!(dataFrame);
    if (this.viewer != null) {
      this.viewer.root.remove();
      this.viewer = undefined;
    }
    this.viewer = this.createViewer(dataFrame);
    this.viewer.setOptions({
      markerMinSize: 10,
      markerMaxSize: 30,
    });
    this.root.appendChild(this.viewer.root);
    this.viewer.root.style.display = 'flex';
    this.loader.classList.remove('ua-wait');
    return dataFrame;
  }
}
