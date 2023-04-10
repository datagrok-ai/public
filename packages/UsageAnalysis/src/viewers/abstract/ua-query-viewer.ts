// import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../../filter';
import {UaViewer} from './ua-viewer';
import ColorHash from 'color-hash';

const colorHash = new ColorHash();

export abstract class UaQueryViewer extends UaViewer {
  queryName: string;
  viewerFunction: Function;
  staticFilter: Object = {};
  filter: Object = {};
  viewer: DG.Viewer | undefined;

  protected constructor(name: string, queryName: string, viewerFunction: Function,
    setStyle?: Function | null, staticFilter?: Object | null, filter?: UaFilter | null) {
    super(name, setStyle);
    this.queryName = queryName;
    this.viewerFunction = viewerFunction;
    // if (staticFilter) this.staticFilter = staticFilter;
    if (filter) this.filter = filter;
    this.root.appendChild(this.loader);
  }

  reloadViewer(staticFilter?: object) {
    this.loader.style.display = 'block';
    if (this.viewer !== undefined)
      this.viewer.root.style.display = 'none';
    const filter = {...this.filter, ...staticFilter};
    grok.data.query('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
      if (dataFrame.columns.byName('count') != null)
        dataFrame.columns.byName('count').tags['format'] = '#';
      const userColumn = dataFrame.columns.byName('user');
      if (userColumn != null) {
        const users: {[key: string]: string} = {};
        userColumn.categories.forEach((u: string) => {users[u] = colorHash.hex(u);});
        userColumn.meta.colors.setCategorical(users);
      }
      this.viewer = this.viewerFunction(dataFrame);
      this.root.appendChild(this.viewer!.root);
      this.viewer!.dataFrame = dataFrame;
      this.viewer!.root.style.display = 'flex';
      this.loader.style.display = 'none';
    });
  }
}
