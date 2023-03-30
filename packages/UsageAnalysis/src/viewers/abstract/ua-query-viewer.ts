// import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../../filter';
import {UaViewer} from './ua-viewer';


export abstract class UaQueryViewer extends UaViewer {
  queryName: string;
  viewerFunction: Function;
  staticFilter: Object = {};
  filter: Object = {};
  viewer: DG.Viewer;

  protected constructor(name: string, queryName: string, viewerFunction: Function,
    setStyle?: Function | null, staticFilter?: Object | null, filter?: UaFilter | null, viewer?: DG.Viewer) {
    super(name, setStyle);
    this.queryName = queryName;
    this.viewerFunction = viewerFunction;
    // if (staticFilter) this.staticFilter = staticFilter;
    if (filter) this.filter = filter;
    this.viewer = viewer as DG.Viewer;
    this.root.appendChild(this.viewer.root);
    this.root.appendChild(this.loader);
  }

  reloadViewer(staticFilter?: object) {
    this.loader.style.display = 'block';
    const filter = {...this.filter, ...staticFilter};
    grok.data.query('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
      if (dataFrame.columns.byName('count') != null)
        dataFrame.columns.byName('count').tags['format'] = '#';
      this.viewerFunction(dataFrame);
      this.viewer.dataFrame = dataFrame;
      this.viewer.setOptions({markerMinSize: 10, markerMaxSize: 30, color: 'user'});
      this.loader.style.display = 'none';
    });
  }
}
