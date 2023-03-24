import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../../filter';
import {UaViewer} from './ua-viewer';


export abstract class UaQueryViewer extends UaViewer {
  queryName: string;
  viewerFunction: Function;
  staticFilter: Object = {};
  filter: Object = {};

  protected constructor(name: string, queryName: string, viewerFunction: Function,
    setStyle?: Function | null, staticFilter?: Object | null, filter?: UaFilter | null) {
    super(name, setStyle);

    this.queryName = queryName;
    this.viewerFunction = viewerFunction;

    if (staticFilter)
      this.staticFilter = staticFilter;
    if (filter)
      this.filter = filter;

    this.init();
  }

  setViewer(loader: any, host: HTMLDivElement) {
    const filter = {...this.filter, ...this.staticFilter};
    grok.data.query('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
      if (dataFrame.columns.byName('count') != null)
        dataFrame.columns.byName('count').tags['format'] = '#';
      let viewer: HTMLElement;
      if (dataFrame.rowCount > 0)
        viewer = this.viewerFunction(dataFrame);
      else
        viewer = ui.divText('No data', {style: {color: 'var(--red-3)', paddingBottom: '25px'}});
      const grid = DG.Viewer.grid(dataFrame);
      grid.props.allowColSelection = false;
      host.appendChild(viewer);
      host.removeChild(loader);
    });
  }

  init(): void {}
}
