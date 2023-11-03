import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import * as DG from 'datagrok-api/dg';


export class LogView extends UaView {
  expanded: { [key: string]: boolean } = {f: true, l: true};

  // @ts-ignore
  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Log';
  }

  async initViewers(): Promise<void> {
    const logViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Log',
      queryName: 'LogTail',
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.grid(t);
        return viewer;
      }});

    this.viewers.push(logViewer);
    this.root.append(logViewer.root);
  }
}
