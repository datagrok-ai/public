import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../filter';
import {UaQueryViewer} from './abstract/ua-query-viewer';
import {BehaviorSubject} from 'rxjs';

export class UaFilterableQueryViewer extends UaQueryViewer {
  protected filterSubscription: BehaviorSubject<UaFilter>;

  public constructor(options: {filterSubscription: BehaviorSubject<UaFilter>, name: string, queryName?: string,
    viewerFunction: (t: DG.DataFrame) => DG.Viewer, setStyle?: Function | null, staticFilter?: Object | null,
    getDataFrame?: () => Promise<DG.DataFrame>}) {
    super(options.name, {queryName: options.queryName, viewerFunction: options.viewerFunction,
      setStyle: options.setStyle, staticFilter: options.staticFilter, getDataFrame: options.getDataFrame});
    this.filterSubscription = options.filterSubscription;
    this.filterSubscription.subscribe((filter) => this.reload(filter));
  }

  reload(filter: UaFilter) {
    this.filter = filter;
    this.reloadViewer();
  }
}
