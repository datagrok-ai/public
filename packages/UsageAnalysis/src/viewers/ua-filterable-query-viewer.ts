import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../filter';
import {UaQueryViewer} from './abstract/ua-query-viewer';
import {BehaviorSubject} from 'rxjs';

export class UaFilterableQueryViewer extends UaQueryViewer {
  protected filterSubscription: BehaviorSubject<UaFilter>;

  public constructor(filterSubscription: BehaviorSubject<UaFilter>, name: string, queryName: string,
    viewerFunction: Function, setStyle?: Function | null, staticFilter?: Object | null, viewer?: DG.Viewer) {
    super(name, queryName, viewerFunction, setStyle, staticFilter, null, viewer);
    this.filterSubscription = filterSubscription;
    this.filterSubscription.subscribe((filter) => this.reload(filter));
  }

  reload(filter: UaFilter) {
    this.filter = filter;
    this.reloadViewer();
  }
}
