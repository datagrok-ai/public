import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../filter';
import {UaQueryViewer} from './abstract/ua-query-viewer';
import {BehaviorSubject} from 'rxjs';

export class UaFilterableQueryViewer extends UaQueryViewer {
  protected filterSubscription: BehaviorSubject<UaFilter>;

  public constructor(filterSubscription: BehaviorSubject<UaFilter>, name: string, queryName: string,
    viewerFunction: (t: DG.DataFrame) => DG.Viewer, setStyle?: Function | null, staticFilter?: Object | null) {
    super(name, queryName, viewerFunction, setStyle, staticFilter, null);
    this.filterSubscription = filterSubscription;
    this.filterSubscription.subscribe((filter) => this.reload(filter));
  }

  reload(filter: UaFilter) {
    this.filter = filter;
    this.reloadViewer();
  }
}
