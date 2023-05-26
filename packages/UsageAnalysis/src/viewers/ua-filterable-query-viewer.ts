import * as DG from 'datagrok-api/dg';

import {UaFilter} from '../filter';
import {UaQueryViewer} from './abstract/ua-query-viewer';
import {BehaviorSubject} from 'rxjs';

export class UaFilterableQueryViewer extends UaQueryViewer {
  protected filterSubscription: BehaviorSubject<UaFilter>;

  public constructor(options: {filterSubscription: BehaviorSubject<UaFilter>, name: string, queryName?: string,
    createViewer: (t: DG.DataFrame) => DG.Viewer, setStyle?: Function | null, staticFilter?: Object | null,
    getDataFrame?: () => Promise<DG.DataFrame>, processDataFrame?: (t: DG.DataFrame) => DG.DataFrame, activated?: boolean}) {
    super(options.name, {queryName: options.queryName, createViewer: options.createViewer,
      setStyle: options.setStyle, staticFilter: options.staticFilter, getDataFrame: options.getDataFrame,
      processDataFrame: options.processDataFrame, activated: options.activated});
    this.filterSubscription = options.filterSubscription;
    this.filterSubscription.subscribe((filter) => this.reload(filter));
  }

  reload(filter: UaFilter) {
    this.filter = filter;
    if (this.activated)
      this.reloadViewer();
  }
}
