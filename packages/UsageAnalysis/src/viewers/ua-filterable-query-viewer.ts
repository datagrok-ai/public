import {UaFilter} from "../filter2";
import {UaQueryViewer} from "./abstract/ua-query-viewer";
import { BehaviorSubject  } from "rxjs";

export class UaFilterableQueryViewer extends UaQueryViewer {
  private filterSubscription: BehaviorSubject<UaFilter>;

  public constructor(filterSubscription: BehaviorSubject<UaFilter>, name: string, queryName: string, viewerFunction: Function, setStyle?: Function | null, staticFilter?: Object, showName?: boolean) {
    super(name, queryName, viewerFunction, setStyle, staticFilter, null, showName);
    this.filterSubscription = filterSubscription;
    this.filterSubscription.subscribe((filter) => this.reload(filter));
    this.reload(this.filterSubscription.getValue());
  }

  reload(filter: UaFilter) {
    this.filter = filter;
    this.reloadViewer();
  }

}