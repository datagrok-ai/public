import {UaQueryViewer} from "./abstract/ua-query-viewer";
import {UaFilter} from "../filter2";

export class UaDataFrameQueryViewer extends UaQueryViewer {

  public constructor(name: string, queryName: string, viewerFunction: Function, setStyle?: Function,
                     staticFilter?: Object, filter?: UaFilter, showName: boolean = true) {
    super(name, queryName, viewerFunction, setStyle, staticFilter, filter, showName);
  }

  init(): void {
    this.reloadViewer();
  }

}