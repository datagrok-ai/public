import * as DG from "datagrok-api/dg";
import {UaFilter} from "../filter2";
import {UaViewer} from "./abstract/ua-viewer";

export class UaDataFrameViewer extends UaViewer {
  dataFrame: DG.DataFrame;
  viewerFunction: Function;

  public constructor(name: string, dataFrame: DG.DataFrame, viewerFunction: Function, showName: boolean = true) {
    super(name, null, showName);
    this.dataFrame = dataFrame;
    this.viewerFunction = viewerFunction;
    this.reloadViewer();
  }

  setViewer(loader: any, host: HTMLDivElement) {
    host.appendChild(DG.Viewer.grid(this.dataFrame).root);
    host.removeChild(loader);
  }

}