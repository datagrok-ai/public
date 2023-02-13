import * as DG from 'datagrok-api/dg';
import {UaViewer} from './abstract/ua-viewer';

export class UaDataFrameViewer extends UaViewer {
  dataFrame: DG.DataFrame;
  viewerFunction: Function;

  public constructor(name: string, dataFrame: DG.DataFrame, viewerFunction: Function, showName: boolean = true) {
    super(name, null, showName);
    this.dataFrame = dataFrame;
    this.viewerFunction = viewerFunction;
    this.reloadViewer();
  }

  setViewer(loader: any, host: HTMLDivElement, nameDiv: HTMLElement) {
    const grid = DG.Viewer.grid(this.dataFrame);
    host.appendChild(grid.root);
    grid.autoSize(600, 400, 200, 100);
    host.removeChild(loader);
  }
}
