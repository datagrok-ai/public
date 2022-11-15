import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';
import {_package} from '../package';

//import {Observable, Subject} from 'rxjs';


export class NglGlService implements bio.NglGlServiceBase {
  readonly nglDiv: HTMLDivElement;
  private readonly ngl: any;

  hostDn: DG.DockNode;

  //private renderQueue: Subject<bio.NglGlRenderTask>;
  private readonly _queue: { key?: keyof any, task: bio.NglGlTask }[];
  private readonly _queueDict: { [key: keyof any]: bio.NglGlTask };

  constructor() {
    const r = window.devicePixelRatio;

    const defaultProps = {source: '(none:1);',};
    this.nglDiv = ui.div([], 'd4-ngl-viewer');
    this.nglDiv.style.width = `${300 / r}px`;
    this.nglDiv.style.height = `${300 / r}px`;

    // const blanket = ui.div();
    // blanket.style.position = 'absolute';
    // blanket.style.width = `${200 / r}px`;
    // blanket.style.height = `${200 / r}px`;

    // const windows = grok.shell.windows;
    // windows.showProperties = false;

    //@ts-ignore
    this.ngl = new NGL.Stage(this.nglDiv);
    // this.ngl.viewer.renderer.domElement.width = 300;
    // this.ngl.viewer.renderer.domElement.height = 300;

    this._queue = [];
    this._queueDict = {};
  }

  /* The flag allows _processQueue() on add item to the queue */
  private _busy: boolean = false;

  render(task: bio.NglGlTask, key?: keyof any): void {
    //console.debug('BSV: NglGlService.render() start ' + `name: ${name}`);

    if (key) {
      if (key in this._queueDict) {
        // remove outdated task from the queue
        const oldTaskI = this._queue.findIndex((item) => item.key == key);
        this._queue.splice(oldTaskI, 1);
      }
      this._queueDict[key] = task;
    }

    this._queue.push({key, task});

    if (!this._busy) {
      this._busy = true;

      const hostDiv = ui.div([this.nglDiv]);
      this.hostDn = grok.shell.tv.dockManager.dock(hostDiv, DG.DOCK_TYPE.RIGHT);
      console.debug('PTV: NglGlService dock()');

      window.setTimeout(async () => { await this._processQueue(); }, 0 /* next event cycle */);
    }
    //console.debug('BSV: NglGlService.render() end ' + `name: ${name}`);
  }

  private async _processQueue() {
    const {key, task} = this._queue.shift()!;
    if (key) delete this._queueDict[key];
    //console.debug('BSV: NglGlService._processQueue() start ' + `name: ${name}`);
    // this.ngl.viewer.scene.onBeforeRender = () => {
    //   const canvas = this.ngl.viewer.renderer.domElement;
    //   let context: WebGLRenderingContext = canvas.getContext('webgl');
    //   // context.clearColor(1, 0.7, 0.7, 1);
    // };

    let emptyPaintingSize: number = -1;
    let taskCompleted: boolean = false;

    const onAfterRenderHandler = async () => {
      console.debug('PTV: NglGlService onAfterRenderHandler()');
      if (emptyPaintingSize == -1) {
        emptyPaintingSize = this.ngl.viewer.renderer.domElement.toDataURL('image/png').length;
        let k = 11;
      } else {
        const currentPaintingSize = this.ngl.viewer.renderer.domElement.toDataURL('image/png').length;
        if (currentPaintingSize > emptyPaintingSize * 2) {
          console.debug('PTV: NglGlService taskCompleted');
          taskCompleted = true;
        }
      }

      if (taskCompleted) {
        this.ngl.viewer.signals.rendered.remove(onAfterRenderHandler);

        const canvas = this.ngl.viewer.renderer.domElement;
        await task.onAfterRender(canvas);

        if (this._queue.length > 0) {
          // window.clearInterval(timer);
          // TODO: Disable window.setInterval
          // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
          window.setTimeout(async () => { await this._processQueue(); }, 0 /* next event cycle */);
        } else {
          // release flag allowing _processQueue on add queue item
          this._busy = false;

          grok.shell.tv.dockManager.close(this.hostDn);
          // this.hostDn.detachFromParent();
          // this.hostDn = null;
          console.debug('PTV: NglGlService undock()');
        }
      }
    };

    // this.ngl.viewer.scene.onAfterRender = () => {
    //
    // };
    this.ngl.viewer.signals.rendered.add(onAfterRenderHandler);
    //const timer = window.setInterval(onAfterRenderHandler, 200);

    let stringBlob = new Blob([task.props.pdb], {type: 'text/plain'});
    // const pdb = await _package.files.readAsText('samples/1bdq.pdb');
    // let stringBlob = new Blob([pdb], {type: 'text/plain'});

    // this.ngl.loadFile(stringBlob, {ext: 'pdb'}).then((o: any) => {
    //   //await this.stage.loadFile(pdbStr).then(function (o: any) {
    //   let canvas = this.ngl.viewer.renderer.domElement;
    //   o.autoView();
    //   canvas = this.ngl.viewer.renderer.domElement;
    //   let a = 5;
    // });
    const r = window.devicePixelRatio;

    this.nglDiv.style.width = `${Math.floor(task.props.width) / r}px`;
    this.nglDiv.style.height = `${Math.floor(task.props.height) / r}px`;

    this.ngl.removeAllComponents();
    await this.ngl.loadFile(stringBlob, {ext: 'pdb'});
    await this.ngl.compList[0].addRepresentation('cartoon');
    await this.ngl.compList[0].autoView();

    let canvas = this.nglDiv.querySelector('canvas');
    canvas!.width = Math.floor(task.props.width);
    canvas!.height = Math.floor(task.props.height);
    await this.ngl.handleResize();

    // await new Promise((r) => setTimeout(r, 4000));
    let k = 11;

    // this.ngl.deck.setProps({
    //   onBeforeRender: ({gl}: { gl: WebGLRenderingContext }) => {
    //     const a: number = ((task.backColor >> 24) & 255) / 255;
    //     const r: number = ((task.backColor >> 16) & 255) / 255;
    //     const g: number = ((task.backColor >> 8) & 255) / 255;
    //     const b: number = ((task.backColor) & 255) / 255;
    //     gl.clearColor(r, g, b, a);
    //   },
    //   onAfterRender: ({gl}: { gl: WebGLRenderingContext }) => {
    //     //console.debug('BSV: NglGlService._processQueue() .onAfterRender() start ' + `name: ${name}`);
    //     task.onAfterRender(gl.canvas);
    //     //console.debug('BSV: NglGlService._processQueue() .onAfterRender() end ' + `name: ${name}`);

    //     if (this._queue.length > 0) {
    //       // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
    //       window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
    //     } else {
    //       // release flag allowing _processQueue on add queue item
    //       this._busy = false;
    //     }
    //   },
    // });

    // //@ts-ignore
    // const canvas: HTMLCanvasElement = this.ngl.deck.canvas!;
    // canvas.width = task.props.size.width;
    // canvas.height = task.props.size.height;

    //this.ngl.setProps(task.props); // should
    //console.debug('BSV: NglGlService._processQueue() end ' + `name: ${name}`);
  }

  public renderOnGridCell(
    gCtx: CanvasRenderingContext2D, bd: DG.Rect, gCell: DG.GridCell, canvas: CanvasImageSource
  ): void {
    gCtx.save();
    try {
      // Correction for vert scrolling happened between task and render, calculate bd.y directly
      bd.y = (gCell.gridRow - Math.floor(gCell.grid.vertScroll.min)) * gCell.grid.props.rowHeight +
        gCell.grid.colHeaderHeight;

      // Correction for horz scrolling happened between task and render, calculate bd.x directly
      let left: number = 0;
      for (let colI = 0; colI < gCell.gridColumn.idx; colI++) {
        const col: DG.GridColumn = gCell.grid.columns.byIndex(colI)!;
        left += col.visible ? col.width : 0;
      }
      if (gCell.grid.props.showRowHeader) {
        left += 20;
      }
      bd.x = left - gCell.grid.horzScroll.min;

      // gCtx.beginPath();
      // gCtx.rect(bd.x + 1, bd.y + 1, bd.width - 2, bd.height - 2);
      // gCtx.clip();

      gCtx.fillStyle = '#E0E0FF';
      gCtx.fillRect(bd.x + 1, bd.y + 1, bd.width - 2, bd.height - 2);

      //@ts-ignore
      gCtx.transform(bd.width / canvas.width, 0, 0, bd.height / canvas.height, bd.x, bd.y);

      gCtx.drawImage(canvas, 0 + 1, 0 + 1);
    } finally {
      gCtx.restore();
    }
  }
}