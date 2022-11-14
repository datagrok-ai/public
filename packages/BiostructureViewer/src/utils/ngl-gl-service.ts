import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';
//import {Observable, Subject} from 'rxjs';



export class NglGlService implements bio.NglGlServiceBase {
  private readonly nglDiv: HTMLDivElement;
  private readonly ngl: any;

  //private renderQueue: Subject<bio.NglGlRenderTask>;
  private readonly _queue: { key?: keyof any, task: bio.NglGlTask }[];
  private readonly _queueDict: { [key: keyof any]: bio.NglGlTask };

  constructor() {
    const defaultProps = {source: '(none:1);',};
    this.nglDiv = ui.div([], 'd4-ngl-viewer');
    this.nglDiv.style.width = '200px';
    this.nglDiv.style.height = '200px';

    //@ts-ignore
    this.ngl = new NGL.Stage(this.nglDiv);
    this.ngl.viewer.renderer.domElement.width = 200;
    this.ngl.viewer.renderer.domElement.height = 200;
    
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
      window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
    }
    //console.debug('BSV: NglGlService.render() end ' + `name: ${name}`);
  }

  private _processQueue() {
    const {key, task} = this._queue.shift()!;
    if (key) delete this._queueDict[key];
    //console.debug('BSV: NglGlService._processQueue() start ' + `name: ${name}`);
    this.ngl.viewer.scene.onBeforeRender = () => {
      const canvas = this.ngl.viewer.renderer.domElement;
      let context: WebGLRenderingContext = canvas.getContext('webgl');
      context.clearColor(1, 0.7, 0.7, 1);
    };
    this.ngl.viewer.scene.onAfterRender = () => {
      const canvas = this.ngl.viewer.renderer.domElement;
      task.onAfterRender(canvas);
      if (this._queue.length > 0) {
        // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
        window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
      } else {
        // release flag allowing _processQueue on add queue item
        this._busy = false;
      }
    };

    let stringBlob = new Blob([task.props.pdb], { type: 'text/plain' });
    this.ngl.loadFile(stringBlob, { ext: "pdb" }).then((o: any) => {
      //await this.stage.loadFile(pdbStr).then(function (o: any) {
        let canvas = this.ngl.viewer.renderer.domElement;
        o.autoView();
        canvas = this.ngl.viewer.renderer.domElement;
        let a= 5;
    });


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
      bd.x = left - gCell.grid.horzScroll.min;

      gCtx.beginPath();
      gCtx.rect(bd.x + 1, bd.y + 1, bd.width - 2, bd.height - 2);
      gCtx.clip();

      // gCtx.fillStyle = '#FFE0E0';
      // gCtx.fillRect(bd.x + 1, bd.y + 1, bd.width - 2, bd.height - 2);

      gCtx.drawImage(canvas, bd.x + 1, bd.y + 1);
    } finally {
      gCtx.restore();
    }
  }
}