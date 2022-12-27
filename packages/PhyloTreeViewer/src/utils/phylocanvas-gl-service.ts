import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable, Subject} from 'rxjs';
import {PhylocanvasGL} from '@phylocanvas/phylocanvas.gl';
import {PhylocanvasGlServiceBase, PhylocanvasGlTask} from '@datagrok-libraries/bio';


export class PhylocanvasGlService implements PhylocanvasGlServiceBase {
  private readonly pcDiv: HTMLDivElement;
  private readonly pc: PhylocanvasGL;

  //private renderQueue: Subject<bio.PhylocanvasGlRenderTask>;
  private readonly _queue: { key?: keyof any, task: PhylocanvasGlTask }[];
  private readonly _queueDict: { [key: keyof any]: PhylocanvasGlTask };

  constructor() {
    const defaultProps = {source: '(none:1);'};
    this.pcDiv = ui.div([], {style: {backgroundColor: '#FFFFF0'}});
    this.pc = new PhylocanvasGL(this.pcDiv, defaultProps);
    this.pc.deck.setProps({useDevicePixels: true});

    this._queue = [];
    this._queueDict = {};
  }

  /* The flag allows _processQueue() on add item to the queue */
  private _busy: boolean = false;

  render(task: PhylocanvasGlTask, key?: keyof any): void {
    //console.debug('PhyloTreeViewer: PhylocanvasGlService.render() start ' + `name: ${name}`);

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
    //console.debug('PhyloTreeViewer: PhylocanvasGlService.render() end ' + `name: ${name}`);
  }

  private _processQueue() {
    const next = () => {
      if (this._queue.length > 0) {
        // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
        window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
      } else {
        // release flag allowing _processQueue on add queue item
        this._busy = false;
      }
    };

    const {key, task} = this._queue.shift()!;
    if (key) delete this._queueDict[key];
    //console.debug('PhyloTreeViewer: PhylocanvasGlService._processQueue() start ' + `name: ${name}`);

    this.pc.deck.setProps({
      onBeforeRender: ({gl}: { gl: WebGLRenderingContext }) => {
        const a: number = ((task.backColor >> 24) & 255) / 255;
        const r: number = ((task.backColor >> 16) & 255) / 255;
        const g: number = ((task.backColor >> 8) & 255) / 255;
        const b: number = ((task.backColor) & 255) / 255;
        gl.clearColor(r, g, b, a);
      },
      onAfterRender: ({gl}: { gl: WebGLRenderingContext }) => {
        // console.debug('PhyloTreeViewer: PhylocanvasGlService._processQueue() .onAfterRender() start ' +
        //   `name: ${name}`);
        task.onAfterRender(gl.canvas);
        // console.debug('PhyloTreeViewer: PhylocanvasGlService._processQueue() .onAfterRender() end ' +
        //   `name: ${name}`);

        next();
      },
    });

    //@ts-ignore
    const canvas: HTMLCanvasElement = this.pc.deck.canvas!;
    canvas.width = task.props.size.width;
    canvas.height = task.props.size.height;

    try {
      this.pc.setProps(task.props); // should
    } catch (err) {
      console.warn('PhyloTreeViewer: PhylocanvasGlService._processQueue() this.pc.setProps() catch -> next();');
      next();
    }
    //console.debug('PhyloTreeViewer: PhylocanvasGlService._processQueue() end ' + `name: ${name}`);
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
