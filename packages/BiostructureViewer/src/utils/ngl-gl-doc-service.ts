import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as NGL from 'NGL';
import {NglGlServiceBase, NglGlTask} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';
import {SignalBinding} from 'signals';

import {_package} from '../package';

const TASK_TIMEOUT: number = 2000;

export class NglGlDocService implements NglGlServiceBase {
  private readonly nglDiv: HTMLDivElement;
  private readonly ngl: NGL.Stage;

  private readonly hostDiv: HTMLDivElement;

  private readonly _queue: { key?: keyof any, task: NglGlTask, dt: number }[];
  private readonly _queueDict: { [key: keyof any]: NglGlTask };
  /** The flag allows {@link _processQueue}() on add item to the queue with {@link render}() */
  private _busy: boolean = false;
  private nglRenderedBinding: SignalBinding<any>;

  constructor() {
    // const r = window.devicePixelRatio;
    this.nglDiv = ui.div([], 'd4-ngl-viewer');
    this.ngl = new NGL.Stage(this.nglDiv);

    // The single NGL component
    this.hostDiv = ui.box(this.nglDiv);
    // this.hostDiv.style.display = 'none'; // Disables drawing at all
    this.hostDiv.style.position = 'absolute';
    this.hostDiv.style.left = '0px';
    this.hostDiv.style.right = '0px';
    this.hostDiv.style.width = '300px';
    this.hostDiv.style.height = '300px';
    this.hostDiv.style.visibility = 'hidden';
    document.body.appendChild(this.hostDiv);

    this._queue = [];
    this._queueDict = {};

    window.setInterval(() => { this._sweepQueue(); }, 200);
  }

  render(task: NglGlTask, key?: keyof any): void {
    _package.logger.debug('NglGlDocService.render() start ' + `key: ${key?.toString()}`);

    if (key !== undefined) {
      if (key in this._queueDict) {
        // remove outdated task from the queue
        const oldTaskI = this._queue.findIndex((item) => item.key === key);
        this._queue.splice(oldTaskI, 1);
      }
      this._queueDict[key] = task;
    }

    this._queue.push({key, task, dt: Date.now()});

    if (!this._busy) {
      this._busy = true;

      // TODO: Use requestAnimationFrame()
      window.setTimeout(async () => { await this._processQueue(); }, 0 /* next event cycle */);
    }
    //_package.logger.debug('NglGlDocService.render() end ' + `name: ${name}`);
  }

  private async _processQueue() {
    const queueItem = this._queue.shift();
    if (!queueItem) {
      _package.logger.error('BsV: NglGlDocService._processQueue() queueItem = undefined ');
      return;
    } // in case of empty queue

    const {key, task} = queueItem!;
    if (key) {
      _package.logger.debug('NglGlDocService._processQueue() ' + `key = ${key.toString()}`);
      delete this._queueDict[key];
    }
    try {
      const dpr = window.devicePixelRatio;

      // TODO: Convert string to Blob once converting PDB string column to Blob
      const stringBlob = new Blob([task.props.pdb], {type: 'text/plain'});

      this.task = undefined;
      this.key = undefined;
      // TODO: Use canvas size switching 0/1 px to required

      const v0 = this.ngl.removeAllComponents();

      this.nglDiv.style.width = `${Math.floor(task.props.width) / dpr}px`;
      this.nglDiv.style.height = `${Math.floor(task.props.height) / dpr}px`;
      this.ngl.viewer.setSize(task.props.width / dpr, task.props.height / dpr);
      const canvas = this.ngl.viewer.renderer.domElement;
      canvas.width = Math.floor(task.props.width);
      canvas.height = Math.floor(task.props.height);

      await this.ngl.loadFile(stringBlob, {ext: 'pdb', defaultRepresentation: true});
      this.ngl.compList[0].addRepresentation('cartoon');
      this.ngl.compList[0].autoView();

      this.task = task;
      this.key = key;
      this.emptyCanvasHash = DG.StringUtils.hashCode(canvas.toDataURL());
      this.nglRenderedBinding = this.ngl.viewer.signals.rendered.add(this.onNglRendered, this);
      this.ngl.viewer.requestRender();
    } catch (err: any) {
      const errMsg: string = err instanceof Error ? err.message : err.toString();
      _package.logger.error(`BsV:NglGlDocService._processQueue() no rethrown error: ${errMsg}`, undefined,
        err instanceof Error ? err.stack : undefined);
      //throw err; // Do not throw to prevent disabling queue handler
    }
  }

  /** Sweep queue for stalled tasks */
  private _sweepQueue(): void {
    const nowDt: number = Date.now();
    let swept: boolean = false;

    for (let qI = this._queue.length - 1; qI >= 0; qI--) {
      const {key, task: _task, dt} = this._queue[qI];
      if ((nowDt - dt) > TASK_TIMEOUT) {
        // stalled task
        this._queue.splice(qI);
        delete this._queueDict[key!];
        swept = true;
      }
    }

    this._busy = this._queue.length > 0;
    if (!this._busy && swept) {
      // Some tasks had been swept
    }
  }

  private emptyCanvasHash?: number = undefined;
  private task?: NglGlTask = undefined;
  private key?: keyof any = undefined;

  private async onNglRendered(): Promise<void> {
    if (this.task === undefined) return;
    _package.logger.debug('NglGlDocService.onNglRendered() ' + `key = ${JSON.stringify(this.key)}`);
    try {
      const canvas = this.ngl.viewer.renderer.domElement;
      const canvasHash = DG.StringUtils.hashCode(canvas.toDataURL());
      if (canvasHash == this.emptyCanvasHash) { // render is not ready yet
        _package.logger.debug(`NglGlDocService.onNglRendered() empty canvas`);
        return;
      }
      _package.logger.debug(`NglGlDocService.onNglRendered() ` +
        `this.emptyCanvasHash = ${this.emptyCanvasHash}, canvasHash = ${canvasHash}`);

      this.nglRenderedBinding.detach();

      await this.task.onAfterRender(canvas);
      this.task = undefined;
      this.key = undefined;
      this.emptyCanvasHash = undefined;

      if (this._queue.length > 0) {
        // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
        window.setTimeout(async () => { await this._processQueue(); }, 0 /* next event cycle */);
      } else {
        // release flag allowing _processQueue on add queue item
        this._busy = false;
      }
    } catch (err: any) {
      const errMsg: string = err instanceof Error ? err.message : err.toString();
      _package.logger.error(`NglGlDocService.onNglRendered() no rethrown error : ${errMsg}`, undefined,
        err instanceof Error ? err.stack : undefined);
      //throw err; // Do not throw to prevent disabling event/signal handler
    }
  }

  public renderOnGridCell(
    gCtx: CanvasRenderingContext2D, bd: DG.Rect, gCell: DG.GridCell, canvas: CanvasImageSource,
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

      // gCtx.fillStyle = '#E0E0FF';
      // gCtx.fillRect(bd.x + 1, bd.y + 1, bd.width - 2, bd.height - 2);

      const cw = 'width' in canvas && (
        canvas.width instanceof SVGAnimatedLength ? canvas.width.baseVal.value : canvas.width as number);
      const ch = 'height' in canvas && (
        canvas.height instanceof SVGAnimatedLength ? canvas.height.baseVal.value : canvas.height as number);
      if (!cw || !ch) throw new Error('NglGlDocService.renderOnGridCell() canvas size is not available');
      gCtx.transform(bd.width / cw, 0, 0, bd.height / ch, bd.x, bd.y);

      gCtx.drawImage(canvas, 0 + 1, 0 + 1);
    } finally {
      gCtx.restore();
    }
  }
}
