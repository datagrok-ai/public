import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {SignalBinding} from 'signals';
import * as NGL from 'NGL';

import {NglGlServiceBase, NglGlTask} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';

import {errInfo} from './err-info';

import {_package} from '../package';
import {delay, testEvent} from '@datagrok-libraries/utils/src/test';
import {Subject} from 'rxjs';
import {awaitNgl} from '../viewers/ngl-viewer-utils';


const TASK_TIMEOUT: number = 2000;
const NGL_ERROR_LIMIT: number = 3;
const NGL_TRY_LIMIT: number = 3;

export class NglGlDocService extends NglGlServiceBase {
  private readonly nglDiv: HTMLDivElement;

  private ngl: NGL.Stage | null = null;
  private nglErrorCount: number = 0;

  private readonly hostDiv: HTMLDivElement;

  private readonly _queue: {
    key?: keyof any,
    task: NglGlTask,
    tryCount: number,
    dt: number
  }[];
  private readonly _queueDict: {
    [key: keyof any]: NglGlTask
  };
  /** The flag allows {@link _processQueue}() on add item to the queue with {@link render}() */
  private _busy: boolean = false;
  private nglRenderedBinding: SignalBinding<any>;

  constructor() {
    super();

    // const r = window.devicePixelRatio;
    this.nglDiv = ui.div([], 'd4-ngl-viewer');

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

  public get errorCount(): number { return this.nglErrorCount; }

  public async reset(): Promise<void> {
    if (!this.ngl) return;
    this.ngl.dispose();
    this.ngl = null;
    $(this.nglDiv).empty();
    this.nglErrorCount = 0;
  }

  render(task: NglGlTask, key?: keyof any, tryCount: number = 0): void {
    _package.logger.debug('NglGlDocService.render() start ' + `key: ${key?.toString()}`);

    if (key !== undefined) {
      if (key in this._queueDict) {
        // remove outdated task from the queue
        const oldTaskI = this._queue.findIndex((item) => item.key === key);
        this._queue.splice(oldTaskI, 1);
      }
      _package.logger.debug(`NglGlDocService.render() _queueDict[ key: ${key?.toString()} ] = <task>`);
      this._queueDict[key] = task;
    }

    _package.logger.debug('NglGlDocService.render() _queue.push(), ' + `key: ${key?.toString()}`);
    this._queue.push({key, task, tryCount, dt: Date.now()});

    if (!this._busy) {
      this._busy = true;

      // TODO: Use requestAnimationFrame()
      _package.logger.debug('NglGlDocService.render(), window.setTimeout() -> this._processQueue() ');
      window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
    }
  }

  private _processQueue(): void {
    if (this._queue.length === 0) return;
    _package.logger.debug(`NglGlDocService._processQueue(), ` +
      `queue: ${JSON.stringify(this._queue.map((t) => t.key))}`);

    let nglRenderBinding: SignalBinding<any> | null = null;

    const finallyProcessQueue = () => {
      if (nglRenderBinding !== null) nglRenderBinding.detach();
      if (this._queue.length > 0) {
        // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
        _package.logger.debug('NglGlDocService._processQueue().finallyProcessQueue(), ' +
          'window.setTimeout() -> this._processQueue() ');
        window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
      } else {
        // release flag allowing _processQueue on add queue item
        _package.logger.debug('NglGlDocService._processQueue().finallyProcessQueue(), this._busy = false');
        this._busy = false;
      }
    };

    const queueItem = this._queue.shift();
    if (!queueItem) {
      _package.logger.error('NglGlDocService._processQueue() queueItem = undefined ');
      finallyProcessQueue();
      return;
    } // in case of empty queue

    const {key, task, tryCount, dt} = queueItem;
    if (key !== undefined) {
      _package.logger.debug('NglGlDocService._processQueue() ' + `key: ${key.toString()}`);
      delete this._queueDict[key];
    }
    if (tryCount > NGL_TRY_LIMIT) {
      _package.logger.warning('NglGlDocService._processQueue(), skip task, ' +
        ` key: ${key?.toString()}, tryCount: ${tryCount}`);
      finallyProcessQueue();
      return;
    }

    let emptyCanvasHash: number;
    let handled: boolean = false;

    const timeoutHandle = window.setTimeout(() => {
      if (!handled) {
        _package.logger.warning('NglGlDocService._processQueue().timeoutHandle() not handled, ' +
          `key: ${key?.toString}`);
        this.nglErrorCount += 1;
        this.render(task, key, tryCount + 1); // return task to the queue
        finallyProcessQueue();
      }
    }, 1000);

    const renderHandler = () => {
      _package.logger.warning('NglGlDocService._processQueue().handler(), ' + `key: ${key?.toString}`);
      if (this.nglOnRendered(key, task, emptyCanvasHash)) {
        handled = true;
        window.clearTimeout(timeoutHandle);
        finallyProcessQueue();
      }
    };

    this.nglRender(key, task, renderHandler).then((res) => {
      [nglRenderBinding, emptyCanvasHash] = res;
    })
      .catch((err: any) => {
        // Not waiting timeout on error
        const [errMsg, errStack] = errInfo(err);
        _package.logger.error(errMsg, undefined, errStack);
        window.clearTimeout(timeoutHandle);
        this.nglErrorCount += 1;
        this.render(task, key, tryCount + 1); // return task to the queue
        finallyProcessQueue();
      });
  }

  private async nglRender(
    key: any, task: NglGlTask, renderHandler: () => void,
  ): Promise<[SignalBinding<any>, number]> {
    const callLogPrefix = 'NglGlDocService.nglRender()';
    _package.logger.debug(`${callLogPrefix}, ` + `key: ${key?.toString()}`);
    const dpr = window.devicePixelRatio;

    // TODO: Convert string to Blob once converting PDB string column to Blob
    const stringBlob = new Blob([task.props.pdb], {type: 'text/plain'});

    // TODO: Use canvas size switching 0/1 px to required

    // if (key === 1) throw new Error('NglGlDocService: Test error');

    if (this.nglErrorCount > NGL_ERROR_LIMIT) {
      _package.logger.warning(`${callLogPrefix}, ` + 'recreate NGL.Stage, ' +
        `nglTimeoutCount = ${this.nglErrorCount} > ${NGL_ERROR_LIMIT}`);
      await this.reset();
    }
    if (!this.ngl) {
      this.ngl = new NGL.Stage(this.nglDiv);
      await awaitNgl(this.ngl, callLogPrefix); // await NGL is ready
    } else
      this.ngl.removeAllComponents();

    this.nglDiv.style.width = `${Math.floor(task.props.width) / dpr}px`;
    this.nglDiv.style.height = `${Math.floor(task.props.height) / dpr}px`;
    this.ngl.viewer.setSize(task.props.width / dpr, task.props.height / dpr);
    const canvas = this.ngl.viewer.renderer.domElement;
    canvas.width = Math.floor(task.props.width);
    canvas.height = Math.floor(task.props.height);

    await this.ngl.loadFile(stringBlob, {ext: 'pdb', defaultRepresentation: true});
    await delay(10); /* Sometimes compList[0] is undefined, without any other error */
    const comp = this.ngl.compList[0];
    if (!comp) throw new Error(`${callLogPrefix}, ` + 'no component added');
    comp.addRepresentation('cartoon');
    comp.autoView();

    const emptyCanvasHash = DG.StringUtils.hashCode(canvas.toDataURL());
    const nglRenderedBinding = this.ngl.viewer.signals.rendered.add(renderHandler, this);
    this.ngl.viewer.requestRender();
    return [nglRenderedBinding, emptyCanvasHash];
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

  private nglOnRendered(key: keyof any | undefined, task: NglGlTask, emptyCanvasHash: number): boolean {
    _package.logger.debug('NglGlDocService.onNglRendered(), ' + `key: ${key?.toString()}`);
    const canvas = this.ngl!.viewer.renderer.domElement;
    const canvasHash = DG.StringUtils.hashCode(canvas.toDataURL());
    if (canvasHash == emptyCanvasHash) { // render is not ready yet
      _package.logger.debug(`NglGlDocService.onNglRendered(), empty canvas`);
      return false;
    }
    _package.logger.debug(`NglGlDocService.onNglRendered(), ` + `key: ${key?.toString()}, ` +
      `this.emptyCanvasHash = ${emptyCanvasHash}, canvasHash = ${canvasHash}`);

    task.onAfterRender(canvas);
    return true;
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
