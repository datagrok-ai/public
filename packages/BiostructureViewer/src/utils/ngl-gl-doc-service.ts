import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as ngl from 'NGL';

import {SignalBinding} from 'signals';
import {delay} from '@datagrok-libraries/utils/src/test';
import {NglGlServiceBase, NglGlTask} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {BiostructureData} from '@datagrok-libraries/bio/src/pdb/types';

import {awaitNgl} from '../viewers/ngl-viewer-utils';

import {_package} from '../package';

// const TASK_TIMEOUT: number = 2000;
const NGL_ERROR_LIMIT: number = 3;
const NGL_TRY_LIMIT: number = 3;

export class NglGlDocService extends NglGlServiceBase {
  private readonly logger: ILogger;

  private nglDiv: HTMLDivElement;

  private ngl: ngl.Stage | null = null;
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

    this.hostDiv = ui.box();
    // this.hostDiv.style.display = 'none'; // Disables drawing at all
    this.hostDiv.style.position = 'absolute';
    this.hostDiv.style.left = '0px';
    this.hostDiv.style.right = '0px';
    this.hostDiv.style.width = '0px';
    this.hostDiv.style.height = '0px';
    this.hostDiv.style.visibility = 'hidden';
    document.body.appendChild(this.hostDiv);

    this._queue = [];
    this._queueDict = {};

    this.logger = _package.logger;
    // window.setInterval(() => { this._sweepQueue(); }, 200);
  }

  public get errorCount(): number { return this.nglErrorCount; }

  public async reset(): Promise<void> {
    if (!this.ngl) return;
    this.ngl.dispose();
    this.ngl = null;
    $(this.nglDiv).empty();
    this.nglDiv.remove();
    this.nglErrorCount = 0;
  }

  render(task: NglGlTask, key?: keyof any, tryCount: number = 0): void {
    this.logger.debug('NglGlDocService.render() start ' + `key: ${key?.toString()}`);

    if (key !== undefined) {
      if (key in this._queueDict) {
        // remove outdated task from the queue
        const oldTaskI = this._queue.findIndex((item) => item.key === key);
        this._queue.splice(oldTaskI, 1);
      }
      this.logger.debug(`NglGlDocService.render() _queueDict[ key: ${key?.toString()} ] = <task>`);
      this._queueDict[key] = task;
    }

    this.logger.debug('NglGlDocService.render() _queue.push(), ' + `key: ${key?.toString()}`);
    this._queue.push({key, task, tryCount, dt: window.performance.now()});

    if (!this._busy) {
      this._busy = true;

      // TODO: Use requestAnimationFrame()
      this.logger.debug('NglGlDocService.render(), window.setTimeout() -> this._processQueue() ');
      window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
    }
  }

  private _processQueue(): void {
    if (this._queue.length === 0) return;
    this.logger.debug(`NglGlDocService._processQueue(), ` +
      `queue: ${JSON.stringify(this._queue.map((t) => t.key))}`);

    let nglRenderBinding: SignalBinding<any> | null = null;

    const finallyProcessQueue = (caller: string) => {
      const logPrefix = `NglGlDocService._processQueue().finallyProcessQueue( <- ${caller} )`;
      if (nglRenderBinding !== null) nglRenderBinding.detach();
      if (this._queue.length > 0) {
        // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
        this.logger.debug(`${logPrefix}, ` + 'window.setTimeout() -> this._processQueue() ');
        window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
      } else {
        // release flag allowing _processQueue on add queue item
        this.logger.debug(`${logPrefix}, ` + 'this._busy = false');
        this._busy = false;
      }
    };

    const queueItem = this._queue.shift();
    if (!queueItem) {
      this.logger.error('NglGlDocService._processQueue() queueItem = undefined ');
      finallyProcessQueue('no queue item');
      return;
    } // in case of empty queue

    const {key, task, tryCount, dt} = queueItem;
    if (key !== undefined) {
      this.logger.debug('NglGlDocService._processQueue() ' + `key: ${key.toString()}`);
      delete this._queueDict[key];
    }
    if (tryCount > NGL_TRY_LIMIT) {
      this.logger.warning('NglGlDocService._processQueue(), skip task, ' +
        ` key: ${key?.toString()}, tryCount: ${tryCount}`);
      finallyProcessQueue('try count limit');
      return;
    }

    let emptyCanvasHash: number;
    let handled: boolean = false;

    const timeoutHandle = window.setTimeout(() => {
      if (!handled) {
        this.logger.warning('NglGlDocService._processQueue().timeoutHandle() not handled, ' +
          `key: ${key?.toString}`);
        this.nglErrorCount += 1;
        this.render(task, key, tryCount + 1); // return task to the queue
        finallyProcessQueue('timeout');
      }
    }, 1000);

    const renderHandler = () => {
      this.logger.debug('NglGlDocService._processQueue().handler(), ' + `key: ${key?.toString()}`);
      if (this.nglOnRendered(key, task, emptyCanvasHash)) {
        handled = true;
        window.clearTimeout(timeoutHandle);
        finallyProcessQueue('render');
      }
    };

    this.nglRender(key, task, renderHandler).then((res) => {
      //[nglRenderBinding, emptyCanvasHash] = res;
      nglRenderBinding = res[0];
      emptyCanvasHash = res[1];
      const triggerRender = res[2];
      if (emptyCanvasHash === undefined)
        console.warn('NglGlDocService._processQueue() this.nglRender.then() emptyCanvasHAsh undefined');
      triggerRender();
    })
      .catch((err: any) => {
        // Not waiting timeout on error
        const [errMsg, errStack] = errInfo(err);
        this.logger.error(errMsg, undefined, errStack);
        window.clearTimeout(timeoutHandle);
        this.nglErrorCount += 1;
        this.render(task, key, tryCount + 1); // return task to the queue
        finallyProcessQueue('nglRender error');
      });
  }

  private async nglRender(
    key: any, task: NglGlTask, renderHandler: () => void,
  ): Promise<[SignalBinding<any>, number, () => void]> {
    const logPrefix = 'NglGlDocService.nglRender()';
    this.logger.debug(`${logPrefix}, ` + `key: ${key?.toString()}`);
    const dpr = window.devicePixelRatio;

    // TODO: Convert string to Blob once converting PDB string column to Blob
    const stringBlob = new Blob([task.props.pdb], {type: 'text/plain'});

    // TODO: Use canvas size switching 0/1 px to required

    // if (key === 1) throw new Error('NglGlDocService: Test error');

    if (this.nglErrorCount > NGL_ERROR_LIMIT) {
      this.logger.warning(`${logPrefix}, ` + 'recreate ngl.Stage, ' +
        `nglTimeoutCount = ${this.nglErrorCount} > ${NGL_ERROR_LIMIT}`);
      await this.reset();
    }
    if (!this.ngl) {
      // re-create ngl div and stage
      this.hostDiv.appendChild(this.nglDiv = ui.div([], 'd4-ngl-viewer'));
      this.ngl = new ngl.Stage(this.nglDiv);
      await awaitNgl(this.ngl, logPrefix); // await NGL is ready
    } else {
      //
      this.ngl.removeAllComponents();
    }

    const lST = window.performance.now();
    this.nglDiv.style.width = `${Math.floor(task.props.width) / dpr}px`;
    this.nglDiv.style.height = `${Math.floor(task.props.height) / dpr}px`;
    this.ngl.viewer.setSize(task.props.width / dpr, task.props.height / dpr);
    const canvas = this.ngl.viewer.renderer.domElement;
    canvas.width = Math.floor(task.props.width);
    canvas.height = Math.floor(task.props.height);

    await this.ngl.loadFile(stringBlob, {ext: 'pdb', defaultRepresentation: true});
    const comp = this.ngl.compList[0];
    comp.addRepresentation('cartoon', undefined);
    comp.autoView();
    // await delay(200); /* Sometimes compList[0] is undefined, without any other error */
    if (!comp) throw new Error('no component added');
    const lET = window.performance.now();

    const echST = window.performance.now();
    const emptyCanvasHash = DG.StringUtils.hashCode(canvas.toDataURL());
    if (emptyCanvasHash === undefined)
      this.logger.warning(`${logPrefix}, emptyCanvasHash undefined at calc`);
    const echET = window.performance.now();

    const rST = window.performance.now();
    const renderHandlerInt = (...params: any): void => {
      const rET: number = window.performance.now();
      this.logger.warning('NglGlDocService.nglRender().renderHandlerInt(), ' +
        `key: ${key?.toString()}, ` +
        `load: ${lET - lST} ms,\n    ` + `emptyCanvasHashET: ${echET - echST} ms, ` + `renderET: ${rET - rST} ms, ` +
        `emptyCanvasHash: ${emptyCanvasHash}, ` + `params: ${JSON.stringify(params)}, `);
      renderHandler();
    };
    const nglRenderedBinding = this.ngl.viewer.signals.rendered.add(renderHandlerInt, this);
    const trigger = () => {
      // this.ngl.viewer.requestRender();
      this.ngl!.viewer.render(false);
    };
    if (emptyCanvasHash === undefined)
      console.warn('NglGlDocService.nglRender() emptyCanvasHash undefined at the end');
    // trigger rendering after returning nglRenderBinding and emptyCanvasHash
    return [nglRenderedBinding, emptyCanvasHash, trigger];
  }

  // /** Sweep queue for stalled tasks */
  // private _sweepQueue(): void {
  //   const logPrefix: string = `NglGlDocService._sweepQueue()`;
  //   const nowDt: number = window.performance.now();
  //   let swept: boolean = false;
  //
  //   for (let qI = this._queue.length - 1; qI >= 0; qI--) {
  //     const {key, task: _task, dt} = this._queue[qI];
  //     if ((nowDt - dt) > TASK_TIMEOUT) {
  //       // stalled task
  //       this.logger.warning(`${logPrefix}, remove task key = ${key?.toString()}`);
  //       this._queue.splice(qI);
  //       delete this._queueDict[key!];
  //       swept = true;
  //     }
  //   }
  //
  //   this._busy = this._queue.length > 0;
  //   if (!this._busy && swept) {
  //     // Some tasks had been swept
  //   }
  // }

  private nglOnRendered(key: keyof any | undefined, task: NglGlTask, emptyCanvasHash: number): boolean {
    const logPrefix = `NglGlDocService.nglOnRendered( key = ${key?.toString()} )`;
    if (emptyCanvasHash === undefined)
      this.logger.warning(`${logPrefix}, emptyCanvasHash undefined`);
    this.logger.debug(`${logPrefix}, start`);
    const canvas = this.ngl!.viewer.renderer.domElement;
    const canvasHash = DG.StringUtils.hashCode(canvas.toDataURL());
    if (canvasHash == emptyCanvasHash) { // render is not ready yet
      this.logger.debug(`NglGlDocService.onNglRendered(), empty canvas`);
      return false;
    }
    this.logger.debug(`${logPrefix}, end, ` + `emptyCanvasHash = ${emptyCanvasHash}, canvasHash = ${canvasHash}`);

    task.onAfterRender(canvas);
    return true;
  }

  public renderOnGridCell(
    gCtx: CanvasRenderingContext2D, bd: DG.Rect, gCell: DG.GridCell, canvas: CanvasImageSource
  ): void {
    const callLog = `renderOnGridCell( gRow = ${gCell.gridRow} )`;
    const logPrefix = `NglGlDocService.${callLog}`;
    const dpr = window.devicePixelRatio;
    gCtx.save();
    gCtx.resetTransform();
    try {
      this.logger.debug(`${logPrefix}, start`);

      const vertScrollMin: number = Math.floor(gCell.grid.vertScroll.min);
      // Correction for vert scrolling happened between task and render, calculate bd.y directly
      bd.y = ((gCell.gridRow - vertScrollMin) * gCell.grid.props.rowHeight +
        gCell.grid.colHeaderHeight) * dpr;

      // Correction for horz scrolling happened between task and render, calculate bd.x directly
      let left: number = 0;
      for (let colI = 0; colI < gCell.gridColumn.idx; colI++) {
        const col: DG.GridColumn = gCell.grid.columns.byIndex(colI)!;
        left += col.visible ? col.width : 0;
      }
      bd.x = (left - gCell.grid.horzScroll.min) * dpr;

      gCtx.beginPath();
      gCtx.rect(bd.x + 1, bd.y + 1, bd.width * dpr - 2, bd.height * dpr - 2);
      gCtx.clip();

      // gCtx.fillStyle = '#E0E0FF';
      // gCtx.fillRect(bd.x + 1, bd.y + 1, bd.width - 2, bd.height - 2);

      const cw = 'width' in canvas && (
        canvas.width instanceof SVGAnimatedLength ? canvas.width.baseVal.value : canvas.width as number);
      const ch = 'height' in canvas && (
        canvas.height instanceof SVGAnimatedLength ? canvas.height.baseVal.value : canvas.height as number);
      if (!cw || !ch) throw new Error('NglGlDocService.renderOnGridCell() canvas size is not available');
      gCtx.transform(1 /* bd.width / cw */, 0, 0, 1 /* bd.height / ch */, bd.x, bd.y);

      // gCtx.fillStyle = '#FFF0F020';
      // gCtx.fillRect(0 + 1, 0 + 1, cw - 2, ch - 2);
      // gCtx.textBaseline = 'top';
      // gCtx.fillStyle = 'green';
      // gCtx.font = 'bold 48px monospace';
      // gCtx.fillText(`t# ${gCell.gridRow}, g# ${gCell.tableRowIndex}\n${tag}`, 3, 3);
      gCtx.drawImage(canvas, 0 + 1, 0 + 1);
    } finally {
      gCtx.restore();
    }
  }
}
