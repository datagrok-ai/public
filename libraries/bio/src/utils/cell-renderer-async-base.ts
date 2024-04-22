import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';
import {Observable, Subject, Unsubscribable} from 'rxjs';

import {testEvent} from '@datagrok-libraries/utils/src/test';

import {ILogger} from './logger';
import {errInfo} from './err-info';

export class PropsBase {
  public constructor(
    public readonly backColor: number,
    public readonly width: number,
    public readonly height: number,
  ) {}
}

export abstract class RenderServiceBase<TProps extends PropsBase> {
  protected readonly _queue: {
    key?: keyof any,
    task: RenderTask<TProps>,
    tryCount: number,
    dt: number
  }[];
  protected readonly _queueDict: {
    [key: keyof any]: RenderTask<TProps>
  };
  /** The flag allows {@link _processQueue}() on add item to the queue with {@link render}() */
  private _busy: boolean = false;

  constructor(
    protected readonly logger: ILogger,
    protected readonly retryLimit: number = 3,
    protected _errorCount: number = 0,
    protected readonly errorLimit: number = 3,
  ) {
    this._queue = [];
    this._queueDict = {};
  }

  private static viewerCounter: number = -1;
  private readonly viewerId: number = RenderServiceBase.viewerCounter;

  protected toLog(): string {
    return `${this.constructor.name}<${this.viewerId}>`;
  }

  public get errorCount(): number { return this._errorCount; }

  /** Disposes MolstarViewer to free WebGL context */
  abstract reset(): Promise<void>;

  /** Queues render tasks
   * @param task Task to render
   * @param key  Specify to skip previously queued tasks with the same key
   */
  render(task: RenderTask<TProps>, key?: keyof any, tryCount: number = 0): void {
    const logPrefix = `${this.toLog()}.render()`;
    this.logger.debug(`${logPrefix}, start ` + `key: ${key?.toString()}`);

    if (key !== undefined) {
      if (key in this._queueDict) {
        // remove outdated task from the queue
        const oldTaskI = this._queue.findIndex((item) => item.key === key);
        this._queue.splice(oldTaskI, 1);
      }
      this.logger.debug(`${logPrefix}, _queueDict[ key: ${key?.toString()} ] = <task>`);
      this._queueDict[key] = task;
    }

    this.logger.debug(`${logPrefix}, _queue.push(), ` + `key: ${key?.toString()}`);
    this._queue.push({key, task, tryCount, dt: window.performance.now()});

    if (!this._busy) {
      this._busy = true;

      // TODO: Use requestAnimationFrame()
      this.logger.debug(`${logPrefix}, window.setTimeout() -> this._processQueue() `);
      window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
    }
  }

  /** Default implementation of rendering tree on grid cell
   * @param gCtx    Context to draw on grid
   * @param bd      Bound rect to clip drawing on task moment
   * @param gCell   Grid cell to draw
   * @param canvas  Image of the tree rendered
   */
  public renderOnGridCell(
    gCtx: CanvasRenderingContext2D, bd: DG.Rect, gCell: DG.GridCell, canvas: CanvasImageSource
  ): void {
    const callLog = `renderOnGridCell( gRow = ${gCell.gridRow} )`;
    const logPrefix = `${this.toLog()}.${callLog}`;
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
      if (!cw || !ch) throw new Error(`${logPrefix}, canvas size is not available`);
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

  private _processQueue(): void {
    const logPrefix = `${this.toLog()}._processQueue()`;
    if (this._queue.length === 0) return;
    this.logger.debug(`${logPrefix}, ` +
      `queue: ${JSON.stringify(this._queue.map((t) => t.key))}`);
    let renderBinding: Unsubscribable | null = null;
    const finallyProcessQueue = (caller: string) => {
      const logPrefixR = `${logPrefix}.finallyProcessQueue( <- ${caller} )`;
      if (renderBinding !== null) renderBinding.unsubscribe();
      if (this._queue.length > 0) {
        // Schedule processQueue the next item only afterRender has asynchronously completed for the previous one
        this.logger.debug(`${logPrefixR}, ` + 'window.setTimeout() -> this._processQueue() ');
        window.setTimeout(() => { this._processQueue(); }, 0 /* next event cycle */);
      } else {
        // release flag allowing _processQueue on add queue item
        this.logger.debug(`${logPrefixR}, ` + 'this._busy = false');
        this._busy = false;
      }
    };

    const queueItem = this._queue.shift();
    if (!queueItem) {
      this.logger.error(`${logPrefix}, queueItem = undefined `);
      finallyProcessQueue('no queue item');
      return;
    } // in case of empty queue

    const {key, task, tryCount, dt} = queueItem;
    if (key !== undefined) {
      this.logger.debug(`${logPrefix}, ` + `key: ${key.toString()}`);
      delete this._queueDict[key];
    }
    if (tryCount > this.retryLimit) {
      this.logger.warning(`${logPrefix}, skip task, ` +
        ` key: ${key?.toString()}, tryCount: ${tryCount}`);
      finallyProcessQueue('try count limit');
      return;
    }

    let emptyCanvasHash: number;
    let handled: boolean = false;

    const timeoutHandle = window.setTimeout(() => {
      if (!handled) {
        this.logger.warning(`${logPrefix}.timeoutHandle(), not handled, ` +
          `key: ${key?.toString}`);
        this._errorCount += 1;
        this.render(task, key, tryCount + 1); // return task to the queue
        finallyProcessQueue('timeout');
      }
    }, 1000);

    const renderHandler = () => {
      this.logger.debug(`${logPrefix}.renderHandler(), ` + `key: ${key?.toString()}`);
      if (this.onRendered(key, task, emptyCanvasHash)) {
        handled = true;
        window.clearTimeout(timeoutHandle);
        finallyProcessQueue('render');
      }
    };

    this.requestRender(key, task, renderHandler).then((res) => {
      //[renderBinding, emptyCanvasHash] = res;
      renderBinding = res[0];
      emptyCanvasHash = res[1];
      const triggerRender = res[2];
      if (emptyCanvasHash === undefined)
        console.warn(`${logPrefix}, this.requestRender.then() emptyCanvasHash undefined`);
      triggerRender();
    })
      .catch((err: any) => {
        // Not waiting timeout on error
        const [errMsg, errStack] = errInfo(err);
        this.logger.error(errMsg, undefined, errStack);
        window.clearTimeout(timeoutHandle);
        this._errorCount += 1;
        this.render(task, key, tryCount + 1); // return task to the queue
        finallyProcessQueue(`${logPrefix} this.requestRender.catch()`);
      });
  }

  // /** Sweep queue for stalled tasks */
  // private _sweepQueue(): void {
  //   const logPrefix: string = `${this.toLog()}._sweepQueue()`;
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

  // -- Actual implementations --

  protected abstract requestRender(
    key: keyof any | undefined, task: RenderTask<TProps>, renderHandler: () => void
  ): Promise<[Unsubscribable, number, () => void]>;

  protected abstract onRendered(
    key: keyof any | undefined, task: RenderTask<TProps>, emptyCanvasHash: number
  ): boolean;
}

export class RenderTask<TProps extends PropsBase> {
  public constructor(
    public readonly name: string,
    public readonly props: TProps,
    public readonly onAfterRender: (canvas: HTMLCanvasElement) => void
  ) {}
}

export abstract class CellRendererBackAsyncBase<TProps extends PropsBase> {
  protected readonly taskQueueMap = new Map<number, RenderTask<TProps>>;
  protected readonly imageCache = new Map<number, HTMLImageElement>();

  protected constructor(
    protected readonly gridCol: DG.GridColumn,
    protected readonly logger: ILogger,
  ) {}

  protected abstract getRenderService(): RenderServiceBase<TProps>;

  protected abstract getRenderTaskProps(gridCell: DG.GridCell, dpr: number): TProps;

  public render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ): void {
    const r = window.devicePixelRatio;
    const service = this.getRenderService();

    if (gridCell.tableRowIndex == null || gridCell.tableColumn == null) return;

    const rowIdx: number = gridCell.tableRowIndex;
    this.logger.debug('PdbRenderer.render() start ' + `rowIdx=${rowIdx}`);

    g.save();
    try {
      const cacheDisabled: boolean = false;
      const image = this.imageCache.has(rowIdx) ? this.imageCache.get(rowIdx) : null;

      if (cacheDisabled || !image ||
        Math.abs(image.width / r - gridCell.gridColumn.width) > 0.5 ||
        Math.abs(image.height / r - gridCell.grid.props.rowHeight) > 0.5
      ) {
        // draw image

        // mark rowIdx as in drawing state
        //imageCache[rowIdx] = null;
        //gridCell.tableColumn.temp.set(PDB_RENDERER_IMAGE_CACHE_KEY, imageCache);

        const task: RenderTask<TProps> = new RenderTask<TProps>(
          gridCell.cell.rowIndex.toString(),
          this.getRenderTaskProps(gridCell, r),
          (canvas: HTMLCanvasElement) => {
            g.save();
            try {
              this.logger.debug('PdbRenderer.render() onAfterRender() ' + `rowIdx = ${rowIdx}`);
              service.renderOnGridCell(g, new DG.Rect(x, y, w, h), gridCell, canvas);

              const imageStr: string = canvas.toDataURL();
              base64ToImg(imageStr).then((image) => {
                this.imageCache.set(rowIdx, image);
              });
            } finally {
              g.restore();
              this.taskQueueMap.delete(rowIdx);
              if (this.taskQueueMap.size === 0)
                this._onRendered.next();
            }
          });

        this.taskQueueMap.set(rowIdx, task);
        service.render(task, rowIdx);
      } else {
        this.logger.debug('PdbRenderer.render(), ' + `from imageCache[${rowIdx}]`);
        service.renderOnGridCell(g, new DG.Rect(x, y, w, h), gridCell, image);
      }
    } catch (err: any) {
      const errMsg: string = err instanceof Error ? err.message : err.toString();
      this.logger.error(`BsV:PdbGridCellRenderer.render() no rethrown error: ${errMsg}`, undefined,
        err instanceof Error ? err.stack : undefined);
      //throw err; // Do not throw to prevent disabling event handler
    } finally {
      g.restore();
    }
  }

  // -- IRenderer --

  private _onRendered: Subject<void> = new Subject<void>();

  get onRendered(): Observable<void> { return this._onRendered; }

  invalidate(caller?: string): void {
    this.gridCol.grid.invalidate();
  }

  async awaitRendered(
    timeout: number = 10000, reason: string = `await rendered ${timeout} timeout`
  ): Promise<void> {
    await testEvent(this._onRendered, () => {}, () => { this.invalidate(); },
      timeout, reason);
  }
}

async function base64ToImg(base64: string): Promise<HTMLImageElement> {
  return new Promise<HTMLImageElement>((resolve, reject) => {
    try {
      const img = new Image();
      img.onload = () => {
        resolve(img);
      };
      img.src = base64;
    } catch (err) {
      reject(err);
    }
  });
}
