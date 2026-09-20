/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MCLComputationInfo, MCLSerializableOptions} from './types';
import {markovCluster, SCLinesRenderer} from './clustering-view';
import {Observable, Subject, timer} from 'rxjs';
import {debounce, tap} from 'rxjs/operators';

export type MCLViewerProps = {
    MCLProps: MCLSerializableOptions;
    lines: string;
    scProps: Partial<DG.IScatterPlotSettings>;
}

export const MAX_MCL_SAVABLE_ROWS = 65534;

// Native handles identify the same table across JS wrappers and replacement MCL viewers.
// Weak keys let a closed table and its computation history be collected together.
const completedComputations = new WeakMap<object, number>();

// depending if the dataframe has a data sync enabled or not, we might want to initialize the viewer in a different way
export class MCLViewer extends DG.JsViewer {
    public sc?: DG.ScatterPlotViewer;
    mclProps: string;
    lines: string; // this will be encoded string in unicode (16bit) limiting number of dataframe rows to 65535. every odd row will be a start point, every even row will be an end point
    scProps: string;
    initPromise = Promise.resolve();
    private initialized = false;
    private resolveInitialization: () => void = () => {};
    private rejectInitialization: (error: unknown) => void = () => {};
    private linesRenderer?: SCLinesRenderer;
    private initializing = false;
    private savingProperties = false;
    private initializationError: string | null = null;
    private initializationGeneration = 0;
    private completedComputation?: MCLComputationInfo;
    private rendered = new Subject<void>();
    scratchCallTimer: ReturnType<typeof setTimeout> | null = null;
    constructor() {
      super();
      this.mclProps = this.string('mclProps', null, {includeInLayout: false, userEditable: false, nullable: true});
      this.lines = this.string('lines', null, {includeInLayout: false, userEditable: false, nullable: true});
      this.scProps = this.string('scProps', null, {includeInLayout: false, userEditable: false, nullable: true});
      this.initPromise = new Promise((resolve, reject) => {
        this.resolveInitialization = resolve;
        this.rejectInitialization = reject;
      });
    }

    get onRendered(): Observable<void> { return this.rendered; }

    get isRenderPending(): boolean {
      return this.scratchCallTimer !== null || this.initializing || this.savingProperties ||
        (this.sc?.isRenderPending ?? false) || (this.linesRenderer?.isRenderPending ?? false);
    }

    get immediateRendering(): boolean { return super.immediateRendering; }

    set immediateRendering(value: boolean) {
      super.immediateRendering = value;
      if (this.sc)
        this.sc.immediateRendering = value;
    }

    getWidgetStatus(): DG.IWidgetStatus {
      const inner = this.sc?.getWidgetStatus();
      return {
        parts: {root: this.root, ...inner?.parts}, hitAreas: inner?.hitAreas ?? {},
        values: {...inner?.values, initialized: this.initialized,
          'completed computations': this.dataFrame ? completedComputations.get(this.dataFrame.dart) ?? 0 : 0,
          ...(this.completedComputation ? {'completed threshold': this.completedComputation.threshold,
            'completed inflation': this.completedComputation.inflation} : {}),
          ...(this.linesRenderer ? {connections: this.linesRenderer.from.length} : {})},
        shortcuts: inner?.shortcuts ?? {}, events: [], description: null,
        error: this.initializationError ?? inner?.error ?? null,
      };
    }

    onFrameAttached(dataFrame: DG.DataFrame): void {
      this.initializationGeneration++;
      this.completedComputation = undefined;
      // if (dataFrame.rowCount > 65535)
      //   throw new Error('MCL viewer supports only dataframes with less than 65535 rows');
      this.sc = dataFrame.plot.scatter({
        showXAxis: false,
        showYAxis: false,
        showXSelector: false,
        showYSelector: false,
        title: 'MCL',
        markerType: DG.MARKER_TYPE.CIRCLE
      });
      this.sc.immediateRendering = this.immediateRendering;
      this.root.appendChild(this.sc.root);
      this.subs.push(this.sc.onAfterDrawScene.subscribe(() => this.rendered.next()));
      this.subs.push(grok.events.onCurrentObjectChanged.subscribe((_) => {
        if (this.sc && grok.shell.o === this) {
          setTimeout(() => {
            grok.shell.o = this.sc;
          });
        }
      }));

      this.subs.push(this.sc.onPropertyValueChanged.pipe(
        tap(() => this.savingProperties = true),
        debounce(() => timer(this.immediateRendering ? 0 : 1000)),
      ).subscribe((_) => {
        this.savingProperties = false;
        if (this.sc)
            this.getProperty('scProps')!.set(this, JSON.stringify(Object.assign({}, this.sc.props)));
      }));
      if (this.mclProps)
        this.scheduleInitialization();
    }

    setScProps() {
      if (!this.sc)
        return;
      const curPropsStr = this.scProps;
      if (!curPropsStr)
        return;
      const scProps = this.sc.props;
      const scPropsStr = JSON.stringify(Object.assign({}, scProps));
      if (curPropsStr === scPropsStr)
        return;
      //its better to chech the changed property and set it that way. otherwise sc internally will set every property one by one
      const newProps = JSON.parse(curPropsStr);

      Object.entries(newProps).forEach(([key, value]: [string, any]) => {
        if (scProps.hasProperty(key) && scProps[key as keyof typeof scProps] !== value)
          (this.sc!.props as any)[key] = value;
      });
    }

    onPropertyChanged(property: DG.Property | null): void {
      if (property === null)
        return;

      if (property.name === 'scProps') {
        this.setScProps();
      } else if (property.name === 'lines') {
        this.decodeLines();
      } else if (property.name === 'mclProps') {
        this.scheduleInitialization();
      }
    }

    private scheduleInitialization(): void {
      if (this.scratchCallTimer !== null)
        clearTimeout(this.scratchCallTimer);
      this.scratchCallTimer = setTimeout(() => {
        this.scratchCallTimer = null;
        const generation = this.initializationGeneration;
        this.initFromScratch().catch((error) => {
          if (generation !== this.initializationGeneration || this.isDetached)
            return;
          this.initializationError = error instanceof Error ? error.message : String(error);
          if (this.sc)
            ui.setUpdateIndicator(this.sc.root, false);
          this.rejectInitialization(error);
          this.rendered.next();
        });
      }, this.immediateRendering ? 0 : 300);
    }

    public isDataFrameSavable(): boolean {
      return this.dataFrame != null && this.dataFrame.rowCount < MAX_MCL_SAVABLE_ROWS;
    }

    async initFromScratch() {
      if (!this.mclProps || !this.sc || !this.dataFrame || this.initialized || this.initializing || this.isDetached)
        return;
      const generation = this.initializationGeneration;
      const dataFrame = this.dataFrame;
      const scatter = this.sc;
      this.initializing = true;
      this.initializationError = null;
      try {
        if (this.lines !== null && this.lines !== undefined) {
          // if lines are already provided, no need to init from scratch
          this.decodeLines();
          return;
        }
        const options: MCLSerializableOptions = JSON.parse(this.mclProps);

        const cols = options.cols.map((colName) => dataFrame.columns.byName(colName));
        const preprocessingFuncs = options.preprocessingFuncs.map((funcName) => funcName ? DG.Func.byName(funcName) : null);

        const res = await markovCluster(dataFrame, cols, options.metrics, options.weights,
          options.aggregationMethod, preprocessingFuncs, options.preprocessingFuncArgs, options.threshold,
          options.maxIterations, options.useWebGPU, options.inflate, options.minClusterSize, scatter);
        if (generation !== this.initializationGeneration || this.isDetached ||
          this.dataFrame?.dart !== dataFrame.dart || this.sc !== scatter)
          return;
        if (!res)
          throw new Error('MCL clustering did not produce a result');
        completedComputations.set(dataFrame.dart, (completedComputations.get(dataFrame.dart) ?? 0) + 1);
        this.completedComputation = res.computation;
        // if dataframe has datasync enabled, we should not save the lines, as they will be saved in the data sync
        if (this.dataFrame.getTag('.script') || !this.isDataFrameSavable()) {
          this.linesRenderer?.destroy();
          this.linesRenderer = new SCLinesRenderer(this.sc!, res.i, res.j, 6, 0.75, '128,128,128');
          this.initialized = true;
          this.resolveInitialization();
          this.sc.invalidateCanvas();
          return;
        }
        this.encodeLines(res.i, res.j);
      } finally {
        if (generation === this.initializationGeneration)
          this.initializing = false;
      }
    }

    decodeLines() {
      if (this.lines === null || this.lines === undefined || !this.sc)
        return;
      const len = this.lines.length;
      if (len % 2 !== 0)
        throw new Error('Invalid lines string');
      const is = new Array(len / 2).fill(null).map((_, i) => this.lines.charCodeAt(2 * i));
      const js = new Array(len / 2).fill(null).map((_, i) => this.lines.charCodeAt(2 * i + 1));
      this.linesRenderer?.destroy();
      this.linesRenderer = new SCLinesRenderer(this.sc!, is, js, 6, 0.75, '128,128,128');
      this.initialized = true;
      this.resolveInitialization();
      this.sc.invalidateCanvas();
    }

    encodeLines(is: ArrayLike<number>, js: ArrayLike<number>) {
      const result = new Array(is.length).fill(null).map((_, i) => `${String.fromCharCode(is[i])}${String.fromCharCode(js[i])}`).join('');
      this.getProperty('lines')!.set(this, result);
    }

    detach(): void {
      this.initializationGeneration++;
      this.initializing = false;
      this.completedComputation = undefined;
      if (this.scratchCallTimer !== null)
        clearTimeout(this.scratchCallTimer);
      this.scratchCallTimer = null;
      this.linesRenderer?.destroy();
      this.resolveInitialization();
      this.rendered.complete();
      super.detach();
    }
}
