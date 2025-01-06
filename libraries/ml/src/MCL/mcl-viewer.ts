/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MCLSerializableOptions} from './types';
import {markovCluster, SCLinesRenderer} from './clustering-view';

export type MCLViewerProps = {
    MCLProps: MCLSerializableOptions;
    lines: string;
    scProps: Partial<DG.IScatterPlotSettings>;
}


// depending if the dataframe has a data sync enabled or not, we might want to initialize the viewer in a different way
export class MCLViewer extends DG.JsViewer {
    public sc?: DG.ScatterPlotViewer;
    mclProps: string;
    lines: string; // this will be encoded string in unicode (16bit) limiting number of dataframe rows to 65535. every odd row will be a start point, every even row will be an end point
    scProps: string;
    initPromise = Promise.resolve();
    private initialized = false;
    private reseolver: () => void = () => {};
    private linesRenderer?: SCLinesRenderer;
    scratchCallTimer: any | null = null;
    constructor() {
      super();
      this.mclProps = this.string('mclProps', null, {includeInLayout: false, userEditable: false, nullable: true});
      this.lines = this.string('lines', null, {includeInLayout: false, userEditable: false, nullable: true});
      this.scProps = this.string('scProps', null, {includeInLayout: false, userEditable: false, nullable: true});
      this.initPromise = new Promise((resolve) => {
        this.reseolver = resolve;
      });
    }

    onFrameAttached(dataFrame: DG.DataFrame): void {
      if (dataFrame.rowCount > 65535)
        throw new Error('MCL viewer supports only dataframes with less than 65535 rows');
      this.sc = dataFrame.plot.scatter({
        showXAxis: false,
        showYAxis: false,
        showXSelector: false,
        showYSelector: false,
        title: 'MCL',
      });
      this.root.appendChild(this.sc.root);
      this.subs.push(grok.events.onCurrentObjectChanged.subscribe((_) => {
        if (this.sc && grok.shell.o === this) {
          setTimeout(() => {
            grok.shell.o = this.sc;
          });
        }
      }));

      this.subs.push(DG.debounce(this.sc.onPropertyValueChanged, 1000).subscribe((_) => {
        if (this.sc)
            this.getProperty('scProps')!.set(this, JSON.stringify(Object.assign({}, this.sc.props)));
      }));
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
        this.scratchCallTimer && clearTimeout(this.scratchCallTimer);
        this.scratchCallTimer = setTimeout(() => { this.initFromScratch(); }, 300);
      }
    }

    async initFromScratch() {
      if (!this.mclProps || !this.sc || !this.dataFrame || this.initialized)
        return;
      if (this.lines) {
        // if lines are already provided, no need to init from scratch
        this.decodeLines();
        return;
      }
      const options: MCLSerializableOptions = JSON.parse(this.mclProps);

      const cols = options.cols.map((colName) => this.dataFrame.columns.byName(colName));
      const preprocessingFuncs = options.preprocessingFuncs.map((funcName) => funcName ? DG.Func.byName(funcName) : null);

      const res = await markovCluster(this.dataFrame, cols, options.metrics, options.weights,
        options.aggregationMethod, preprocessingFuncs, options.preprocessingFuncArgs, options.threshold,
        options.maxIterations, options.useWebGPU, options.inflate, options.minClusterSize, this.sc);
      if (!res)
        return;
      // if dataframe has datasync enabled, we should not save the lines, as they will be saved in the data sync
      if (this.dataFrame.getTag('.script')) {
        this.linesRenderer?.destroy();
        this.linesRenderer = new SCLinesRenderer(this.sc!, res.i, res.j, 6, 0.75, '128,128,128');
        this.initialized = true;
        this.reseolver();
        return;
      }
      this.encodeLines(res.i, res.j);
    }

    decodeLines() {
      if (!this.lines)
        return;
      const len = this.lines.length;
      if (len % 2 !== 0)
        throw new Error('Invalid lines string');
      const is = new Array(len / 2).fill(null).map((_, i) => this.lines.charCodeAt(2 * i));
      const js = new Array(len / 2).fill(null).map((_, i) => this.lines.charCodeAt(2 * i + 1));
      this.linesRenderer?.destroy();
      this.linesRenderer = new SCLinesRenderer(this.sc!, is, js, 6, 0.75, '128,128,128');
      this.initialized = true;
      this.reseolver();
    }

    encodeLines(is: ArrayLike<number>, js: ArrayLike<number>) {
      const result = new Array(is.length).fill(null).map((_, i) => `${String.fromCharCode(is[i])}${String.fromCharCode(js[i])}`).join('');
        this.getProperty('lines')!.set(this, result);
    }
}
