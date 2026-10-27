import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Subscription} from 'rxjs';

const CORRELATIONS_TITLE = 'Correlations';
const RAW_TITLE = 'Raw';
const PANEL_TITLE = 'Weights';
const FINAL_TITLE = 'Final';
const WEIGHTS_TITLE = 'Weights';

export type WeightStatus = 'kept' | 'pinned' | 'dropped';

// Display-only stub so a zero bar keeps its row on the chart; scores and tooltips use the real value.
const ZERO_BAR = 0.004;

const SIGN_COLORS = {positive: 0xFF1F77B4, negative: 0xFFFF7F0E};
const STATUS_COLORS: Record<WeightStatus, number> = {kept: 0xFF1F77B4, pinned: 0xFF2CA02C, dropped: 0xFF8A8F99};

export function weightsFrame(names: string[], raw: Float64Array, final: Float64Array,
  status: WeightStatus[]): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('property', names),
    DG.Column.fromFloat32Array('raw', Float32Array.from(raw, (r) => Number.isNaN(r) ? DG.FLOAT_NULL : r)),
    DG.Column.fromStrings('sign', Array.from(raw, (r) => r < 0 ? 'negative' : 'positive')),
    DG.Column.fromFloat32Array('weight', Float32Array.from(final, (w) => w === 0 ? ZERO_BAR : w)),
    DG.Column.fromStrings('status', status),
  ]);
  df.col('sign')!.meta.colors.setCategorical(SIGN_COLORS);
  df.col('status')!.meta.colors.setCategorical(STATUS_COLORS);
  return df;
}

/** [valueMin, valueMax] of a chart over the finite values of `a`, rounded outward to 0.1 with room for a bar that
 *  ends on a round value (the chart reads it as float32, a hair above the float64 bound). */
export function axisRange(a: Float64Array, withNegative: boolean): [number, number] {
  let lo = 0;
  let hi = 0.1;
  for (const v of a) {
    if (Number.isFinite(v)) {
      lo = Math.min(lo, v);
      hi = Math.max(hi, v);
    }
  }
  return [withNegative ? Math.floor(Math.min(lo, -0.1) * 10 - 1e-6) / 10 : 0, Math.ceil(hi * 10 + 1e-6) / 10];
}

function createCorrPlot(df: DG.DataFrame, names: string[]): DG.Viewer {
  return DG.Viewer.fromType(DG.VIEWER.CORR_PLOT, df, {
    xColumnNames: names,
    yColumnNames: names,
    correlationType: 'Pearson',
    title: CORRELATIONS_TITLE,
  });
}

// Re-applied on every data update: a bar chart given a new data frame picks its own value column.
function barChartOptions(value: 'raw' | 'weight', stack: 'sign' | 'status', title: string,
  [valueMin, valueMax]: [number, number]): object {
  return {
    splitColumnName: 'property',
    valueColumnName: value,
    valueAggrType: 'sum',
    stackColumnName: stack,
    orientation: 'Horizontal',
    barSortType: 'by value',
    barSortOrder: 'desc',
    valueMin, valueMax,
    showValueSelector: false,
    showCategorySelector: false,
    showStackSelector: false,
    title, showTitle: false,
    legendVisibility: 'Always', legendPosition: 'RightBottom',
  };
}

function showBarChart(chart: DG.Viewer | null, df: DG.DataFrame, options: object): DG.Viewer {
  if (!chart)
    return DG.Viewer.fromType(DG.VIEWER.BAR_CHART, df, options);
  chart.dataFrame = df;
  chart.setOptions(options);
  return chart;
}

const RAW_TOOLTIP = 'Weights straight from the correlations between the desirabilities of the properties, before ' +
  'any adjustment. A property that the others already cover gets less; a negative weight means it is redundant.';
const FINAL_TOOLTIP = 'Weights used in the score: negative raw weights set to 0, then priorities and pinned weights ' +
  'applied. They sum to 1; the dashed line marks the equal weight.';
const WEIGHTS_TOOLTIP = 'Weights used in the score, from the correlations between the desirabilities of the ' +
  'properties, with priorities and pinned weights applied. A property that the others already cover gets less. ' +
  'They sum to 1; the dashed line marks the equal weight.';

/** A chart under its own title, which, unlike the viewer's title, can carry a tooltip. */
function chartPane(chart: DG.Viewer, title: string, tooltip: string): HTMLElement {
  const header = ui.divText(title, 'chem-mpo-campo-chart-title');
  ui.tooltip.bind(header, tooltip);
  return ui.divV([header, ui.div([chart.root], 'chem-mpo-campo-slot')], 'chem-mpo-campo-chart');
}

/** Two panes with a draggable divider, sized by flex shares so they always fill the container. `ui.splitV` and
 *  `ui.splitH` fix pixel sizes when first inserted, which is before the dock has laid the panel out. */
function flexSplit(a: HTMLElement, b: HTMLElement, vertical: boolean): HTMLElement {
  a.style.flex = b.style.flex = '1 1 0';
  const divider = ui.div([], vertical ? 'chem-mpo-campo-divider-v' : 'chem-mpo-campo-divider-h');
  const root = ui.div([a, divider, b], vertical ? 'chem-mpo-campo-split-v' : 'chem-mpo-campo-split-h');
  divider.addEventListener('mousedown', (e) => {
    e.preventDefault();
    const box = root.getBoundingClientRect();
    const move = (m: MouseEvent) => {
      const share = vertical ? (m.clientY - box.top) / box.height : (m.clientX - box.left) / box.width;
      const s = Math.min(0.9, Math.max(0.1, share));
      a.style.flexGrow = `${s}`;
      b.style.flexGrow = `${1 - s}`;
    };
    const up = () => {
      document.removeEventListener('mousemove', move);
      document.removeEventListener('mouseup', up);
    };
    document.addEventListener('mousemove', move);
    document.addEventListener('mouseup', up);
  });
  return root;
}

/** The right half of the profile view: a detached table view of the dataset. In Correlation-aware mode one panel is
 *  docked as a tab over its grid: the Correlation plot above, the Computed and Applied charts side by side below. */
export class CampoViewers {
  readonly view: DG.TableView;
  private corr: DG.Viewer | null = null;
  private rawChart: DG.Viewer | null = null;
  private finalChart: DG.Viewer | null = null;
  private corrSlot = ui.div([], 'chem-mpo-campo-slot');
  private chartsSlot = ui.div([], 'chem-mpo-campo-slot');
  private panel = flexSplit(this.corrSlot, this.chartsSlot, true);
  private docked = false;
  private closedByUser = false;
  private added = false;
  // the latest show of each part requested before the view was added
  private pending: Map<string, () => void> = new Map();
  private subs: Subscription[] = [];

  constructor(df: DG.DataFrame) {
    this.view = DG.TableView.create(df, false);
    this.view.root.classList.add('chem-mpo-campo-view');
    // the dock fills the view, which a `ui-box` child rule would otherwise give scrollbars over a pixel of rounding
    this.view.root.style.setProperty('overflow', 'hidden', 'important');
    // a detached view has no dock manager until it is added, once its root is in the document
    this.subs.push(ui.onSizeChanged(this.view.root).subscribe(() => {
      if (this.added)
        this.resize();
      else if (this.view.root.isConnected) {
        this.view._onAdded();
        this.added = true;
        this.subs.push(this.view.dockManager.onClosed.subscribe((el: HTMLElement) => {
          if (el === this.panel && this.docked) {
            this.closedByUser = true;
            this.clear();
          }
        }));
        for (const show of this.pending.values())
          show();
        this.pending.clear();
      }
    }));
  }

  private whenAdded(key: string, show: () => void): void {
    if (this.added)
      show();
    else
      this.pending.set(key, show);
  }

  /** The mode was switched on, or RESET restored it: a panel the user closed comes back. */
  open(): void {
    this.closedByUser = false;
  }

  private dock(): void {
    if (this.docked)
      return;
    const dm = this.view.dockManager;
    dm.dock(this.panel, DG.DOCK_TYPE.FILL, dm.findNode(this.view.grid.root) ?? null, PANEL_TITLE);
    this.docked = true;
  }

  /** Lays the dock out to the view's size. The JS TableView lacks `_handleResize`, its Dart view is a DockView. */
  private resize(): void {
    DG.DockView.prototype._handleResize.call(this.view);
  }

  showCorr(frame: DG.DataFrame, names: string[]): void {
    this.whenAdded('corr', () => {
      if (this.closedByUser)
        return;
      if (this.corr) {
        this.corr.dataFrame = frame;
        this.corr.setOptions({xColumnNames: names, yColumnNames: names});
      } else {
        this.corr = createCorrPlot(frame, names);
        this.corrSlot.append(this.corr.root);
      }
      this.dock();
    });
  }

  showWeights(df: DG.DataFrame, raw: Float64Array, final: Float64Array, clipped: boolean): void {
    this.whenAdded('weights', () => {
      if (this.closedByUser)
        return;
      this.finalChart = showBarChart(this.finalChart, df,
        barChartOptions('weight', 'status', FINAL_TITLE, axisRange(final, false)));
      this.finalChart.meta.formulaLines.clear().addLine({
        formula: '${weight} = ' + (1 / final.length).toFixed(3),
        title: 'Equal', color: '#d62728', style: 'dashed',
      });
      const hadRaw = !!this.rawChart;
      this.rawChart = clipped ? showBarChart(this.rawChart, df,
        barChartOptions('raw', 'sign', RAW_TITLE, axisRange(raw, true))) : null;
      if (hadRaw !== !!this.rawChart || !this.chartsSlot.firstChild) {
        ui.empty(this.chartsSlot);
        this.chartsSlot.append(this.rawChart ? flexSplit(chartPane(this.rawChart, RAW_TITLE, RAW_TOOLTIP),
          chartPane(this.finalChart, FINAL_TITLE, FINAL_TOOLTIP), false) :
          chartPane(this.finalChart, WEIGHTS_TITLE, WEIGHTS_TOOLTIP));
      }
      this.dock();
    });
  }

  private clear(): void {
    ui.empty(this.corrSlot);
    ui.empty(this.chartsSlot);
    this.docked = false;
    this.corr = this.rawChart = this.finalChart = null;
  }

  /** The mode was switched off, or stage A failed: only the grid stays. */
  hide(): void {
    this.pending.clear();
    if (this.docked) {
      this.docked = false;
      this.view.dockManager.close(this.panel);
    }
    this.clear();
  }

  dispose(): void {
    this.hide();
    for (const s of this.subs)
      s.unsubscribe();
    this.subs = [];
    this.view.detach();
  }
}
