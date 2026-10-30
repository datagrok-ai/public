import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as echarts from 'echarts';
import {RenderSignals, echartsFramePending, unsubscribeAll} from '../../utils/utils';

export class EChartViewer extends DG.JsViewer {
  private _chart: echarts.ECharts | null = null;
  option: any;
  /** Held by the base's own subscriptions too, for the viewers that report them. */
  protected signals = new RenderSignals(() => this.frameOwed());

  top?: string;
  bottom?: string;
  animationDuration?: number;
  animationDurationUpdate?: number;
  tableName?: string;
  viewerFilter: DG.BitSet | null = null;

  constructor() {
    super();

    //common properties
    this.tableName = this.string('table', null, {fieldName: 'tableName', category: 'Data', editor: 'table'});
    this.addRowSourceAndFormula();
    const chartDiv = ui.div([], {style: {position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}});
    chartDiv.style.cssText += 'overflow: hidden!important;';
    this.root.appendChild(chartDiv);
    const warn = console.warn.bind(console);
    console.warn = () => {};
    try {
      this._chart = echarts.init(chartDiv);
    } finally {
      console.warn = warn;
    }
    this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => {
      this.chart.resize();
      this.signals.settle();
    }));
  }

  protected frameOwed(): boolean {
    return echartsFramePending(this._chart);
  }

  get chart(): echarts.ECharts {
    if (this._chart === null) {
      const chartDiv = this.root.children[0] as HTMLDivElement;
      this._chart = echarts.init(chartDiv);
    }
    return this._chart;
  }

  set chart(value: echarts.ECharts | null) {
    this._chart = value;
  }

  initCommonProperties() {
    /**
     * Removed the 'right' and 'left' properties from the configuration
     * as they disrupt the layout and cause visualization issues in the ECharts viewer.
     */
    this.top = this.string('top', '5px');
    this.bottom = this.string('bottom', '5px');

    this.animationDuration = this.int('animationDuration', 500);
    this.animationDurationUpdate = this.int('animationDurationUpdate', 750);
  }

  onTableAttached() {
    this.addSelectionOrDataSubs();
    this.render();
  }

  onSourceRowsChanged() {
    this.render();
  }

  addSelectionOrDataSubs() {
    this.subs.push(this.signals.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(this.signals.debounce(this.dataFrame.onDataChanged, 50).subscribe((_) => this.render()));
  }

  protected resubscribe(add: () => void): void {
    unsubscribeAll(this.subs);
    add();
  }

  protected filterByClick(match: (index: number) => boolean, event: MouseEvent): void {
    this.viewerFilter ??= DG.BitSet.create(this.dataFrame.rowCount);
    this.viewerFilter.handleClick(match, event);
    if (this.viewerFilter.trueCount === 0)
      this.viewerFilter = null;
    this.dataFrame.rows.requestFilter();
  }

  protected addViewerFilterSub(): void {
    this.subs.push(this.dataFrame.onRowsFiltering.subscribe(() => {
      if (this.viewerFilter)
        this.dataFrame.filter.and(this.viewerFilter);
    }));
  }

  protected clearViewerFilter(): void {
    if (this.viewerFilter === null)
      return;
    this.viewerFilter = null;
    this.dataFrame.rows.requestFilter();
  }

  detach(): void {
    this.clearViewerFilter();
    this.signals.reset();
    super.detach();
  }

  prepareOption() {}

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    const properties = p !== null ? [p] : this.props.getProperties();

    for (const p of properties)
      this.option.series[0][p.name] = p.get(this);

    if (render)
      this.chart.setOption(this.option);
  }

  render() {
    this.option.series[0].data = this.getSeriesData();
    this.chart.setOption(this.option);
  }

  getSeriesData() {
    throw new Error('Method not implemented.');
  }

  updateTable() {
    const dataFrame = grok.shell.tables.find((df: DG.DataFrame) => df.name === this.tableName);
    if (dataFrame) {
      this.clearViewerFilter();
      this.dataFrame = dataFrame;
    }
  }
}
