import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {_package} from '../package';
import {LAB_RES_N, LAB_TEST, SUBJECT_ID, VS_TEST, VS_RES_N, VISIT_DAY_STR,
  BW_TEST, BW_RES_N, BG_TEST, BG_RES_N,
  VISIT} from '../constants/columns-constants';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {CDISC_STANDARD} from '../utils/types';
import {studies} from '../utils/app-utils';


export class MatrixesView extends ClinicalCaseViewBase {
  plotViewer: DG.Viewer | null = null;
  corrPlotDiv = ui.box();
  matrixDataframe: DG.DataFrame;
  matrixTableView: DG.TableView | null = null;
  matrixFiltersGroup: HTMLElement | null = null;
  domains = ['lb', 'vs', 'bw', 'bg'];
  domainFields = {'lb': {'test': LAB_TEST, 'res': LAB_RES_N}, 'vs': {'test': VS_TEST, 'res': VS_RES_N},
    'bw': {'test': BW_TEST, 'res': BW_RES_N}, 'bg': {'test': BG_TEST, 'res': BG_RES_N}};
  initialDataframe: DG.DataFrame;
  isSend = false;
  xColumns: DG.Column[] = [];
  yColumns: DG.Column[] = [];
  plotType: 'correlation' | 'matrix' = 'correlation';

  constructor(name, studyId) {
    super(name, studyId);
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/correlation_matrix.md`;
  }

  loaded = false;

  createView(): void {
    this.isSend = studies[this.studyId].config.standard === CDISC_STANDARD.SEND;
    this.domains = this.domains.filter((it) => studies[this.studyId].domains[it] !== null &&
      !this.optDomainsWithMissingCols.includes(it));
    this.domains.forEach((it) => {
      if (this.isSend)
        createVisitDayStrCol(studies[this.studyId].domains[it]);

      const df = studies[this.studyId].domains[it].clone(null, [SUBJECT_ID,
        this.isSend ? VISIT_DAY_STR : VISIT,
        this.domainFields[it]['test'], this.domainFields[it]['res']]);
      df.getCol(this.domainFields[it]['test']).name = 'test';
      df.getCol(this.domainFields[it]['res']).name = 'res';
      if (!this.initialDataframe)
        this.initialDataframe = df;
      else
        this.initialDataframe.append(df, true);
    });
    this.createCorrelationMatrixDataframe(this.initialDataframe);

    // Set default values for xColumns and yColumns as first 5 columns (or less if fewer available)
    const availableColumnNames = this.matrixDataframe.columns.names()
      .filter((name) => name !== SUBJECT_ID && name !== (this.isSend ? VISIT_DAY_STR : VISIT));
    const defaultColumnNames = availableColumnNames.slice(0, Math.min(5, availableColumnNames.length));
    this.xColumns = defaultColumnNames.map((name) => this.matrixDataframe.col(name)).filter((col) => col !== null);
    this.yColumns = defaultColumnNames.map((name) => this.matrixDataframe.col(name)).filter((col) => col !== null);

    const plotTypeChoice = ui.input.choice('Plot', {
      value: this.plotType,
      items: ['correlation', 'matrix'],
      onValueChanged: () => {
        this.plotType = plotTypeChoice.value as 'correlation' | 'matrix';
        this.switchPlot();
      },
    });

    const filterIcon = ui.icons.filter(() => {
      if (!this.matrixFiltersGroup) {
        grok.shell.warning('Filters are not available yet');
        return;
      }
      ui.showPopup(ui.div(this.matrixFiltersGroup), filterIcon, {vertical: true});
    }, 'Matrix filters');

    this.root.className = 'grok-view ui-box';
    this.root.append(this.corrPlotDiv);
    // this.root.style.marginTop = '15px';
    this.setRibbonPanels([
      [
        plotTypeChoice.root,
        filterIcon,
      ],
    ]);
    this.switchPlot();
    grok.shell.o = this.propertyPanel();
  }


  private getFilters() {
    const filters = DG.Viewer.fromType(DG.VIEWER.FILTERS, this.matrixDataframe, {
      'columnNames': this.matrixDataframe.columns.names(),
      'showContextMenu': true,
    }).root;
    return ui.div(filters, 'clinical-case-correlation-filters');
  }

  private createCorrelationMatrixDataframe(df: DG.DataFrame) {
    this.matrixDataframe = df
      .groupBy([SUBJECT_ID, this.isSend ? VISIT_DAY_STR : VISIT])
      .pivot('test')
      .avg('res')
      .aggregate();

    // Create table view without adding to workspace
    this.matrixTableView = DG.TableView.create(this.matrixDataframe, false);

    this.matrixFiltersGroup = this.getFilters();
  }

  private switchPlot() {
    // Clear existing plot
    ui.empty(this.corrPlotDiv);
    ui.setUpdateIndicator(this.corrPlotDiv, true, `Loading ${this.plotType} plot...`);

    const viewerType = this.plotType === 'correlation' ? DG.VIEWER.CORR_PLOT : DG.VIEWER.MATRIX_PLOT;
    const plotOptions: any = {
      xColumnNames: this.xColumns.map((col) => col.name),
      yColumnNames: this.yColumns.map((col) => col.name),
    };

    this.matrixDataframe.plot
      .fromType(viewerType, plotOptions).then((v: any) => {
        this.plotViewer = v;
        this.root.className = 'grok-view ui-box';
        this.corrPlotDiv.append(this.plotViewer.root);
      }).finally(() => ui.setUpdateIndicator(this.corrPlotDiv, false));
  }

  override async propertyPanel() {
    if (!this.matrixDataframe)
      return;

    const acc = this.createAccWithTitle(this.name);

    const xColumnsInput = ui.input.columns('X', {
      table: this.matrixDataframe,
      value: this.xColumns,
    });

    const yColumnsInput = ui.input.columns('Y', {
      table: this.matrixDataframe,
      value: this.yColumns,
    });

    const updatePlotColumns = () => {
      this.xColumns = xColumnsInput.value;
      this.yColumns = yColumnsInput.value;
      this.plotViewer.setOptions({
        xColumnNames: this.xColumns.map((col) => col.name),
        yColumnNames: this.yColumns.map((col) => col.name),
      });
    };

    xColumnsInput.onChanged.subscribe(() => updatePlotColumns());
    yColumnsInput.onChanged.subscribe(() => updatePlotColumns());

    acc.addPane('Columns', () => ui.divV([
      xColumnsInput.root,
      yColumnsInput.root,
    ]));

    return acc.root;
  }
}
