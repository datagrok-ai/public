import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {USE_AUTO_SELECTION, USE_CUSTOM_AXES, SCATTER_ROW_LIM, SIZE, COL_NAME, OPT_TYPE,
  NumericArray,
  LABEL,
  DIFFERENCE} from './defs';
import {getParetoMask} from './pareto-computations';

export class ParetoFrontViewer extends DG.JsViewer {
  private title: string;
  private showTitle: boolean;
  private description: string;
  private descriptionPosition: string;
  private descriptionVisibilityMode: string;
  private minimizeColumnNames: string[];
  private maximizeColumnNames: string[];
  private labelColumnNames: string[];
  private useAutoSelection: boolean = true;
  private xAxisColumnName: string;
  private yAxisColumnName: string;
  private useCustomAxes: boolean;
  private toChangeUseCustomAxes: boolean = true;
  private legendVisibility: string;
  private legendPosition: string;
  private toChangeScatterMarkerSize = false;
  private colorColumnName: string | null = null;

  private scatter = DG.Viewer.scatterPlot(grok.data.demo.demog(), {
    showColorSelector: false,
    showSizeSelector: false,
    autoLayout: false,
  });

  private numCols: DG.Column[] = [];
  private numColNames: string[] = [];
  private numColsCount: number = 0;
  private rowCount: number = 0;
  private isApplicable: boolean = false;
  private errMsg: string = '';
  private resultColName: string = '';
  private sizeColName: string = '';
  private optimizedColNames: string[] = [];

  get type(): string {
    return 'ParetoFrontViewer';
  }

  constructor() {
    super();

    this.title = this.string('title', 'Pareto front');
    this.showTitle = this.bool('showTitle', false, {category: 'Description'});

    this.description = this.string('description');

    this.descriptionPosition = this.string('descriptionPosition', 'Top', {choices: ['Left', 'Right', 'Top', 'Bottom']});

    this.descriptionVisibilityMode = this.string('descriptionVisibilityMode', 'Auto', {
      choices: ['Auto', 'Always', 'Never']},
    );

    this.minimizeColumnNames = this.addProperty('minimizeColumnNames', DG.TYPE.COLUMN_LIST, null, {
      columnTypeFilter: DG.TYPE.NUMERICAL,
      category: 'Data',
      description: 'Columns with features to be minimized during Pareto optimization.',
    });

    this.maximizeColumnNames = this.addProperty('maximizeColumnNames', DG.TYPE.COLUMN_LIST, null, {
      columnTypeFilter: DG.TYPE.NUMERICAL,
      category: 'Data',
      description: 'Columns with features to be maximized during Pareto optimization.',
    });

    this.xAxisColumnName = this.string('xAxisColumnName', null, {
      category: 'Axes',
      description: 'A column to be used on the X axis of the scatter plot.',
      nullable: false,
    });

    this.yAxisColumnName = this.string('yAxisColumnName', null, {
      category: 'Axes',
      description: 'A column to be used on the Y axis of the scatter plot.',
      nullable: false,
    });

    this.useCustomAxes = this.bool('useCustomAxes', USE_CUSTOM_AXES ? true : null, {
      defaultValue: USE_CUSTOM_AXES ? true : null,
      category: 'Axes',
      // eslint-disable-next-line max-len
      description: 'If checked, custom coordinate axes are used. If unchecked, axes are selected automatically based on the optimized features.',
    });

    this.labelColumnNames = this.addProperty('labelColumnsColumnNames', DG.TYPE.COLUMN_LIST, null, {
      category: 'Labels',
      description: 'Label columns to show next to the markers.',
    });

    this.useAutoSelection = this.bool('useAutoSelection', USE_AUTO_SELECTION, {
      category: 'Labels',
      description: 'Automatically select the most relevant columns for the legend.',
      defaultValue: USE_AUTO_SELECTION,
    });

    this.legendVisibility = this.string('legendVisibility', 'Auto', {choices: ['Auto', 'Always', 'Never']});
    this.legendPosition = this.string('legendPosition', 'Top', {
      choices: ['Auto', 'Left', 'Right', 'Top', 'Bottom', 'RightTop', 'RightBottom', 'LeftTop', 'LeftBottom'],
    });

    this.root.append(this.scatter.root);
  } // constructor

  private initializeData() {
    this.isApplicable = this._testColumns();

    if (!this.isApplicable)
      return;

    const cols = this.dataFrame.columns;
    const colList = cols.toList();
    this.numCols = colList.filter((col) => col.isNumerical);
    this.numColNames = this.numCols.map((col) => col.name);
    this.numColsCount = this.numCols.length;
    this.rowCount = this.dataFrame.rowCount;
    this.toChangeScatterMarkerSize = this.rowCount > SCATTER_ROW_LIM;
    this.resultColName = cols.getUnusedName(COL_NAME.OPT);
    this.sizeColName = cols.getUnusedName(COL_NAME.SIZE);
  }

  private computeParetoFront(): void {
    if (!this.isApplicable)
      return;

    const data: NumericArray[] = [];
    const sense: OPT_TYPE[] = [];

    if (this.minimizeColumnNames != null) {
      this.minimizeColumnNames.forEach((name) => {
        data.push(this.dataFrame.col(name)!.getRawData());
        sense.push(OPT_TYPE.MIN);
      });
    }

    if (this.maximizeColumnNames != null) {
      this.maximizeColumnNames.forEach((name) => {
        data.push(this.dataFrame.col(name)!.getRawData());
        sense.push(OPT_TYPE.MAX);
      });
    }

    if (data.length > 0) {
      const mask = getParetoMask(data, sense, this.rowCount);
      const colOpt = DG.Column.fromStrings(this.resultColName, mask.map((res) => res ? LABEL.OPTIMAL : LABEL.NON_OPT));

      let resCol = this.dataFrame.col(this.resultColName);

      if (resCol == null) {
        this.dataFrame.columns.add(colOpt);
        resCol = colOpt;
      } else {
        const newCats = colOpt.categories;
        const newRaw = colOpt.getRawData();

        const prevRaw = resCol.getRawData();
        const prevCats = resCol.categories;
        const indeces = newCats.map((cat) => prevCats.indexOf(cat));

        for (let k = 0; k < this.rowCount; ++k)
          prevRaw[k] = indeces[newRaw[k]];
      }

      this.markResColWithColor(resCol);

      if (this.toChangeScatterMarkerSize) {
        const sizeCol = this.dataFrame.col(this.sizeColName);

        if (sizeCol == null) {
          this.dataFrame.columns.add(DG.Column.fromInt32Array(
            this.sizeColName,
            new Int32Array(mask.map((res) => res ? SIZE.OPTIMAL : SIZE.NON_OPT))),
          );
        } else {
          const raw = sizeCol.getRawData();
          for (let k = 0; k < this.rowCount; ++k)
            raw[k] = mask[k] ? SIZE.OPTIMAL : SIZE.NON_OPT;
        }
      }
    }
  } // computeParetoFront

  private markResColWithColor(col: DG.Column): void {
    col.colors.setCategorical({
      'optimal': '#2ca02c',
      'non-optimal': '#e3e3e3',
    });
  } // markResColWithColor

  _showErrorMessage(msg: string) {this.root.appendChild(ui.divText(msg, 'd4-viewer-error'));}

  _testColumns(): boolean {
    if (this.dataFrame.rowCount < 1) {
      this.errMsg = 'Cannot compute Pareto front: the table is empty.';
      return false;
    }

    if (this.dataFrame.columns.length < 2) {
      this.errMsg = 'Cannot compute Pareto front: at least two numeric columns are required.';
      return false;
    }

    return true;
  } // isApplicable

  private setScatterOptions() {
    if (this.toChangeScatterMarkerSize)
      this.scatter.setOptions({markerMinSize: SIZE.NON_OPT, markerMaxSize: SIZE.OPTIMAL});

    this.scatter.setOptions({
      title: this.title,
      showTitle: this.showTitle,
      description: this.description,
      descriptionPosition: this.descriptionPosition,
      descriptionVisibilityMode: this.descriptionVisibilityMode,
      //labelColumnNames: this.labelColumnNames,
      xColumnName: this.xAxisColumnName,
      yColumnName: this.yAxisColumnName,
      legendVisibility: this.legendVisibility,
      legendPosition: this.legendPosition,
      colorColumnName: this.colorColumnName,
    });
  }

  private init() {
    this.initializeData();
    if (this.isApplicable) {
      this.scatter.dataFrame = this.dataFrame;

      const initColNames = this.numColNames.filter((_, idx) => this.numColsCount - idx - 1 < DIFFERENCE);
      this.setOptions({
        maximizeColumnNames: [],
        minimizeColumnNames: initColNames,
        xAxisColumnName: initColNames[0],
        yAxisColumnName: initColNames[1],
        useCustomAxes: USE_CUSTOM_AXES ? true : null,
      });
    }
  }

  onTableAttached() {
    this.init();

    // Stream subscriptions
    // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    // this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    // this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));

    this.render();
  }

  // Cancel subscriptions when the viewer is detached
  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
    super.detach();
  }

  // Override to handle property changes
  onPropertyChanged(property: DG.Property) {
    if (!this.isApplicable)
      return;

    switch (property.name) {
    case 'xAxisColumnName':
    case 'yAxisColumnName':
      if (this.toChangeUseCustomAxes)
        this.setOptions({useCustomAxes: true});
      this.toChangeUseCustomAxes = true;
      break;

    case 'useCustomAxes':
      this.updateAxesColumnOptions();
      break;
    }

    this.render((property.name === 'minimizeColumnNames') || (property.name === 'maximizeColumnNames'));
  }

  render(computeData = false) {
    if (!this.isApplicable) {
      this.scatter.root.hidden = true;
      this._showErrorMessage(this.errMsg);
      return;
    }

    if (computeData) {
      this.computeParetoFront();
      this.updateOptimizedColNames();
      this.updateAxesColumnOptions();
      this.colorColumnName = (this.optimizedColNames.length < 1) ? null : COL_NAME.OPT;
    }

    this.setScatterOptions();
    this.scatter.root.hidden = false;
  } // render

  private updateOptimizedColNames() {
    this.optimizedColNames = [];

    if (this.minimizeColumnNames != null)
      this.optimizedColNames.push(...this.minimizeColumnNames);

    if (this.maximizeColumnNames != null)
      this.optimizedColNames.push(...this.maximizeColumnNames);
  } // updateOptimizedColNames

  private updateAxesColumnOptions(): void {
    if (this.useCustomAxes)
      return;

    const length = this.optimizedColNames.length;

    if (length < 1)
      return;

    const xIdx = this.optimizedColNames.indexOf(this.xAxisColumnName);
    const yIdx = this.optimizedColNames.indexOf(this.yAxisColumnName);

    if (length > 1) {
      if (xIdx < 0) {
        this.toChangeUseCustomAxes = false;
        this.setOptions({xAxisColumnName: this.optimizedColNames[yIdx !== 0 ? 0 : 1]});
      }

      if (yIdx < 0) {
        this.toChangeUseCustomAxes = false;
        this.setOptions({yAxisColumnName: this.optimizedColNames[xIdx !== 1 ? 1 : 0]});
      }
    } else {
      if ((xIdx < 0) && (yIdx < 0)) {
        this.toChangeUseCustomAxes = false;
        this.setOptions({xAxisColumnName: this.optimizedColNames[0]});
      }
    }
  } // updateAxesColumnOptions
} // ParetoFrontViewer
