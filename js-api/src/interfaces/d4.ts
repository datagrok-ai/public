/// this file was generated automatically from d4 classes declarations


export interface IBarChartSettings {
  /// Determines the rows shown on the scatter plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or (${WEIGHT} / 2) > 100"
  filter: string;

  /// Determines what happens when you click on a bar.
  onClick: keyof typeof RowGroupAction;

  /// Value column. See *Value Aggr Type* for aggregation options.
  value: string;
  valueColumnName: string;

  /// Value aggregation.
  valueAggrType: string;

  /// When true, each outermost bar is of the same width.
  /// This mode is useful for comparing relative value frequency when the *Stack* column is specified.
  relativeValues: boolean;

  /// Indicates whether the "no data" bar should appear
  /// when the *Split* value is not present.
  includeNulls: boolean;

  /// Whether to sort bars *by category* or *by value*.
  /// See also *Bar Sort Order*
  barSortType: string;

  /// Whether the bars should be sorted in ascending or descending order.
  /// See also *Bar Sort Type*.
  barSortOrder: string;

  axisType: keyof typeof AxisType;

  valueMin: number;

  valueMax: number;

  showValueAxis: boolean;

  showValueSelector: boolean;

  orientation: string;

  /// A categorical column to split data on (each bar represents a category)
  split: string;
  splitColumnName: string;

  /// Time unit map function for *Split* (applicable to dates only).
  splitMap: string;

  showCategoryValues: boolean;

  showValuesInsteadOfCategories: boolean;

  showCategorySelector: boolean;

  /// A categorical column to further split data on.
  /// Each category would become a part of the bar resulting from *Split*.
  /// Shows stacked bars only when *Value Aggr Type* is additive.
  stack: string;
  stackColumnName: string;

  /// Time unit map function for *Stack* (applicable to dates only).
  stackMap: string;

  showStackSelector: boolean;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  color: string;
  colorColumnName: string;

  /// Color aggregation type.
  colorAggrType: string;

  invertColorScheme: boolean;

  linearColorScheme: Array<number>;

  /// Whether the selected rows are indicated.
  /// Only works for cumulative aggregations such as count.
  showSelectedRows: boolean;

  showMouseOverRect: boolean;

  /// Show which part is filtered
  /// Only works with RowSource = All
  showFilteredRows: boolean;

  showMouseOverRows: boolean;

  autoLayout: boolean;

  maxCategoryWidth: number;

  categoryValueWidth: number;

  showValueAxisLine: boolean;

  barBorderLineMouseOverWidth: number;

  barBorderLineWidth: number;

  maxBarHeight: number;

  barCornerRadius: number;

  verticalAlign: keyof typeof VerticalAlignType;

  font: string;

  axisFont: string;

  minTextHeight: number;

  backColor: number;

  axisColor: number;

  barColor: number;

  categoryColor: number;

  valueTextColor: number;

  barBorderLineMouseOverColor: number;

  barBorderLineFilteredColor: number;

  barBorderLineColor: number;

  outerMarginLeft: number;

  outerMarginRight: number;

  outerMarginTop: number;

  outerMarginBottom: number;

  showAllCats: boolean;

  useSplitColors: boolean;

  showEmptyBars: boolean;

  showLabels: string;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  // Properties common for all viewers
  // todo: use code generation
  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  /// Namespace-qualified function that gets executed when a viewer is initialized
  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum RowSet {
  All = 'All',
  Filtered = 'Filtered',
  Selected = 'Selected',
  SelectedOrCurrent = 'SelectedOrCurrent',
  FilteredSelected = 'FilteredSelected',
  MouseOverGroup = 'MouseOverGroup',
  CurrentRow = 'CurrentRow',
  MouseOverRow = 'MouseOverRow',
}

export enum RowGroupAction {
  Select = 'Select',
  Filter = 'Filter',
  None = 'None',
}

export enum AxisType {
  linear = 'linear',
  logarithmic = 'logarithmic',
}

export enum VerticalAlignType {
  Center = 'Center',
  Top = 'Top',
  Bottom = 'Bottom',
}

export enum VisibilityMode {
  Auto = 'Auto',
  Always = 'Always',
  Never = 'Never',
}

export enum FlexAutoPosition {
  Auto = 'Auto',
  Left = 'Left',
  Right = 'Right',
  Top = 'Top',
  Bottom = 'Bottom',
  RightTop = 'RightTop',
  RightBottom = 'RightBottom',
  LeftTop = 'LeftTop',
  LeftBottom = 'LeftBottom',
}

export enum FlexPosition {
  Left = 'Left',
  Right = 'Right',
  Top = 'Top',
  Bottom = 'Bottom',
}

export interface IBoxPlotSettings {
  categoryColumnNames: Array<string>;

  /// If *Category 2* is not defined, sets *Markers Color* the same as the currently selected column.
  category1: string;
  category1ColumnName: string;

  /// Time unit map function for *Category 1 Column Names* (applicable to dates only).
  category1Map: string;

  /// If defined, sets *Markers Color* the same as the currently selected column.
  category2: string;
  category2ColumnName: string;

  /// Time unit map function for *Category 2 Column Names* (applicable to dates only).
  category2Map: string;

  showStatistics: boolean;

  showCategoryAxis: boolean;

  showCategorySelector: boolean;

  labelOrientation: keyof typeof TextOrientation;

  /// Display subcategories - category combinations in the x axis table.
  showMinorCategories: boolean;

  value: string;
  valueColumnName: string;

  axisType: keyof typeof AxisType;

  valueMin: number;

  valueMax: number;

  invertYAxis: boolean;

  showValueAxis: boolean;

  showValueSelector: boolean;

  /// If checked, focuses on filtered points (on value axis) as the filter changes.
  /// If current value column is filtered, zoom is not applied.
  zoomValuesByFilter: boolean;

  /// Include plots, which are empty or have null values.
  showEmptyCategories: boolean;

  /// Column to color-code boxes (Q2-Q3 region) or inner violin shapes.
  /// See also *Bin Color Aggr Type*.
  binColor: string;
  binColorColumnName: string;

  /// Aggregation function for color-coding.
  /// See also *Bin Color*.
  binColorAggrType: string;

  showColorSelector: boolean;

  /// Column to color-code markers. Changing *Category 1* or *Category 2*
  /// sets the color scheme to categorical (same as selected category column).
  markerColor: string;
  markerColorColumnName: string;

  /// Categorical coloring time unit map function (applicable to dates only).
  markerColorMap: string;

  colorAxisType: keyof typeof AxisType;

  invertColorScheme: boolean;

  colorMin: number;

  colorMax: number;

  markers: string;
  markersColumnName: string;

  /// Marker category time unit map function (applicable to dates only).
  markersMap: string;

  markerMinSize: number;

  markerMaxSize: number;

  showSizeSelector: boolean;

  markerSizeColumnName: string;

  markerType: string;

  markerSize: number;

  markerOpacity: number;

  showMeanCross: boolean;

  showLowerDash: boolean;

  showUpperDash: boolean;

  showMedianDash: boolean;

  /// Points are not shown if the number of rows is greater than *Show Values Limit*.
  showValuesLimit: number;

  /// Show points inside the interquartile range (Q3 - Q1).
  showInsideValues: boolean;

  /// Show points outside the interquartile range (Q3 - Q1).
  showOutsideValues: boolean;

  /// Show p-value. Press T to toggle.
  /// Currently works only when there are two categories.
  /// Welch's t-test is used for calculating the p-value.
  showPValue: boolean;

  showMouseOverPoint: boolean;

  showMouseOverRowGroup: boolean;

  statistics: Array<string>;

  showTotalCount: boolean;

  /// Shown values count inside lower and upper bounds, where:
  /// IQR = Q3 - Q1.
  /// Lower Bound = Q1 - (1.5 * IQR).
  /// Upper Bound = Q3 + (1.5 * IQR).
  showInliersCount: boolean;

  /// Shown values count outside lower and upper bounds, where:
  /// IQR = Q3 - Q1.
  /// Lower Bound = Q1 - (1.5 * IQR).
  /// Upper Bound = Q3 + (1.5 * IQR).
  showOutliersCount: boolean;

  showMin: boolean;

  showMax: boolean;

  showAvg: boolean;

  showMed: boolean;

  showStdev: boolean;

  showQ1: boolean;

  showQ3: boolean;

  viewport: string;

  autoLayout: boolean;

  /// When checked, both categories from **Category 1** and **Category 2** columns are shown ignoring axis size adjustments.
  showAllCategories: boolean;

  plotStyle: string;

  axisFont: string;

  categoryFont: string;

  statisticsFont: string;

  whiskerLineWidth: number;

  interquartileLineWidth: number;

  whiskerWidthRatio: number;

  axisUseColumnFormat: boolean;

  /// Number of KDE bins to display a violin plot.
  bins: number;

  whiskerColor: number;

  violinWhiskerColor: number;

  backColor: number;

  filteredRowsColor: number;

  filteredOutRowsColor: number;

  selectedRowsColor: number;

  missingValueColor: number;

  defaultBoxColor: number;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  /// Controls box plot tooltip visibility
  showTooltip: string;

  showLabels: keyof typeof VisibilityMode;

  /// Newline-separated list of column names to be used in a tooltip.
  /// Requires *showTooltip* to be enabled.
  rowTooltip: string;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum TextOrientation {
  Auto = 'Auto',
  Horz = 'Horz',
  Vert = 'Vert',
}

export interface ICalendarSettings {
  date: string;
  dateColumnName: string;

  controlsFont: string;

  showHeader: boolean;

  redWeekends: boolean;

  /// Determines what happens when you click a date.
  onClick: keyof typeof RowGroupAction;

  showFilteredOnly: boolean;

  backColor: number;

  oddMonthColor: number;

  evenMonthColor: number;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ICardSettings {
  caption: string;

  valueSourceType: keyof typeof CardValueSourceType;

  /// Source-type specific value.
  value: string;

  format: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum CardValueSourceType {
  Constant = 'Constant',
  TableMeta = 'TableMeta',
  ColumnMeta = 'ColumnMeta',
  Markup = 'Markup',
  Formula = 'Formula',
  External = 'External',
}

export interface IConfusionMatrixSettings {
  /// Column to be put on the X axis
  x: string;
  xColumnName: string;

  /// Column to be put on the Y axis
  y: string;
  yColumnName: string;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ICorrelationPlotSettings {
  controlsFont: string;

  /// Columns to be put on the X axis
  xColumnNames: Array<string>;

  /// Columns to be put on the Y axis
  yColumnNames: Array<string>;

  correlationType: keyof typeof CorrelationType;

  /// Shows the Pearson correlation coefficient inside the corresponding cell.
  showPearsonR: boolean;

  /// Shows the tooltip with the corresponding scatter plot inside.
  showTooltip: boolean;

  /// Ignores double click behavior on the grid cells.
  ignoreDoubleClick: boolean;

  backColor: number;

  defaultCellFont: string;

  colHeaderFont: string;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum CorrelationType {
  Pearson = 'Pearson',
  Spearman = 'Spearman',
}

export interface IDensityPlotSettings {
  binShape: string;

  /// Columns to be put on the X axis
  x: string;
  xColumnName: string;

  /// Columns to be put on the Y axis
  y: string;
  yColumnName: string;

  autoLayout: boolean;

  axisFont: string;

  showColorScale: boolean;

  invertColorScheme: boolean;

  colorTransformType: keyof typeof AxisType;

  linearColorScheme: Array<number>;

  showXAxis: boolean;

  showYAxis: boolean;

  xAxisType: keyof typeof AxisType;

  yAxisType: keyof typeof AxisType;

  invertXAxis: boolean;

  invertYAxis: boolean;

  showXSelector: boolean;

  showYSelector: boolean;

  xMin: number;

  yMin: number;

  xMax: number;

  yMax: number;

  bins: number;

  allowZoom: boolean;

  binToRange: boolean;

  showBinSelector: boolean;

  backColor: number;

  viewport: string;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IFiltersSettings {
  active: boolean;

  histogramLook: any;

  showFilterCountsIndication: boolean;

  showFilterIndication: boolean;

  /// Indicate the proportion of the selected rows inside each category.
  showSelectionIndication: boolean;

  showHeader: boolean;

  showHistogram: boolean;

  showMinMax: boolean;

  showSearchBox: boolean;

  // shows the current mouse over row in the table using grey vertical bar in corresponding row
  showMouseOverRow: boolean;

  // shows the current row in the table using green vertical bar in corresponding row
  showCurrentRow: boolean;

  // shouws the mouse over group porportion in the filter (similar to how selection proportion is shown).
  showMouseOverGroupRow: boolean;

  /// Show a filter that represents all boolean columns in the table.
  showBoolCombinedFilter: boolean;

  columnNames: Array<string>;

  filters: Array<{[index: string]: any}>;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IFormSettings {
  /// Determines what gets shown on the form.
  syncMode: string;

  showNavigation: boolean;

  showPrevRowArrow: boolean;

  showNextRowArrow: boolean;

  showRowSelector: boolean;

  showFieldEditor: boolean;

  showDesignEditor: boolean;

  showColumnSelector: boolean;

  showSaveFile: boolean;

  showOpenFile: boolean;

  sketchState: {[index: string]: any};

  columnNames: Array<string>;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IGridSettings {
  controlsFont: string;

  /// Indicates whether the grid is editable.
  /// See also *Show Add New Row Icon*
  allowEdit: boolean;

  /// When [allowEditable] is true, shows the last virtual row that the user can edit.
  /// This row gets appended to the underlying table as soon as any value is entered.
  /// The grid should also be in the editable mode
  showAddNewRowIcon: boolean;

  /// Automatically adds a new row in the end of the dataframe when the last row is edited
  /// The grid should also be in the editable mode
  addNewRowOnLastRowEdit: boolean;

  /// When [allowEditable] is true, allows user to remove the mouse over row.
  /// The grid should also be in the editable mode
  showRemoveRowIcon: boolean;

  showColumnLabels: boolean;

  /// Column header height. If not specified, it is calculated automatically.
  /// See also *Col Labels Orientation*, *Horz Col Labels Height*
  colHeaderHeight: number;

  /// Height of the column labels when the orientation is vertical,
  /// and *Col Header Height* is not specified.
  vertColLabelsHeight: number;

  /// Height of the column labels when the orientation is horizontal,
  /// and *Col Header Height* is not specified.
  horzColLabelsHeight: number;

  /// Applicable only to grid
  rowHeight: number;

  /// Indicates mouse-over row by drawing a vertical stripe on the row header
  showMouseOverRowIndicator: boolean;

  /// Indicates current row with the *Current Row Color*.
  showCurrentRowIndicator: boolean;

  sortByColumnNames: Array<string>;

  sortTypes: Array<boolean>;

  pinnedRowColumnNames: Array<string>;

  pinnedRowValues: Array<string>;

  /// Indicates whether the control is in the grid or heatmap mode.
  /// Typically, you don't need to change it manually.
  isGrid: boolean;

  /// When set to false, default menu appears under the 'Grid' submenu.
  topLevelDefaultMenu: boolean;

  /// Whether items applicable to all viewers (such as Pickup Style) should
  /// be shown in a popup menu. Also requires *Show Context Menu*.
  showDefaultPopupMenu: boolean;

  /// Mouse drag on the data cells selects both rows and columns
  allowBlockSelection: boolean;

  /// Shift+click on a header to select a column
  /// Shift+mouse drag on the headers to select multiple columns
  /// Ctrl+click to invert selection
  /// Ctrl+Shift+click to deselect
  allowColSelection: boolean;

  /// Mouse drag on the rows header selects rows
  /// Reorder rows by dragging them
  allowRowReordering: boolean;

  /// Whether to sort when user double-clicks on the column header
  allowSorting: boolean;

  /// Mouse drag on the rows headers selects rows
  /// Ctrl+click to invert selection
  /// Shift+mouse drag to select multiple rows
  /// Ctrl+Shift+mouse drag to unselect
  allowRowSelection: boolean;

  /// Right-click and drag to pan content
  allowContentPanning: boolean;

  showColumnGroups: boolean;

  showRowHeader: boolean;

  showRowGridlines: boolean;

  /// Whether the "hamburger" menu should be shown for a column
  /// when the mouse is over its header
  allowColumnMenu: boolean;

  /// Automatically scroll column into view when this column becomes current
  autoScrollColumnIntoView: boolean;

  /// Automatically scroll current row into view when it is set from outside
  /// (for instance, as a result of clicking on a point in a scatter plot)
  autoScrollRowIntoView: boolean;

  /// Automatically resize column widths when row height is resized
  autoResizeColumnWidths: boolean;

  showColumnGridlines: boolean;

  /// Reordering columns by dragging the header
  allowColReordering: boolean;

  /// Whether the current object (shown in the context panel) is changed
  /// when you click on a column header.
  allowChangeCurrentObject: boolean;

  /// Whether row (rows) can be dragged out of the grid.
  allowRowDragging: boolean;

  extendLastColumn: boolean;

  /// Resize rows by dragging the border between rows on a row header.
  /// Applicable only to grid.
  allowRowResizing: boolean;

  /// Indicates the way colors are sampled in the heatmap mode when there is not enough
  /// pixels on the screen for each row:
  /// True: each row is draws (but the result is blended and the resulting color might not represent any row)
  /// False: a row is sampled and then drawn as one pixel (but non-sampled rows do not get drawn at all)
  /// Applicable only to heatmap.
  drawEveryRow: boolean;

  /// Whether the context menu is shown
  showContextMenu: boolean;

  /// Whether to show notifications when the user tries to edit a read-only table
  showReadOnlyNotifications: boolean;

  /// Whether to show scrollbars in the heatmap mode
  /// Note that scrollbars will still be visible if they are not "expanded"
  showHeatmapScrollbars: boolean;

  frozenColumns: number;

  showCurrentCellOutline: boolean;

  /// Color-coding that applies to all columns.
  /// Additionally, each column can be individually color-coded.
  colorCoding: keyof typeof GridColorCodingType;

  defaultCellFont: string;

  maxFontSize: number;

  colHeaderFont: string;

  /// Orientation of the column header text.
  /// In spreadsheet mode, it defaults to horizontal no matter how small the columns are.
  /// In heat map mode, it depends on whether the text can fit in the area.
  colLabelsOrientation: keyof typeof TextOrientation;

  /// Resizing column header by dragging the border between the header and the first row
  allowColHeaderResizing: boolean;

  /// Resizing columns by dragging the border between column headers
  allowColResizing: boolean;

  missingValueColor: number;

  selectedRowsColor: number;

  selectedColsColor: number;

  currentRowColor: number;

  mouseOverRowColor: number;

  mouseOverRowStripeColor: number;

  backColor: number;

  colHeaderTextColor: number;

  colHeaderBackColor: number;

  colHeaderMouseOverTextColor: number;

  cellTextColor: number;

  currentCellTextColor: number;

  rowHeaderBackColor: number;

  /// true: colors are scaled based on the global min/max in all numerical columns
  /// false: colors are scaled based on the column min/max.
  /// Applicable only to heatmap.
  globalColorScaling: boolean;

  /// Whether the heatmap should be color-coded (you might want to do it programmatically)
  /// See also [globalColorScaling]
  heatmapColors: boolean;

  /// Controls grid tooltip visibility
  showTooltip: string;

  showLabels: keyof typeof VisibilityMode;

  showCellTooltip: boolean;

  /// Include currently visible columns in a tooltip
  showVisibleColumnsInTooltip: boolean;

  showColumnTooltip: boolean;

  /// Newline-separated list of column names to be used in a tooltip.
  /// Requires *showTooltip* to be enabled.
  rowTooltip: string;

  marginLeft: number;

  marginTop: number;

  marginRight: number;

  marginBottom: number;

  allowStickyMeta: keyof typeof AllowStickyMetaType;

  /// Heatmap horizontal scroll positions (maxRangeValue, minValue, maxValue)
  heatmapHorzScroll: Array<number>;

  /// Heatmap vertical scroll positions (maxRangeValue, minValue, maxValue)
  heatmapVertScroll: Array<number>;

  /// Determines whether newly added columns are added to the grid
  syncNewColumns: boolean;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  columnHeaderTypes: Array<string>;

  cellStyle: any;

  currentRowCellStyle: any;

  columns: Array<any>;

  isHeatmap: boolean;

  maxHeatmapColumns: number;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  allowDynamicMenus: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum GridColorCodingType {
  Auto = 'Auto',
  All = 'All',
  None = 'None',
}

export enum AllowStickyMetaType {
  Auto = 'Auto',
  Always = 'Always',
  Never = 'Never',
}

export interface IGridCellStyle {
  font: string;

  horzAlign: string;

  vertAlign: string;

  /// When defined, overrides the default cell tooltip
  tooltip: string;

  cursor: string;

  textWrap: string;

  /// Marker to be shown when the value does not fit in the cell
  marker: string;

  textColor: number;

  backColor: number;

  marginLeft: number;

  marginRight: number;

  marginTop: number;

  marginBottom: number;

  textVertical: boolean;

  /// Applies to image columns only
  imageScale: number;

  /// Applies to image columns only
  opacity: number;

  clip: boolean;

  /// For 'html' cell types only
  element: any;

  /// When defined, the cell editor becomes a combo box with the specified values
  choices: Array<string>;

}

export interface IHistogramSettings {
  /// Whether the filtered out rows should be shown with the semi-transparent color
  /// See also *Filtered Out Color*
  showFilteredOutRows: boolean;

  /// Allows to filter table using the range slider on the bottom.
  filteringEnabled: boolean;

  /// A numerical column used to calculate the distribution of values.
  value: string;
  valueColumnName: string;

  showXAxis: boolean;

  allowColumnSelection: boolean;

  /// Show bin selector in the left top panel when the mouse is over the histogram
  showBinSelector: boolean;

  ///Show bin values on top of the bars
  showValues: boolean;

  /// Number of bins on the histogram
  bins: number;

  valueMin: number;

  valueMax: number;

  /// A categorical column to split data on (each bar represents a category)
  split: string;
  splitColumnName: string;

  /// Whether the values should be normalized when multiple histograms are shown.
  /// If true, you are comparing distributions; if false, you are comparing absolute values.
  /// Requires *Split Column Name* to be set.
  normalizeValues: boolean;

  /// If true, split are shown as stacked bins
  splitStack: boolean;

  /// Spline tension in case multiple histograms are shown.
  /// Requires *Split Column Name* to be set.
  splineTension: number;

  showYAxis: boolean;

  /// Whether the horizontal axis should be zoomed to the range of the visible bins.
  zoomToRange: boolean;

  /// Whether the values should be normalized to the filter or globally.
  normalizeToFilter: boolean;

  /// Bin the values that are in the filter range.
  binToRange: boolean;

  /// Whether markers should be drown when multiple histograms are shown.
  /// Requires *Split Column Name* to be set.
  showMarkers: boolean;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  color: string;
  colorColumnName: string;

  colorAggrType: string;

  invertColorScheme: boolean;

  linearColorScheme: Array<number>;

  /// Indicates current row as a dot on the horizontal axis
  showCurrentRow: boolean;

  /// Indicates current row as a dot on the horizontal axis
  showMouseOverRow: boolean;

  /// Show the distribution of the values that the mouse is currently over in another viewer.
  showMouseOverRowGroup: boolean;

  /// Whether the distribution should be rendered as bars or as a spline.
  /// When *Split* is defined, histogram always shows splines.
  spline: boolean;

  /// Whether the area below the spline should be filled with the corresponding color.
  /// Only applicable when *spline* is true and *split* is empty
  fillSpline: boolean;

  /// Shows distribution lines for each stack or spline category in the histogram (see also **Split** and  **Split Stack**).
  showDistributionLines: boolean;

  showColumnSelector: boolean;

  showSplitSelector: boolean;

  showRangeSlider: boolean;

  /// Visibility of the free-text inputs for the filter range
  showRangeInputs: boolean;

  /// How much space does bin occupy (1 = no margins, 0 = no bin)
  binWidthRatio: number;

  showHistogram: boolean;

  /// Shows the context menu.
  showContextMenu: boolean;

  //style
  autoLayout: boolean;

  xAxisHeight: number;

  yAxisWidth: number;

  filterHeight: number;

  rowIndicatorSize: number;

  backColor: number;

  axisFont: string;

  filteredBinsColor: number;

  selectedBinsColor: number;

  valuesColor: number;

  filteredOutColor: number;

  showCharts: boolean;

  viewport: string;

  marginLeft: number;

  marginTop: number;

  marginRight: number;

  marginBottom: number;

  filterMarginTop: number;

  filterMarginBottom: number;

  aggTooltipColumns: string;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ILineChartSettings {
  /// Deprecated, use splitColumnNames instead
  split: string;
  splitColumnName: string;

  /// A categorical column by which lines are split
  splitColumnNames: Array<string>;

  /// Defines a Y column for the chart on the bottom used for zooming
  overview: string;
  overviewColumnName: string;

  /// Aggregation types for all columns
  yAggrTypes: Array<string>;

  /// When true, X axis is synchronized with the corresponding filter's range values.
  /// Otherwise, when the filter is changed points are filtered out on a chart but the min-max stays.
  axesFollowFilter: boolean;

  packCategories: boolean;

  /// When true, multiple *Y Columns* charts get rendered on top of each other,
  /// otherwise they are stacked
  multiAxis: boolean;

  /// Column to be used on the X axis
  x: string;
  xColumnName: string;

  /// Time unit map function for x column (applicable to dates only).
  xMap: string;

  xAxisType: keyof typeof AxisType;

  /// When defined, background is colored according to the segment column.
  /// Example: time series data with the "stimuli" column
  segment: string;
  segmentColumnName: string;

  invertXAxis: boolean;

  showXAxis: boolean;

  showXSelector: boolean;

  showVerticalGridLines: boolean;

  xAxisLabelOrientation: string;

  xAxisTickmarksMode: keyof typeof AxisTickmarksMode;

  xMin: number;

  xMax: number;

  yMin: number;

  yMax: number;

  /// Numerical columns to be used on Y axes.
  /// Depending on the *
  yColumnNames: Array<string>;

  yAxisType: keyof typeof AxisType;

  showYAxis: boolean;

  yGlobalScale: boolean;

  showHorizontalGridLines: boolean;

  /// Axis title to be shown on the left axis in multi-axis mode
  yAxisTitle: string;

  /// Axis title to be shown on the left axis in multi-axis mode
  y2AxisTitle: string;

  yAxisTickmarksMode: keyof typeof AxisTickmarksMode;

  showYSelectors: boolean;

  showAggrSelectors: boolean;

  showSplitSelector: boolean;

  interpolation: keyof typeof LineInterpolationMode;

  splineTension: number;

  /// A categorical column that determines the shape of the markers.
  markers: string;
  markersColumnName: string;

  /// Marker category time unit map function (applicable to dates only).
  markersMap: string;

  markerType: string;

  markerSize: number;

  /// A boolean column that determines whether to show markers.
  showMarkers: keyof typeof VisibilityMode;

  markersVisibility: string;
  markersVisibilityColumnName: string;

  /// Show vertical line reflecting the position of the current row
  /// See also *Current Line Color*
  showCurrentRowLine: boolean;

  /// Determines whether the line is highlighted when you hover over the corresponding category.
  /// Example: "Split by" = "SEX" and you hover over the "Male" category in the filter.
  showMouseOverCategory: boolean;

  overviewAggrType: string;

  /// Shows/hides upper and lower control limits, and [Western Electric rules](https://sentient.cloud/what-are-western-electric-rules-2/).
  showStaticalProcessControl: boolean;

  /// Shows/hides upper and lower control limits.
  showControlLimits: boolean;

  showSigma1: boolean;

  showSigma2: boolean;

  showAverage: boolean;

  lowerControlLimit: number;

  upperControlLimit: number;

  /// Rile 1: One point is more than 3 standard deviations from the mean. When sample(s) is grossly out of control.
  showOutlier: boolean;

  /// Rule 2: Nine (or more) points in a row are on the same side of the mean. Some prolonged "bias" exists.
  showBias: boolean;

  /// Rule 3: Six (or more) points in a row are continually increasing (or decreasing). When a "trend" exists.
  showConsistentTrend: boolean;

  /// Rule 4: Fourteen (or more) points in a row alternate in direction, increasing then decreasing.
  /// This much "oscillation" is beyond noise. Note that the rule is considered with directionality only.
  /// The position of the mean and the size of the standard deviation have no bearing.
  showOscillation: boolean;

  /// Rule 5:Two (or three) out of the three points in a row are more than 2 standard deviations from the mean
  /// in the same direction. There is a medium tendency for samples to be mediumly out of control.
  /// The side of the mean for the third point is unspecified.
  showMediumShift: boolean;

  /// Rule 6: Four (or five) out of five points in a row are more than 1 standard deviation from the mean
  /// in the same direction. There is strong tendency for samples to be slightly out of control.
  /// The side of the mean for the fifth point is unspecified.
  showSustainedShift: boolean;

  /// Rule 7: Fifteen points in a row are all within 1 standard deviation of the mean on either side of the mean.
  /// With 1 standard deviation, greater variation would be expected.
  showSuppressedVariation: boolean;

  /// Show vertical line reflecting the position of the mouse-over row
  /// See also *Mouse Over Line Color*
  showMouseOverRowLine: boolean;

  /// Use column format for axis labels, where possible
  axesUseColumnFormat: boolean;

  showAggrTypeSelector: boolean;

  /// Marker type for showing the distribution of the aggregated values
  /// when multiple values have the same X value
  whiskersType: string;

  overviewType: string;

  /// Show additional chart on the left
  leftPanel: string;

  autoLayout: boolean;

  segmentsFont: string;

  lineWidth: number;

  whiskerWidth: number;

  lineTransparency: number;

  /// Height of the overview chart
  overviewHeight: number;

  histogramWidth: number;

  /// If true, *X Axis Height* is calculated automatically to fit the required precision.
  /// If false, the specified *X Axis Height*
  autoAxisSize: boolean;

  annotationFont: string;

  /// Requires *Auto Axis Size* to be turned off.
  xAxisHeight: number;

  chartTypes: Array<string>;

  lineColoringType: string;

  lineColor: number;

  backColor: number;

  axisLineColor: number;

  axisTextColor: number;

  axisFont: string;

  markerColor: number;

  mouseOverLineColor: number;

  currentLineColor: number;

  gridLineColor: number;

  selectedRowsColor: number;

  statisticalProcessLineColor: number;

  statisticalProcessAreaColor: number;

  statisticalProcessRuleColor: number;

  xAxisMin: number;

  xAxisMax: number;

  yAxisMin: number;

  yAxisMax: number;

  aggrType: string;

  /// Shows top panel with the "Split by" selector
  showTopPanel: boolean;

  /// Show the "x" close icon for each chart
  showCloseLink: boolean;

  xAxisCustomTickmarks: Array<number>;

  yAxisCustomTickmarks: Array<number>;

  /// Controls scatter plot tooltip visibility
  showTooltip: string;

  showLabels: keyof typeof VisibilityMode;

  /// Newline-separated list of column names to be used in a tooltip.
  /// Requires *showTooltip* to be enabled.
  rowTooltip: string;

  rowGroupTooltip: string;

  /// When selected, column names are added to the legend.
  /// Requires *Multi Axis* to be enabled.
  addYColumnsToLegend: boolean;

  autoAdjustMultiAxisLegendPosition: boolean;

  multiAxisLegendPosition: keyof typeof FlexExtendedPosition;

  annotationRegions: string;

  showViewerAnnotationRegions: boolean;

  showDataframeAnnotationRegions: boolean;

  /// Enables lasso region drawing mode (instead of polygon drawing default one).
  lassoTool: boolean;

  categoryCustomColorIndices: Array<number>;

  categoryCustomColors: Array<number>;

  innerChartMarginTop: number;

  innerChartMarginBottom: number;

  outerChartMarginLeft: number;

  outerChartMarginTop: number;

  outerChartMarginRight: number;

  outerChartMarginBottom: number;

  viewport: string;

  formulaLines: string;

  /// Control the visibility of dataframe-originated formula lines.
  /// Edit formula lines by right-clicking and selecting "Tools | Formula Lines" from the popup menu.
  /// Requires the PowerPack plugin.
  showDataframeFormulaLines: boolean;

  /// Control the visibility of dataframe-originated formula lines.
  /// Edit formula lines by right-clicking and selecting "Tools | Formula Lines" from the popup menu.
  /// Requires the PowerPack plugin.
  showViewerFormulaLines: boolean;

  aggTooltipColumns: string;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum AxisTickmarksMode {
  Auto = 'Auto',
  MinMax = 'MinMax',
  Custom = 'Custom',
  AutoAndCustom = 'AutoAndCustom',
}

export enum LineInterpolationMode {
  None = 'None',
  Spline = 'Spline',
}

export enum FlexExtendedPosition {
  LeftTop = 'LeftTop',
  LeftCenter = 'LeftCenter',
  LeftBottom = 'LeftBottom',
  CenterTop = 'CenterTop',
  CenterCenter = 'CenterCenter',
  CenterBottom = 'CenterBottom',
  RightTop = 'RightTop',
  RightCenter = 'RightCenter',
  RightBottom = 'RightBottom',
}

export interface IMapViewerSettings {
  region: string;
  regionColumnName: string;

  color: string;
  colorColumnName: string;

  colorAggrType: string;

  map: string;

  showColorScale: boolean;

  showColorSelector: boolean;

  allowPanZoom: boolean;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IMarkupViewerSettings {
  stretch: boolean;

  content: string;

  mode: keyof typeof TextInterpretationMode;

  /// Whether the rendered html is passed through Grok's [Markup] engine (don't confuse it
  /// with the Markup that might be used for html rendering)
  markupEnabled: boolean;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum TextInterpretationMode {
  None = 'None',
  Html = 'Html',
  Markup = 'Markup',
  Auto = 'Auto',
}

export interface IMatrixPlotSettings {
  controlsFont: string;

  /// Columns to use on the X axis
  xColumnNames: Array<string>;

  /// Column to use on the Y axis
  yColumnNames: Array<string>;

  font: string;

  cellPlotType: string;

  autoLayout: boolean;

  showXAxes: boolean;

  showYAxes: boolean;

  backColor: number;

  innerViewerLook: any;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface INetworkDiagramSettings {
  node1: string;
  node1ColumnName: string;

  node2: string;
  node2ColumnName: string;

  edgeColorColumnName: string;

  edgeColorAggrType: string;

  edgeWidthColumnName: string;

  edgeWidthAggrType: string;

  node1Size: string;
  node1SizeColumnName: string;

  node1SizeAggrType: string;

  node2Size: string;
  node2SizeColumnName: string;

  node2SizeAggrType: string;

  node1ColorColumnName: string;

  node1ColorAggrType: string;

  node2ColorColumnName: string;

  node2ColorAggrType: string;

  node1Image: string;
  node1ImageColumnName: string;

  node2Image: string;
  node2ImageColumnName: string;

  node1Label: string;
  node1LabelColumnName: string;

  node2Label: string;
  node2LabelColumnName: string;

  autoLayout: boolean;

  node1Shape: keyof typeof ShapeType;

  node1Color: number;

  ///put url, or url with placeholders ... to apply individual img
  node1Img: string;

  node1Physics: boolean;

  node2Shape: keyof typeof ShapeType;

  node2Color: number;

  ///put url, or url with placeholders ... to apply individual img
  node2Img: string;

  node2Physics: boolean;

  /// Merge same values from different columns into one node
  mergeNodes: boolean;

  showFilteredOutNodes: boolean;

  filteredOutNodesColor: number;

  Physics: boolean;

  suspendSimulation: boolean;

  missSimulation: boolean;

  improvedLayout: boolean;

  useCommonProperties: boolean;

  useGoogleImage: boolean;

  nodeShape: keyof typeof ShapeType;

  nodeImg: string;

  nodeColor: number;

  edgeColor: number;

  edgeWidth: number;

  showArrows: keyof typeof ArrowType;

  edgesPhysics: boolean;

  selectedColor: number;

  hoverColor: number;

  edgeLinearColorScheme: Array<number>;

  edgeCategoricalColorScheme: Array<number>;

  node1LinearColorScheme: Array<number>;

  node1CategoricalColorScheme: Array<number>;

  node2LinearColorScheme: Array<number>;

  node2CategoricalColorScheme: Array<number>;

  showColumnSelectors: boolean;

  selectRowsOnClick: boolean;

  selectEdgesOnClick: boolean;

  /// Loads external data; called when you double-click on a node.
  /// The specified function gets called with the node value as a single argument.
  /// Its signature: `dataframe expand(dynamic nodeId)`.
  onNodeExpandFunction: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export enum ShapeType {
  none = 'none',
  ellipse = 'ellipse',
  circle = 'circle',
  database = 'database',
  box = 'box',
  dot = 'dot',
  diamond = 'diamond',
  circularImage = 'circularImage',
  image = 'image',
  icon = 'icon',
}

export enum ArrowType {
  none = 'none',
  to = 'to',
  from = 'from',
  middle = 'middle',
  to_middle_from = 'to_middle_from',
  to_from = 'to_from',
  to_middle = 'to_middle',
  middle_from = 'middle_from',
}

export interface IPcPlotSettings {
  /// Whether the filtered out values are shown.
  /// See also *Filtered Out Line Color*
  showFilteredOutLines: boolean;

  /// Columns to use
  columnNames: Array<string>;

  /// Columns where logarithmic axis is used.
  /// Should be a subset of *Column Names*.
  logColumnsColumnNames: Array<string>;

  color: string;
  colorColumnName: string;

  /// Categorical coloring time unit map function (applicable to dates only).
  colorMap: string;

  showColorSelector: boolean;

  colorAxisType: keyof typeof AxisType;

  invertColorScheme: boolean;

  colorMin: number;

  colorMax: number;

  /// Determines the way a value is mapped to the vertical scale.
  /// TRUE: bottom is column minimum, top is column maximum. Use when columns contain values in different units
  /// FALSE: uses the same scale. This lets you compare values across columns
  /// if units are the same (for instance, use it for tracking change over time).'
  normalizeEachColumn: boolean;

  showCurrentLine: boolean;

  showMouseOverLine: boolean;

  showMouseOverRowGroup: boolean;

  /// Either all lines are shown or only current line, mouse over line, selected ones.
  showAllLines: boolean;

  /// Whether the in-chart filters are visible
  showFilters: boolean;

  currentLineWidth: number;

  autoLayout: boolean;

  lineWidth: number;

  mouseOverLineWidth: number;

  minMaxHeight: number;

  transformation: string;

  labelsOrientation: keyof typeof TextOrientation;

  minMaxOrientation: keyof typeof TextOrientation;

  axisFont: string;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  backColor: number;

  selectedRowsColor: number;

  filteredOutLineColor: number;

  missingValueColor: number;

  lineColor: number;

  currentLineColor: number;

  mouseOverLineColor: number;

  whiskerColor: number;

  showDensity: boolean;

  showMinMax: boolean;

  showLabels: boolean;

  maxCategories: number;

  horzMargin: number;

  densityStyle: string;

  showInterquartileRange: boolean;

  showUpperDash: boolean;

  showLowerDash: boolean;

  showMeanCross: boolean;

  showMedian: boolean;

  showCircles: boolean;

  whiskerLineWidth: number;

  interquartileLineWidth: number;

  bins: number;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IPieChartSettings {
  category: string;
  categoryColumnName: string;

  /// Category time unit map function (applicable to dates only).
  categoryMap: string;

  pieSortType: string;

  pieSortOrder: string;

  includeNulls: boolean;

  labelPosition: string;

  autoLayout: boolean;

  segmentAngle: string;
  segmentAngleColumnName: string;

  segmentAngleAggrType: string;

  segmentLength: string;
  segmentLengthColumnName: string;

  segmentLengthAggrType: string;

  /// Action to be performed when you click on a pie
  onClick: keyof typeof RowGroupAction;

  startAngle: number;

  maxRadius: number;

  shift: number;

  outlineLineWidth: number;

  backColor: number;

  outlineColor: number;

  mouseOverOutlineColor: number;

  innerLabelColor: number;

  missingValueColor: number;

  /// Show the label (text or image) inside the pie.
  showLabel: boolean;

  /// Show the percent value inside the pie.
  showPercentage: boolean;

  /// Show the number of values inside the pie.
  showValue: boolean;

  showColumnSelector: boolean;

  /// Highlight part of the pie that corresponds to the mouse-over rows
  showMouseOverRowGroup: boolean;

  /// Highlight selected rows
  showSelectedRows: boolean;

  marginLeft: number;

  marginTop: number;

  marginRight: number;

  marginBottom: number;

  aggTooltipColumns: string;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IPivotViewerSettings {
  showHeader: boolean;

  showCommandBar: boolean;

  pivotColumnNames: Array<string>;

  groupByColumnNames: Array<string>;

  aggregateColumnNames: Array<string>;

  aggregateAggTypes: Array<string>;

  viewerSettings: Array<any>;

  /// Filters dataframe based on current row. works only if rowSource is set to All. See also "rowSource".
  filteringEnabled: boolean;

  gridLook: any;

  allowViewers: boolean;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IPointsViewerSettings {
  /// List of columns to show
  columnNames: Array<string>;

  /// Additional columns to display in the tooltip
  tooltipColumnNames: Array<string>;

  colorScheme: Array<number>;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IRocCurveSettings {
  /// Columns to be put on the X axis
  predictionColumnNames: Array<string>;

  /// Column to be put on the Y axis
  targetColumn: string;

  /// Positive class name
  positiveClass: string;

  /// Select to draw thresholds
  showThreshold: boolean;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IScatterPlotSettings {
  /// Invalid are null values and not positive numbers if axis is logarithmic.
  filterOutInvalid: boolean;

  /// When true, filtered out points are rendered using *Filtered Out Rows Color*.
  showFilteredOutPoints: boolean;

  /// When true, scatter plot will zoom to an area defined by the range filters for X and Y columns,
  /// even if *Zoom And Filter* property is not set to "Zoom by Filter".
  axesFollowFilter: boolean;

  /// Determines the relationship between table filter and scatter plot area:
  /// * No action: they are disconnected
  /// * Filter by zoom: scatter plot acts as a filter; as you zoom in, points get filtered out
  /// * Zoom by filter: scatter plot focuses on the filtered points as the filter changes
  /// * Pack and zoom by filter: removes filtered out categories and focuses on the filtered points as the filter changes.
  zoomAndFilter: string;

  /// A column to use on the X axis. Could be numerical or categorical.
  x: string;
  xColumnName: string;

  /// Time unit map function for x column (applicable to dates only).
  xMap: string;

  /// A column to use on the Y axis. Could be numerical or categorical.
  y: string;
  yColumnName: string;

  /// Time unit map function for y column (applicable to dates only).
  yMap: string;

  xAxisType: keyof typeof AxisType;

  yAxisType: keyof typeof AxisType;

  invertXAxis: boolean;

  invertYAxis: boolean;

  xMin: number;

  yMin: number;

  xMax: number;

  yMax: number;

  showVerticalGridLines: boolean;

  showHorizontalGridLines: boolean;

  showXAxis: boolean;

  showYAxis: boolean;

  showXSelector: boolean;

  showYSelector: boolean;

  /// Shows a distribution histogram along the X axis (at the top)
  showXHistogram: boolean;

  /// Shows a distribution histogram along the Y axis (on the right)
  showYHistogram: boolean;

  /// Number of bins for axis histograms
  histogramBins: number;

  /// Point lower bound for x axis whiskers. Selecting it disables *X Whisker Range*.
  xWhiskerMin: string;
  xWhiskerMinColumnName: string;

  /// Point upper bound for x axis whiskers. Selecting it disables *X Whisker Range*.
  xWhiskerMax: string;
  xWhiskerMaxColumnName: string;

  /// Point range for x axis whiskers. Applied only if *X Whisker Min* and *X Whisker Max* are not set.
  xWhiskerRange: string;
  xWhiskerRangeColumnName: string;

  /// Point lower bound for y axis whiskers. Selecting it disables *Y Whisker Range*.
  yWhiskerMin: string;
  yWhiskerMinColumnName: string;

  /// Point upper bound for y axis whiskers. Selecting it disables *Y Whisker Range*.
  yWhiskerMax: string;
  yWhiskerMaxColumnName: string;

  /// Point range for y axis whiskers. Applied only if *Y Whisker Min* and *Y Whisker Max* are not set.
  yWhiskerRange: string;
  yWhiskerRangeColumnName: string;

  xAxisLabelOrientation: string;

  /// A column to be used for color-coding. Could be numerical or categorical.
  /// If not set, *Filtered Rows Color* is used for markers that pass the filter.
  /// Color palettes could defined either for columns in the column context panel,
  /// or via *Linear Color Scheme* and *Categorical Color Scheme* properties.
  color: string;
  colorColumnName: string;

  /// Categorical coloring time unit map function (applicable to dates only).
  colorMap: string;

  showColorSelector: boolean;

  colorAxisType: keyof typeof AxisType;

  invertColorScheme: boolean;

  colorMin: number;

  colorMax: number;

  /// A numerical column to use for size-coding markers.
  /// See also *Marker Min Size* and *Marker Max Size*.
  size: string;
  sizeColumnName: string;

  showSizeSelector: boolean;

  /// A categorical column that determines the shape of the markers.
  markers: string;
  markersColumnName: string;

  /// Marker category time unit map function (applicable to dates only).
  markersMap: string;

  markerType: string;

  // By default - automatic sizing based on current dataframe
  markerDefaultSize: number;

  markerOpacity: number;

  /// Randomly shift (x, y) marker position up to the *Jitter Size* pixels.
  /// Useful when multiple points fall on the same exact position.
  /// If *Jitter Size Y* is defined, then *Jitter Size* shifts x only.
  jitterSize: number;

  /// Randomly shift y marker position up to the *Jitter Size Y* pixels.
  jitterSizeY: number;

  markerDrawBorder: boolean;

  markerBorderWidth: number;

  markerMinSize: number;

  markerMaxSize: number;

  /// When defined, a line would be drawn for each series (defined by the categorical color column)
  /// using the order specified by "Lines Order"
  linesOrder: string;
  linesOrderColumnName: string;

  /// Defines the width of the lines connecting the markers. See **Lines Width**.
  linesWidth: number;

  /// Label columns to show next to the markers.
  labelColumnNames: Array<string>;

  /// Determines the rows shown on the scatter plot.
  showLabelsFor: keyof typeof RowSet;

  /// Determines how to show marker label:
  /// * Always - show labels for all visible markers
  /// * Auto - show labels only for markers where enough space is available
  /// * Never - show no labels.
  displayLabels: keyof typeof VisibilityMode;

  /// Determines whether to show column names next to label values.
  showLabelNamedColumns: keyof typeof VisibilityMode;

  /// If checked, display a label content as marker.
  useLabelAsMarker: boolean;

  /// To display labels separately or as markers (works for non-text labels).
  labelColorAsMarker: boolean;

  /// Marker size in which label is inscribed.
  labelAsMarkerSize: number;

  /// Label inner content size.
  labelContentSize: number;

  /// Regression line visibility (toggle by pressing R)
  showRegressionLine: boolean;

  showRegressionLineEquation: boolean;

  showSpearmanCorrelation: boolean;

  showPearsonCorrelation: boolean;

  showMeanAbsoluteError: boolean;

  showRootMeanSquareError: boolean;

  regressionPerCategory: boolean;

  /// Control the visibility of dataframe-originated formula lines.
  /// Edit formula lines by right-clicking and selecting "Tools | Formula Lines" from the popup menu.
  /// Requires the PowerPack plugin.
  showDataframeFormulaLines: boolean;

  /// Control the visibility of dataframe-originated formula lines.
  /// Edit formula lines by right-clicking and selecting "Tools | Formula Lines" from the popup menu.
  /// Requires the PowerPack plugin.
  showViewerFormulaLines: boolean;

  /// Controls the indication of the current row
  showCurrentPoint: boolean;

  /// Controls the indication of the mouse-over row
  showMouseOverPoint: boolean;

  /// Highlight 'mouse-over' rows (such as the ones that fall into a histogram bin that
  /// the mouse is currently hovering over).
  showMouseOverRowGroup: boolean;

  /// When true, selected markers are highlighted using the selected rows color.
  /// When false, selected markers use their regular color coding.
  showSelectedRows: boolean;

  /// Shows tickmarks and labels for minimum and maximum value on each axis.
  showMinMaxTickmarks: boolean;

  /// Shows exact X and Y coordinates for the mouse cursor.
  showDropLines: boolean;

  mouseDrag: string;

  /// When true, lasso area selector is used instead of the rectangular one.
  /// Toggle this option by pressing L.
  lassoTool: boolean;

  allowZoom: boolean;

  autoLayout: boolean;

  backColor: number;

  filteredRowsColor: number;

  filteredOutRowsColor: number;

  selectedRowsColor: number;

  missingValueColor: number;

  labelColor: number;

  axisLineColor: number;

  axisTextColor: number;

  gridLineColor: number;

  regressionLineColor: number;

  whiskerColor: number;

  regressionLineTransparency: number;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  /// Determines whether the axes should follow the non-precision-related format (such as "money")
  /// set for the corresponding column.
  axesUseColumnFormat: boolean;

  formulaLines: string;

  annotationRegions: string;

  showViewerAnnotationRegions: boolean;

  showDataframeAnnotationRegions: boolean;

  viewport: string;

  /// Controls scatter plot tooltip visibility
  showTooltip: string;

  showLabels: keyof typeof VisibilityMode;

  /// Controls whether columns on X and Y axes are displayed in tooltip
  /// * Do not add: they are not shown
  /// * Data values only: only they are shown
  /// * Merge: standard behavior
  dataValues: string;

  /// Newline-separated list of column names to be used in a tooltip.
  /// Requires *showTooltip* to be enabled.
  rowTooltip: string;

  rowGroupTooltip: string;

  /// If true, *X Axis Height* and *Y Axis Width* are calculated automatically to fit the required precision.
  /// If false, the specified *X Axis Height* and *Y Axis Width* properties are used.
  autoAxisSize: boolean;

  /// Requires *Auto Axis Size* to be turned off.
  xAxisHeight: number;

  /// Requires *Auto Axis Size* to be turned off.
  yAxisWidth: number;

  axisFont: string;

  labelFont: string;

  formulaFont: string;

  annotationFont: string;

  defaultRenderer: boolean;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IScatterPlot3dSettings {
  x: string;
  xColumnName: string;

  y: string;
  yColumnName: string;

  z: string;
  zColumnName: string;

  size: string;
  sizeColumnName: string;

  color: string;
  colorColumnName: string;

  label: string;
  labelColumnName: string;

  showAxes: boolean;

  xAxisType: keyof typeof AxisType;

  yAxisType: keyof typeof AxisType;

  zAxisType: keyof typeof AxisType;

  backColor: number;

  filteredRowsColor: number;

  filteredOutRowsColor: number;

  selectedRowsColor: number;

  missingValueColor: number;

  axisLineColor: number;

  axisTextColor: number;

  gridLineColor: number;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  dynamicCameraMovement: boolean;

  showVerticalGridLines: boolean;

  showHorizontalGridLines: boolean;

  showFilteredOutPoints: boolean;

  /// Highlight 'mouse-over' rows (such as the ones that fall into a histogram bin that
  /// the mouse is currently hovering over).
  showMouseOverRowGroup: boolean;

  markerType: string;

  markerOpacity: number;

  markerRandomRotation: boolean;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ISparklinesSettings {
  /// List of columns to show aggregations on
  columnNames: Array<string>;

  /// List of aggregations for the columns
  aggregations: Array<string>;

  sparklineType: string;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  color: string;
  colorColumnName: string;

  /// Color aggregation type.
  colorAggrType: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface IStatsViewerSettings {
  controlsFont: string;

  columnNames: Array<string>;

  stats: Array<string>;

  backColor: number;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ISummarySettings {
  /// List of columns to show aggregations on
  columnNames: Array<string>;

  /// List of aggregations for the columns
  aggregations: Array<string>;

  /// Controls the source of the data comparison
  /// * Row: shows vertical bars based on each row category
  /// * Column: shows horizontal bars based on each column category
  /// * Global: shows horizontal bars based on all selected categories
  normalization: string;

  /// Visualization type (text, circles or bars)
  visualization: string;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  color: string;
  colorColumnName: string;

  /// Color aggregation type.
  colorAggrType: string;

  /// Whether to apply color coding to the background or to the text.
  applyTo: string;

  /// Custom color scheme for the color-coding.
  colorSchemes: Array<Array<number>>;

  invertColorScheme: boolean;

  /// If true - show bars on different sides of zero axes for negative and positive values
  zeroAxis: boolean;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ITileViewerSettings {
  controlsFont: string;

  lanesColumnName: string;

  cardMarkup: string;

  allowDragBetweenLanes: boolean;

  /// Whether the form auto-generates whenever columns change
  autoGenerate: boolean;

  sketchState: {[index: string]: any};

  columnsJson: string;

  tilesFont: string;

  lanes: Array<string>;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ITreeMapSettings {
  splitByColumnNames: Array<string>;

  color: string;
  colorColumnName: string;

  colorAggrType: string;

  size: string;
  sizeColumnName: string;

  autoLayout: boolean;

  sizeAggrType: string;

  defaultColor: number;

  showColumnSelectionPanel: boolean;

  outerMarginLeft: number;

  outerMarginRight: number;

  outerMarginTop: number;

  outerMarginBottom: number;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

export interface ITrellisPlotSettings {
  xColumnNames: Array<string>;

  yColumnNames: Array<string>;

  viewerType: string;

  //if false, full screen icon will not be shown on inner viewer hover
  allowViewerFullScreen: boolean;

  yLabelsOrientation: keyof typeof TextOrientation;

  xLabelsOrientation: keyof typeof TextOrientation;

  categoryLabelFont: string;

  innerViewerLook: any;

  innerViewerLooks: {[index: string]: any};

  showGridlines: string;

  showXSelectors: boolean;

  showYSelectors: boolean;

  /// Shows horizontal axis for each inner viewer row.
  /// Depending on selected inner viewer, axis may be hidden automatically.
  showXAxes: keyof typeof VisibilityMode;

  /// Shows vertical axis for each inner viewer column.
  /// Depending on selected inner viewer, axis may be hidden automatically.
  showYAxes: keyof typeof VisibilityMode;

  /// If checked, all inner viewers use the same axis scales.
  globalScale: boolean;

  /// Shows range sliders for inner viewer axis. Available only if corresponding axes (X or Y) are shown.
  showRangeSliders: boolean;

  showXLabels: boolean;

  showYLabels: boolean;

  showControlPanel: boolean;

  syncMouseOverRow: boolean;

  /// Action to be performed when you click on a trellis cell
  onClick: keyof typeof RowGroupAction;

  packCategories: boolean;

  useTiledView: boolean;

  tilesPerRow: number;

  autoLayout: boolean;

  backColor: number;

  viewport: string;

  legendVisibility: keyof typeof VisibilityMode;

  legendPosition: keyof typeof FlexAutoPosition;

  /// Determines the rows shown on the plot.
  rowSource: keyof typeof RowSet;

  /// Formula that filters out rows to show.
  /// Examples:
  ///   ${AGE} > 20 or ${WEIGHT / 2)} > 100,
  ///   ${SEVERITY} == 'Medium',
  ///   ${RACE}.endsWith('sian')
  filter: string;

  /// Viewer controls elements font.
  controlsFont: string;

  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  /// Viewer description that gets shown at the *Descriptor Position*.
  /// Markup is supported.
  description: string;

  /// Help to be shown when user clicks on the '?' icon on top.
  /// Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

  /// Namespace-qualified function that gets executed when a viewer is initialized
  initializationFunction: string;

  /// JavaScript that gets executed after a viewer is initialized and added to the TableView
  onInitializedScript: string;

  descriptionPosition: keyof typeof FlexPosition;

  descriptionVisibilityMode: keyof typeof VisibilityMode;

}

