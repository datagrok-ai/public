/// this file was generated automatically from d4 classes declarations


export interface IScatterPlot3dLookSettings {
  xColumnName: string;

  yColumnName: string;

  zColumnName: string;

  sizeColumnName: string;

  colorColumnName: string;

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

  showXAxis: boolean;

  showYAxis: boolean;

  showZAxis: boolean;

  showFilteredOutPoints: boolean;

  /// Highlight 'mouse-over' rows (such as the ones that fall into a histogram bin that
  /// the mouse is currently hovering over).
  showMouseOverRowGroup: boolean;

  markerType: string;

  markerOpacity: number;

  markerRandomRotation: boolean;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface ITreeMapLookSettings {
  splitByColumnNames: Array<string>;

  colorColumnName: string;

  colorAggrType: string;

  sizeColumnName: string;

  sizeAggrType: string;

  defaultColor: number;

  showColumnSelectionPanel: boolean;

  outerMarginLeft: number;

  outerMarginRight: number;

  outerMarginTop: number;

  outerMarginBottom: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IHistogramLookSettings {
  /// A numerical column used to calculate the distribution of values.
  valueColumnName: string;

  /// A categorical column to split data on (each bar represents a category)
  splitColumnName: string;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  colorColumnName: string;

  colorAggrType: string;

  /// Number of bins on the histogram
  bins: number;

  backColor: number;

  filteredBinsColor: number;

  selectedBinsColor: number;

  filteredOutColor: number;

  /// Allows to filter table using the range slider on the bottom.
  /// This option also controls slider visibility.
  filteringEnabled: boolean;

  showHistogram: boolean;

  /// Indicates current row as a dot on the horizontal axis
  showCurrentRow: boolean;

  /// Indicates current row as a dot on the horizontal axis
  showMouseOverRow: boolean;

  /// Whether the filtered out rows should be shown with the semi-transparent color
  /// See also *Filtered Out Color*
  showFilteredOutRows: boolean;

  /// Show bin selector in the left top panel when the mouse is over the histogram
  showBinSelector: boolean;

  showColumnSelector: boolean;

  showRangeSlider: boolean;

  /// Visibility of the free-text inputs for the filter range
  showRangeInputs: boolean;

  /// Show the distribution of the values that the mouse is currently over in another viewer.
  showMouseOverRowGroup: boolean;

  showXAxis: boolean;

  showYAxis: boolean;

  /// Shows the context menu.
  showContextMenu: boolean;

  showCharts: boolean;

  /// Spline tension in case multiple histograms are shown.
  /// Requires *Split Column Name* to be set.
  splineTension: number;

  allowColumnSelection: boolean;

  rowIndicatorSize: number;

  marginLeft: number;

  marginTop: number;

  marginRight: number;

  marginBottom: number;

  filterMarginTop: number;

  filterMarginBottom: number;

  /// How much space does bin occupy (1 = no margins, 0 = no bin)
  binWidthRatio: number;

  filterHeight: number;

  xAxisHeight: number;

  yAxisWidth: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IFiltersLookSettings {
  showFilterCountsIndication: boolean;

  showFilterIndication: boolean;

  /// Indicate the proportion of the selected rows inside each category.
  showSelectionIndication: boolean;

  showHeader: boolean;

  showHistogram: boolean;

  showSearchBox: boolean;

  /// Show a filter that represents all boolean columns in the table.
  showBoolCombinedFilter: boolean;

  columnNames: Array<string>;

  filters: Array<Map<string, any>>;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IScatterPlotLookSettings {
  /// A column to use on the X axis. Could be numerical or categorical.
  xColumnName: string;

  /// A column to use on the Y axis. Could be numerical or categorical.
  yColumnName: string;

  /// A numerical column to use for size-coding markers.
  /// See also *Marker Min Size* and *Marker Max Size*.
  sizeColumnName: string;

  /// A column to be used for color-coding. Could be numerical or categorical.
  /// If not set, *Filtered Rows Color* is used for markers that pass the filter.
  /// Color palettes could defined either for columns in the column properties panel,
  /// or via *Linear Color Scheme* and *Categorical Color Scheme* properties.
  colorColumnName: string;

  /// Labels to show next to the markers.
  labelsColumnName: string;

  /// A categorical column that determines the shape of the markers.
  markersColumnName: string;

  backColor: number;

  filteredRowsColor: number;

  filteredOutRowsColor: number;

  selectedRowsColor: number;

  missingValueColor: number;

  labelColor: number;

  axisLineColor: number;

  axisTextColor: number;

  gridLineColor: number;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  /// When true, scatter plot will zoom to an area defined by the range filters for X and Y columns,
  /// even if *Zoom And Filter* property is not set to "Zoom by Filter".
  axesFollowFilter: boolean;

  /// Determines the relationship between table filter and scatter plot area:
  /// * No action: they are disconnected
  /// * Filter by zoom: scatter plot acts as a filter; as you zoom in, points get filtered out
  /// * Zoom by filter: scatter plot focuses on the filtered points as the filter changes
  /// * Zoom by filter: removes filtered out categories and focuses on the filtered points as the filter changes.
  zoomAndFilter: string;

  showVerticalGridLines: boolean;

  showHorizontalGridLines: boolean;

  labelColorAsMarker: boolean;

  /// Determines whether the axes should follow the non-precision-related format (such as "money")
  /// set for the corresponding column.
  axesUseColumnFormat: boolean;

  showXAxis: boolean;

  showYAxis: boolean;

  /// Shows tickmarks and labels for minimum and maximum value on each axis.
  showMinMaxTickmarks: boolean;

  showXSelector: boolean;

  showYSelector: boolean;

  showSizeSelector: boolean;

  showColorSelector: boolean;

  /// When true, filtered out points are rendered using *Filtered Out Rows Color*.
  showFilteredOutPoints: boolean;

  /// Controls the indication of the mouse-over row
  showMouseOverPoint: boolean;

  /// Controls the indication of the current row
  showCurrentPoint: boolean;

  /// Regression line visibility (toggle by pressing R)
  showRegressionLine: boolean;

  /// Control the visibility of dataframe-originated formula lines.
  /// Edit formula lines by right-clicking and selecting "Tools | Formula Lines" from the popup menu.
  /// Requires the PowerPack plugin.
  showDataframeFormulaLines: boolean;

  /// Control the visibility of dataframe-originated formula lines.
  /// Edit formula lines by right-clicking and selecting "Tools | Formula Lines" from the popup menu.
  /// Requires the PowerPack plugin.
  showViewerFormulaLines: boolean;

  /// When true, lasso area selector is used instead of the rectangular one.
  /// Toggle this option by pressing L.
  lassoTool: boolean;

  /// Determines the rows shown on the scatter plot.
  // Commented until we make UI for formula editor:
  // @Prop(editor: 'editFormula')
  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or ${WEIGHT / 2) > 100"
  filter: string;

  formulaLines: string;

  viewport: string;

  /// Newline-separated list of column names to be used in a tooltip
  rowTooltip: string;

  rowGroupTooltip: string;

  invertXAxis: boolean;

  invertYAxis: boolean;

  xAxisLabelOrientation: string;

  /// If true, *X Axis Height* and *Y Axis Width* are calculated automatically to fit the required precision.
  /// If false, the specified *X Axis Height* and *Y Axis Width* properties are used.
  autoAxisSize: boolean;

  /// Requires *Auto Axis Size* to be turned off.
  xAxisHeight: number;

  /// Requires *Auto Axis Size* to be turned off.
  yAxisWidth: number;

  /// Highlight 'mouse-over' rows (such as the ones that fall into a histogram bin that
  /// the mouse is currently hovering over).
  showMouseOverRowGroup: boolean;

  markerType: string;

  markerOpacity: number;

  markerDrawBorder: boolean;

  markerBorderWidth: number;

  markerMinSize: number;

  markerMaxSize: number;

  markerDefaultSize: number;

  /// Randomly shift (x, y) marker position up to the *Jitter Size* pixels.
  /// Useful when multiple points fall on the same exact position.
  jitterSize: number;

  axisFont: string;

  labelFont: string;

  /// Shows exact X and Y coordinates for the mouse cursor.
  showDropLines: boolean;

  defaultRenderer: boolean;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface ILineChartLookSettings {
  /// Column to be used on the X axis
  xColumnName: string;

  /// Numerical columns to be used on Y axes.
  /// Depending on the *
  yColumnNames: Array<string>;

  /// Aggregation types for all columns
  yAggrTypes: Array<string>;

  /// A categorical column by which lines are split
  splitColumnName: string;

  /// When defined, background is colored according to the segment column.
  /// Example: time series data with the "stimuli" column
  segmentColumnName: string;

  /// Defines a Y column for the chart on the bottom used for zooming
  overviewColumnName: string;

  filter: string;

  /// When true, multiple *Y Columns* charts get rendered on top of each other,
  /// otherwise they are stacked
  multiAxis: boolean;

  lineWidth: number;

  splineTension: number;

  chartTypes: Array<string>;

  lineColoringType: string;

  lineColor: number;

  backColor: number;

  axisLineColor: number;

  axisTextColor: number;

  markerColor: number;

  mouseOverLineColor: number;

  currentLineColor: number;

//  int innerChartMarginLeft = 5;
//  int innerChartMarginRight = 0;
  innerChartMarginTop: number;

  innerChartMarginBottom: number;

  outerChartMarginLeft: number;

  outerChartMarginTop: number;

  outerChartMarginRight: number;

  outerChartMarginBottom: number;

  xAxisMin: number;

  xAxisMax: number;

  yAxisMin: number;

  yAxisMax: number;

  aggrType: string;

  /// Marker type for showing the distribution of the aggregated values
  /// when multiple values have the same X value
  whiskersType: string;

  /// When true, X axis is synchronized with the corresponding filter's range values.
  /// Otherwise, when the filter is changed points are filtered out on a chart but the min-max stays.
  axesFollowFilter: boolean;

  /// Use column format for axis labels, where possible
  axesUseColumnFormat: boolean;

  yGlobalScale: boolean;

  /// Axis title to be shown on the left axis in multi-axis mode
  yAxisTitle: string;

  /// Axis title to be shown on the left axis in multi-axis mode
  y2AxisTitle: string;

  /// Height of the overview chart
  overviewHeight: number;

  overviewType: string;

  /// Shows top panel with the "Split by" selector
  showTopPanel: boolean;

  /// Show vertical line reflecting the position of the mouse-over row
  /// See also *Mouse Over Line Color*
  showMouseOverRowLine: boolean;

  /// Show vertical line reflecting the position of the current row
  /// See also *Current Line Color*
  showCurrentRowLine: boolean;

  /// Determines whether the line is highlighted when you hover over the corresponding category.
  /// Example: "Split by" = "SEX" and you hover over the "Male" category in the filter.
  showMouseOverCategory: boolean;

  showXSelector: boolean;

  showYSelectors: boolean;

  showAggrSelectors: boolean;

  /// Show the "x" close icon for each chart
  showCloseLink: boolean;

  showSplitSelector: boolean;

  showXAxis: boolean;

  showYAxis: boolean;

  invertXAxis: boolean;

  xAxisCustomTickmarks: Array<number>;

  yAxisCustomTickmarks: Array<number>;

  segmentsFont: string;

  columnSelectorsFont: string;

  markerType: string;

  markerSize: number;

  /// Show additional chart on the left
  leftPanel: string;

  histogramWidth: number;

  showMarkers: string;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IBarChartLookSettings {
  /// A categorical column to split data on (each bar represents a category)
  splitColumnName: string;

  /// Aggregation function (applicable to dates only).
  splitFunction: string;

  /// A categorical column to further split data on.
  /// Each category would become a part of the bar resulting from *Split*.
  stackColumnName: string;

  /// Value column. See *Value Aggr Type* for aggregation options.
  valueColumnName: string;

  /// Value aggregation.
  valueAggrType: string;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  colorColumnName: string;

  /// Color aggregation type.
  colorAggrType: string;

  /// Indicates whether the "no data" bar should appear
  /// when the *Split* value is not present.
  includeNulls: boolean;

  /// When true, each outermost bar is of the same width.
  /// This mode is useful for comparing relative value frequency when the *Stack* column is specified.
  relativeValues: boolean;

  /// Determines the rows shown on the scatter plot.
  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or ${WEIGHT / 2) > 100"
  filter: string;

  /// Whether to sort bars *by category* or *by value*.
  /// See also *Bar Sort Order*
  barSortType: string;

  /// Whether the bars should be sorted in ascending or descending order.
  /// See also *Bar Sort Type*.
  barSortOrder: string;

  outerMarginLeft: number;

  outerMarginRight: number;

  outerMarginTop: number;

  outerMarginBottom: number;

  backColor: number;

  axisColor: number;

  barColor: number;

  categoryColor: number;

  valueTextColor: number;

  barBorderLineMouseOverColor: number;

  barBorderLineMouseOverWidth: number;

  barBorderLineFilteredColor: number;

  barBorderLineColor: number;

  barBorderLineWidth: number;

  barCornerRadius: number;

  font: string;

  maxBarHeight: number;

  maxCategoryWidth: number;

  categoryValueWidth: number;

  minTextHeight: number;

  showValueAxis: boolean;

  showValueAxisLine: boolean;

  showValueSelector: boolean;

  showCategoryValues: boolean;

  showCategorySelector: boolean;

  showStackSelector: boolean;

  /// Whether the selected rows are indicated.
  /// Only works for cumulative aggregations such as count.
  showSelectedRows: boolean;

  /// Determines what happens when you click on a bar.
  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IDensityPlotLookSettings {
  /// Columns to be put on the X axis
  xColumnName: string;

  /// Columns to be put on the Y axis
  yColumnName: string;

  showXAxis: boolean;

  showYAxis: boolean;

  xBins: number;

  yBins: number;

  backColor: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IBoxPlotLookSettings {
  valueColumnName: string;

  categoryColumnName: string;

  /// Column to color-code boxes (Q2-Q3 region).
  /// See also *Bin Color Aggr Type*.
  binColorColumnName: string;

  /// Aggregation function for color-coding.
  /// See also *Bin Color*.
  binColorAggrType: string;

  /// Column to color-code markers.
  markerColorColumnName: string;

  /// Determines the rows shown on the scatter plot.
  // Commented until we make UI for formula editor:
  // @Prop(editor: 'editFormula')
  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or ${WEIGHT / 2) > 100"
  filter: string;

  markerType: string;

  markerSize: number;

  markerOpacity: number;

  /// Points are not shown if the number of rows is greater than *Show Values Limit*.
  showValuesLimit: number;

  /// Show points inside the Q2-Q3 bar
  showInsideValues: boolean;

  /// Show points outside Q2-Q3
  showOutsideValues: boolean;

  showValueAxis: boolean;

  showCategoryAxis: boolean;

  showLowerDash: boolean;

  showUpperDash: boolean;

  showMedianDash: boolean;

  showMeanCross: boolean;

  showStatistics: boolean;

  statistics: Array<string>;

  /// Show p-value. Press T to toggle.
  /// Currently works only when there are two categories.
  /// Welch's t-test is used for calculating the p-value.
  showPValue: boolean;

  borderColor: number;

  borderLineWidth: number;

  backColor: number;

  filteredRowsColor: number;

  filteredOutRowsColor: number;

  selectedRowsColor: number;

  missingValueColor: number;

  defaultBoxColor: number;

  showCategorySelector: boolean;

  showValueSelector: boolean;

  rowTooltip: string;

  axisFont: string;

  categoryFont: string;

  statisticsFont: string;

  axisUseColumnFormat: boolean;

  invertYAxis: boolean;

  whiskerWidthRatio: number;

  maxBinWidth: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IPieChartLookSettings {
  categoryColumnName: string;

  categoryFunction: string;

  pieSortType: string;

  pieSortOrder: string;

  includeNulls: boolean;

  segmentAngleColumnName: string;

  segmentAngleAggrType: string;

  segmentLengthColumnName: string;

  segmentLengthAggrType: string;

  /// Action to be performed when you click on a pie
  startAngle: number;

  shift: number;

  outlineLineWidth: number;

  backColor: number;

  outlineColor: number;

  mouseOverOutlineColor: number;

  innerLabelColor: number;

  showInnerPercent: boolean;

  showInnerLabel: boolean;

  showColumnSelector: boolean;

  /// Highlight part of the pie that corresponds to the mouse-over rows
  showMouseOverRowGroup: boolean;

  /// Highlight selected rows
  showSelectedRows: boolean;

  marginLeft: number;

  marginTop: number;

  marginRight: number;

  marginBottom: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IMatrixPlotLookSettings {
  /// Columns to use on the X axis
  xColumnNames: Array<string>;

  /// Column to use on the Y axis
  yColumnNames: Array<string>;

  font: string;

  cellPlotType: string;

  backColor: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IGridLookSettings {
  /// Indicates whether the control is in the grid or heatmap mode.
  /// Typically, you don't need to change it manually.
  isGrid: boolean;

  frozenColumns: number;

  missingValueColor: number;

  selectedRowsColor: number;

  selectedColsColor: number;

  currentRowColor: number;

  mouseOverRowColor: number;

  mouseOverRowStripeColor: number;

  colHeaderFont: string;

  defaultCellFont: string;

  /// Color-coding that applies to all columns.
  /// Additionally, each column can be individually color-coded.
  backColor: number;

  colHeaderTextColor: number;

  colHeaderBackColor: number;

  colHeaderMouseOverTextColor: number;

  cellTextColor: number;

  currentCellTextColor: number;

  rowHeaderBackColor: number;

  marginLeft: number;

  marginTop: number;

  marginRight: number;

  marginBottom: number;

  /// Orientation of the column header text.
  /// In spreadsheet mode, it defaults to horizontal no matter how small the columns are.
  /// In heat map mode, it depends on whether the text can fit in the area.
  /// Column header height. If not specified, it is calculated automatically.
  /// See also *Col Labels Orientation*, *Horz Col Labels Height*
  colHeaderHeight: number;

  /// Height of the column labels when the orientation is horizontal,
  /// and *Col Header Height* is not specified.
  horzColLabelsHeight: number;

  /// Height of the column labels when the orientation is vertical,
  /// and *Col Header Height* is not specified.
  vertColLabelsHeight: number;

  rowHeight: number;

  /// Indicates the way colors are sampled in the heatmap mode when there is not enough
  /// pixels on the screen for each row:
  /// True: each row is draws (but the result is blended and the resulting color might not represent any row)
  /// False: a row is sampled and then drawn as one pixel (but non-sampled rows do not get drawn at all)
  drawEveryRow: boolean;

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

  /// Reordering columns by dragging the header
  allowColReordering: boolean;

  /// Resizing columns by dragging the border between column headers
  allowColResizing: boolean;

  /// Resizing column header by dragging the border between the header and the first row
  allowColHeaderResizing: boolean;

  /// Shift+click on a header to select a column
  /// Ctrl+click to invert selection
  /// Ctrl+Shift+click to deselect
  allowColSelection: boolean;

  /// Reorder rows by drag-and-dropping
  allowRowReordering: boolean;

  /// Resize rows by dragging the border between rows on a row header
  allowRowResizing: boolean;

  /// Mouse drag on the row header selects rows
  /// Ctrl+click to invert selection
  /// Shift+mouse drag to select multiple rows
  /// Ctrl+Shift+mouse drag to unselect
  allowRowSelection: boolean;

  /// Mouse drag on the data cells selects both rows and columns
  allowBlockSelection: boolean;

  /// Whether the "hamburger" menu should be shown for a column
  /// when the mouse is over its header
  allowColumnMenu: boolean;

  /// Whether row rows can be dragged out of the grid.
  allowRowDragging: boolean;

  /// Automatically scroll current row into view when it is set from outside
  /// (for instance, as a result of clicking on a point in a scatter plot)
  autoScrollRowIntoView: boolean;

  /// Automatically scroll current row into view when this column becomes current
  autoScrollColumnIntoView: boolean;

  /// Indicates current row with the *Current Row Color*.
  showCurrentRowIndicator: boolean;

  /// Indicates mouse-over row by drawing a vertical stripe on the row header
  showMouseOverRowIndicator: boolean;

  showMouseOverRowStripe: boolean;

  showCurrentCellOutline: boolean;

  showColumnLabels: boolean;

  showRowHeader: boolean;

  showCellTooltip: boolean;

  showColumnTooltip: boolean;

  /// Include currently visible columns in a tooltip
  showVisibleColumnsInTooltip: boolean;

  showColumnGridlines: boolean;

  showRowGridlines: boolean;

  /// Whether items applicable to all viewers (such as Pickup Style) should
  /// be shown in a popup menu. Also requires *Show Context Menu*.
  showDefaultPopupMenu: boolean;

  /// Whether the context menu is shown
  showContextMenu: boolean;

  /// Whether the spreadsheet should only show rows that pass the filter
  showFilteredRowsOnly: boolean;

  /// When set to false, default menu appears under the 'Grid' submenu.
  topLevelDefaultMenu: boolean;

  extendLastColumn: boolean;

  maxFontSize: number;

  /// Whether the current object (shown in the property panel) is changed
  /// when you click on a column header.
  allowChangeCurrentObject: boolean;

  colorScheme: Array<number>;

  columnHeaderTypes: Array<string>;

  columns: Array<any>;

  isHeatmap: boolean;

  maxHeatmapColumns: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface ICalendarLookSettings {
  dateColumnName: string;

  showHeader: boolean;

  redWeekends: boolean;

  clickFiltersRows: boolean;

  showFilteredOnly: boolean;

  backColor: number;

  oddMonthColor: number;

  evenMonthColor: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface ITrellisPlotLookSettings {
  xColumnNames: Array<string>;

  yColumnNames: Array<string>;

  viewerType: string;

  categoryLabelFont: string;

  //Map innerViewerLookState;
  showControlPanel: boolean;

  syncMouseOverRow: boolean;

  backColor: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IPcPlotLookSettings {
  transformation: string;

  /// Columns to use
  columnNames: Array<string>;

  /// Columns where logarithmic axis is used.
  /// Should be a subset of *Column Names*.
  logColumnsColumnNames: Array<string>;

  colorColumnName: string;

  backColor: number;

  selectedRowsColor: number;

  filteredOutLineColor: number;

  missingValueColor: number;

  lineColor: number;

  currentLineColor: number;

  mouseOverLineColor: number;

  lineWidth: number;

  currentLineWidth: number;

  mouseOverLineWidth: number;

  showCurrentPoint: boolean;

  showMouseOverPoint: boolean;

  showMouseOverRowGroup: boolean;

  /// Whether the in-chart filters are visible
  showFilters: boolean;

  /// Whether the filtered out values are shown.
  /// See also *Filtered Out Line Color*
  showFilteredOutLines: boolean;

  showMin: boolean;

  showMax: boolean;

  minMaxHeight: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IMapViewerLookSettings {
  regionColumnName: string;

  colorColumnName: string;

  colorAggrType: string;

  map: string;

  showColorScale: boolean;

  showColorSelector: boolean;

  allowPanZoom: boolean;

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IGMapViewerLookSettings {
  latitudeColumnName: string;

  longitudeColumnName: string;

  mapOpacity: number;

  plotOpacity: number;

  plotPointerEvents: boolean;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IStatsViewerLookSettings {
  columnNames: Array<string>;

  stats: Array<string>;

  backColor: number;

  rows: string;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface ICorrelationPlotLookSettings {
  /// Columns to be put on the X axis
  xColumnNames: Array<string>;

  /// Columns to be put on the Y axis
  yColumnNames: Array<string>;

  /// Shows the Pearson correlation coefficient inside the corresponding cell.
  showPearsonR: boolean;

  /// Shows the tooltip with the corresponding scatter plot inside.
  showTooltip: boolean;

  backColor: number;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IFormLookSettings {
  /// Determines what gets shown on the form.
  syncMode: string;

  showNavigation: boolean;

  showPrevRowArrow: boolean;

  showNextRowArrow: boolean;

  showRowSelector: boolean;

  showFieldEditor: boolean;

  showDesignEditor: boolean;

  showColumnSelector: boolean;

  sketchState: Map<any, any>;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IMarkupViewerLookSettings {
  stretch: boolean;

  content: string;

  /// Whether the rendered html is passed through Grok's [Markup] engine (don't confuse it
  /// with the Markup that might be used for html rendering)
  markupEnabled: boolean;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IWordCloudLookSettings {
  wordColumnName: string;

  minSize: number;

  maxSize: number;

  backColor: number;

  autoFontSize: boolean;

  font: string;

  maxWords: number;

  sizeColumnName: string;

  sizeColumnAggrType: string;

  colorColumnName: string;

  colorColumnAggrType: string;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface INetworkDiagramLookSettings {
  node1ColumnName: string;

  node2ColumnName: string;

  edgeColorColumnName: string;

  edgeColorAggrType: string;

  edgeWidthColumnName: string;

  edgeWidthAggrType: string;

  node1SizeColumnName: string;

  node1SizeAggrType: string;

  node2SizeColumnName: string;

  node2SizeAggrType: string;

  node1ColorColumnName: string;

  node1ColorAggrType: string;

  node2ColorColumnName: string;

  node2ColorAggrType: string;

  node1ImageColumnName: string;

  node2ImageColumnName: string;

  node1LabelColumnName: string;

  node2LabelColumnName: string;

  node1Color: number;

  ///put url, or url with placeholders ... to apply individual img
  node1Img: string;

  node1Physics: boolean;

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

  nodeImg: string;

  nodeColor: number;

  edgeColor: number;

  edgeWidth: number;

  edgesPhysics: boolean;

  selectedColor: number;

  hoverColor: number;

  showColumnSelectors: boolean;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface ICardLookSettings {
  caption: string;

  /// Source-type specific value.
  value: string;

  format: string;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface ITileViewerLookSettings {
  cardMarkup: string;

  /// Whether the form auto-generates whenever columns change
  autoGenerate: boolean;

  sketchState: Map<any, any>;

  columnsJson: string;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

export interface IPivotViewerLookSettings {
  showHeader: boolean;

  columnsColumnNames: Array<string>;

  rowsColumnNames: Array<string>;

  measures: Array<string>;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

}

