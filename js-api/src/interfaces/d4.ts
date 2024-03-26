/// this file was generated automatically from d4 classes declarations


export interface IScatterPlot3dLookSettings {
  xColumnName: string;

  yColumnName: string;

  zColumnName: string;

  sizeColumnName: string;

  colorColumnName: string;

  labelColumnName: string;

  showAxes: boolean;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface ITreeMapLookSettings {
  splitByColumnNames: Array<string>;

  colorColumnName: string;

  colorAggrType: string;

  sizeColumnName: string;

  autoLayout: boolean;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IHistogramLookSettings {
  /// Whether the filtered out rows should be shown with the semi-transparent color
  /// See also *Filtered Out Color*
  showFilteredOutRows: boolean;

  /// Allows to filter table using the range slider on the bottom.
  /// This option also controls slider visibility.
  filteringEnabled: boolean;

  /// A numerical column used to calculate the distribution of values.
  valueColumnName: string;

  showXAxis: boolean;

  allowColumnSelection: boolean;

  /// Show bin selector in the left top panel when the mouse is over the histogram
  showBinSelector: boolean;

  /// Number of bins on the histogram
  bins: number;

  /// A categorical column to split data on (each bar represents a category)
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
  colorColumnName: string;

  colorAggrType: string;

  invertColorScheme: boolean;

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

  filteredOutColor: number;

  showCharts: boolean;

  marginLeft: number;

  marginTop: number;

  marginRight: number;

  marginBottom: number;

  filterMarginTop: number;

  filterMarginBottom: number;

  aggTooltipColumns: string;

  //StreamController _changes;
  allowDynamicMenus: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IFiltersLookSettings {
  active: boolean;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IScatterPlotLookSettings {
  /// Determines the rows shown on the scatter plot.
  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or ${WEIGHT / 2) > 100"
  filter: string;

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
  xColumnName: string;

  /// A column to use on the Y axis. Could be numerical or categorical.
  yColumnName: string;

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

  xAxisLabelOrientation: string;

  /// A column to be used for color-coding. Could be numerical or categorical.
  /// If not set, *Filtered Rows Color* is used for markers that pass the filter.
  /// Color palettes could defined either for columns in the column context panel,
  /// or via *Linear Color Scheme* and *Categorical Color Scheme* properties.
  colorColumnName: string;

  showColorSelector: boolean;

  invertColorScheme: boolean;

  /// A numerical column to use for size-coding markers.
  /// See also *Marker Min Size* and *Marker Max Size*.
  sizeColumnName: string;

  showSizeSelector: boolean;

  /// A categorical column that determines the shape of the markers.
  markersColumnName: string;

  markerType: string;

  // -1 default - automatic sizing based on current dataframe
  markerDefaultSize: number;

  markerOpacity: number;

  /// Randomly shift (x, y) marker position up to the *Jitter Size* pixels.
  /// Useful when multiple points fall on the same exact position.
  jitterSize: number;

  markerDrawBorder: boolean;

  markerBorderWidth: number;

  markerMinSize: number;

  markerMaxSize: number;

  /// Labels to show next to the markers.
  labelsColumnName: string;

  /// Determines the rows shown on the scatter plot.
  labelColorAsMarker: boolean;

  /// Regression line visibility (toggle by pressing R)
  showRegressionLine: boolean;

  showRegressionLineEquation: boolean;

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

  linearColorScheme: Array<number>;

  categoricalColorScheme: Array<number>;

  regressionLineColor: number;

  regressionLineTransparency: number;

  /// Determines whether the axes should follow the non-precision-related format (such as "money")
  /// set for the corresponding column.
  axesUseColumnFormat: boolean;

  formulaLines: string;

  viewport: string;

  /// Controls scatter plot tooltip visibility
  showTooltip: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface ILineChartLookSettings {
  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or ${WEIGHT / 2) > 100"
  filter: string;

  /// Deprecated, use splitColumnNames instead
  splitColumnName: string;

  /// A categorical column by which lines are split
  splitColumnNames: Array<string>;

  /// Defines a Y column for the chart on the bottom used for zooming
  overviewColumnName: string;

  /// Aggregation types for all columns
  yAggrTypes: Array<string>;

  /// When true, X axis is synchronized with the corresponding filter's range values.
  /// Otherwise, when the filter is changed points are filtered out on a chart but the min-max stays.
  axesFollowFilter: boolean;

  /// When true, multiple *Y Columns* charts get rendered on top of each other,
  /// otherwise they are stacked
  multiAxis: boolean;

  /// Column to be used on the X axis
  xColumnName: string;

  /// When defined, background is colored according to the segment column.
  /// Example: time series data with the "stimuli" column
  segmentColumnName: string;

  invertXAxis: boolean;

  showXAxis: boolean;

  showXSelector: boolean;

  xAxisLabelOrientation: string;

  xMin: number;

  xMax: number;

  /// Numerical columns to be used on Y axes.
  /// Depending on the *
  yColumnNames: Array<string>;

  showYAxis: boolean;

  yGlobalScale: boolean;

  /// Axis title to be shown on the left axis in multi-axis mode
  yAxisTitle: string;

  /// Axis title to be shown on the left axis in multi-axis mode
  y2AxisTitle: string;

  showYSelectors: boolean;

  showAggrSelectors: boolean;

  showSplitSelector: boolean;

  splineTension: number;

  markerType: string;

  markerSize: number;

  showMarkers: string;

  /// Show vertical line reflecting the position of the current row
  /// See also *Current Line Color*
  showCurrentRowLine: boolean;

  /// Determines whether the line is highlighted when you hover over the corresponding category.
  /// Example: "Split by" = "SEX" and you hover over the "Male" category in the filter.
  showMouseOverCategory: boolean;

  overviewAggrType: string;

  /// Show vertical line reflecting the position of the mouse-over row
  /// See also *Mouse Over Line Color*
  showMouseOverRowLine: boolean;

  /// Use column format for axis labels, where possible
  axesUseColumnFormat: boolean;

  /// Marker type for showing the distribution of the aggregated values
  /// when multiple values have the same X value
  whiskersType: string;

  overviewType: string;

  /// Show additional chart on the left
  leftPanel: string;

  autoLayout: boolean;

  segmentsFont: string;

  columnSelectorsFont: string;

  lineWidth: number;

  lineTransparency: number;

  /// Height of the overview chart
  overviewHeight: number;

  histogramWidth: number;

  /// If true, *X Axis Height* is calculated automatically to fit the required precision.
  /// If false, the specified *X Axis Height*
  autoAxisSize: boolean;

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

  /// Newline-separated list of column names to be used in a tooltip.
  /// Requires *showTooltip* to be enabled.
  rowTooltip: string;

  rowGroupTooltip: string;

  innerChartMarginTop: number;

  innerChartMarginBottom: number;

  outerChartMarginLeft: number;

  outerChartMarginTop: number;

  outerChartMarginRight: number;

  outerChartMarginBottom: number;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IBarChartLookSettings {
  /// Determines the rows shown on the scatter plot.
  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or ${WEIGHT / 2) > 100"
  filter: string;

  /// Determines what happens when you click on a bar.
  /// Value column. See *Value Aggr Type* for aggregation options.
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

  showValueAxis: boolean;

  showValueSelector: boolean;

  /// A categorical column to split data on (each bar represents a category)
  splitColumnName: string;

  /// Aggregation function (applicable to dates only).
  splitFunction: string;

  showCategoryValues: boolean;

  showValuesInsteadOfCategories: boolean;

  showCategorySelector: boolean;

  /// A categorical column to further split data on.
  /// Each category would become a part of the bar resulting from *Split*.
  stackColumnName: string;

  showStackSelector: boolean;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  colorColumnName: string;

  /// Color aggregation type.
  colorAggrType: string;

  invertColorScheme: boolean;

  /// Whether the selected rows are indicated.
  /// Only works for cumulative aggregations such as count.
  showSelectedRows: boolean;

  showMouseOverRect: boolean;

  autoLayout: boolean;

  maxCategoryWidth: number;

  categoryValueWidth: number;

  showValueAxisLine: boolean;

  barBorderLineMouseOverWidth: number;

  barBorderLineWidth: number;

  maxBarHeight: number;

  barCornerRadius: number;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IDensityPlotLookSettings {
  /// Columns to be put on the X axis
  xColumnName: string;

  /// Columns to be put on the Y axis
  yColumnName: string;

  autoLayout: boolean;

  axisFont: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IBoxPlotLookSettings {
  /// Determines the rows shown on the box plot.
  /// Formula that filters out rows to show.
  /// Example: "${AGE} > 20 or ${WEIGHT / 2) > 100"
  filter: string;

  showStatistics: boolean;

  categoryColumnName: string;

  showCategoryAxis: boolean;

  showCategorySelector: boolean;

  valueColumnName: string;

  invertYAxis: boolean;

  showValueAxis: boolean;

  showValueSelector: boolean;

  /// Column to color-code boxes (Q2-Q3 region).
  /// See also *Bin Color Aggr Type*.
  binColorColumnName: string;

  /// Aggregation function for color-coding.
  /// See also *Bin Color*.
  binColorAggrType: string;

  /// Column to color-code markers.
  markerColorColumnName: string;

  markerType: string;

  markerSize: number;

  markerOpacity: number;

  showMeanCross: boolean;

  showLowerDash: boolean;

  showUpperDash: boolean;

  showMedianDash: boolean;

  showValuesLimit: number;

  /// Show points inside the Q2-Q3 bar
  showInsideValues: boolean;

  /// Show points outside Q2-Q3
  showOutsideValues: boolean;

  /// Show p-value. Press T to toggle.
  /// Currently works only when there are two categories.
  /// Welch's t-test is used for calculating the p-value.
  showPValue: boolean;

  showMouseOverPoint: boolean;

  statistics: Array<string>;

  autoLayout: boolean;

  axisFont: string;

  categoryFont: string;

  statisticsFont: string;

  borderLineWidth: number;

  whiskerWidthRatio: number;

  maxBinWidth: number;

  axisUseColumnFormat: boolean;

  borderColor: number;

  backColor: number;

  filteredRowsColor: number;

  filteredOutRowsColor: number;

  selectedRowsColor: number;

  missingValueColor: number;

  defaultBoxColor: number;

  /// Controls box plot tooltip visibility
  showTooltip: string;

  /// Newline-separated list of column names to be used in a tooltip.
  /// Requires *showTooltip* to be enabled.
  rowTooltip: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IPieChartLookSettings {
  categoryColumnName: string;

  categoryFunction: string;

  pieSortType: string;

  pieSortOrder: string;

  includeNulls: boolean;

  autoLayout: boolean;

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

  aggTooltipColumns: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IMatrixPlotLookSettings {
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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface ISummaryLookSettings {
  /// Column to show aggregation on
  aggrColumnName: string;

  /// List of columns to show aggregations on
  columnNames: Array<string>;

  /// Aggregation that will be used for the columns
  aggregation: string;

  /// List of aggregations for the columns (temporarily unavailable from UI)
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
  colorColumnName: string;

  /// Color aggregation type.
  colorAggrType: string;

  /// Whether to apply color coding to the background or to the text.
  applyTo: string;

  /// Custom color scheme for the color-coding.
  colorSchemes: Array<Array<number>>;

  invertColorScheme: boolean;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface ISparklinesLookSettings {
  /// List of columns to show aggregations on
  columnNames: Array<string>;

  /// Aggregation that will be used for the columns
  aggregation: string;

  /// List of aggregations for the columns
  aggregations: Array<string>;

  sparklineType: string;

  /// Numerical column to be used for color-coding.
  /// The values in the bin get aggregated using the *Color Aggr Type* property.
  colorColumnName: string;

  /// Color aggregation type.
  colorAggrType: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IGridLookSettings {
  /// Determines the rows shown in grid
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

  /// Mouse drag on the rows headers selects rows
  /// Ctrl+click to invert selection
  /// Shift+mouse drag to select multiple rows
  /// Ctrl+Shift+mouse drag to unselect
  allowRowSelection: boolean;

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

  showColumnGridlines: boolean;

  /// Reordering columns by dragging the header
  allowColReordering: boolean;

  /// Whether the current object (shown in the context panel) is changed
  /// when you click on a column header.
  allowChangeCurrentObject: boolean;

  /// Whether row (rows) can be dragged out of the grid.
  allowRowDragging: boolean;

  extendLastColumn: boolean;

  /// Resize rows by dragging the border between rows on a row header
  allowRowResizing: boolean;

  /// Indicates the way colors are sampled in the heatmap mode when there is not enough
  /// pixels on the screen for each row:
  /// True: each row is draws (but the result is blended and the resulting color might not represent any row)
  /// False: a row is sampled and then drawn as one pixel (but non-sampled rows do not get drawn at all)
  drawEveryRow: boolean;

  /// Whether the context menu is shown
  showContextMenu: boolean;

  /// Whether to show notifications when the user tries to edit a read-only table
  showReadOnlyNotifications: boolean;

  frozenColumns: number;

  showCurrentCellOutline: boolean;

  /// Color-coding that applies to all columns.
  /// Additionally, each column can be individually color-coded.
  defaultCellFont: string;

  maxFontSize: number;

  colHeaderFont: string;

  /// Orientation of the column header text.
  /// In spreadsheet mode, it defaults to horizontal no matter how small the columns are.
  /// In heat map mode, it depends on whether the text can fit in the area.
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
  globalColorScaling: boolean;

  /// Controls grid tooltip visibility
  showTooltip: string;

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

  /// Determines whether newly added columns are added to the grid
  syncNewColumns: boolean;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface ITrellisPlotLookSettings {
  xColumnNames: Array<string>;

  yColumnNames: Array<string>;

  viewerType: string;

  filter: string;

  categoryLabelFont: string;

  innerViewerLooks: Map<string, any>;

  globalScale: boolean;

  showGridlines: string;

  showXSelectors: boolean;

  showYSelectors: boolean;

  showXLabels: boolean;

  showYLabels: boolean;

  showControlPanel: boolean;

  syncMouseOverRow: boolean;

  packCategories: boolean;

  useTiledView: boolean;

  tilesPerRow: number;

  autoLayout: boolean;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IPcPlotLookSettings {
  /// Whether the filtered out values are shown.
  /// See also *Filtered Out Line Color*
  showFilteredOutLines: boolean;

  /// Columns to use
  columnNames: Array<string>;

  /// Columns where logarithmic axis is used.
  /// Should be a subset of *Column Names*.
  logColumnsColumnNames: Array<string>;

  colorColumnName: string;

  /// Determines the way a value is mapped to the vertical scale.
  /// TRUE: bottom is column minimum, top is column maximum. Use when columns contain values in different units
  /// FALSE: uses the same scale. This lets you compare values across columns
  /// if units are the same (for instance, use it for tracking change over time).'
  normalizeEachColumn: boolean;

  showCurrentLine: boolean;

  showMouseOverLine: boolean;

  showMouseOverRowGroup: boolean;

  showMin: boolean;

  showMax: boolean;

  /// Whether the in-chart filters are visible
  showFilters: boolean;

  currentLineWidth: number;

  autoLayout: boolean;

  lineWidth: number;

  mouseOverLineWidth: number;

  minMaxHeight: number;

  transformation: string;

  backColor: number;

  selectedRowsColor: number;

  filteredOutLineColor: number;

  missingValueColor: number;

  lineColor: number;

  currentLineColor: number;

  mouseOverLineColor: number;

  showDensity: boolean;

  showMinMax: boolean;

  showLabels: boolean;

  maxCategories: number;

  horzMargin: number;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

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

  showSaveFile: boolean;

  showOpenFile: boolean;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IMarkupViewerLookSettings {
  stretch: boolean;

  content: string;

  /// Whether the rendered html is passed through Grok's [Markup] engine (don't confuse it
  /// with the Markup that might be used for html rendering)
  markupEnabled: boolean;

  //StreamController _changes;
  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

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

  autoLayout: boolean;

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

  selectRowsOnClick: boolean;

  selectEdgesOnClick: boolean;

  /// Loads external data; called when you double-click on a node.
  /// The specified function gets called with the node value as a single argument.
  /// Its signature: `dataframe expand(dynamic nodeId)`.
  onNodeExpandFunction: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface ITileViewerLookSettings {
  lanesColumnName: string;

  cardMarkup: string;

  allowDragBetweenLanes: boolean;

  /// Whether the form auto-generates whenever columns change
  autoGenerate: boolean;

  sketchState: Map<any, any>;

  columnsJson: string;

  lanes: Array<string>;

  //StreamController _changes;
  //StreamController _changes;
  allowDynamicMenus: boolean;

  // Properties common for all viewers
  // todo: use code generation
  // Properties common for all viewers
  // todo: use code generation
  showContextMenu: boolean;

  title: string;

  showTitle: boolean;

  table: string;

  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  // Viewer description that gets shown at the *Descriptor Position*.
  // Markup is supported.
  description: string;

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

export interface IPivotViewerLookSettings {
  /// Determines the rows used in pivot viewer.
  showHeader: boolean;

  pivotColumnNames: Array<string>;

  groupByColumnNames: Array<string>;

  aggregateColumnNames: Array<string>;

  aggregateAggTypes: Array<string>;

  viewerSettings: Array<any>;

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

  // Help to be shown when user clicks on the '?' icon on top.
  // Could either be in markdown, or a URL (starting with '/' or 'http').
  help: string;

}

