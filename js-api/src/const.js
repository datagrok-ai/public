/**
 * @typedef {string} AggregationType
 * @typedef {string} SyncType
 * @typedef {string} JoinType
 * @typedef {string} ColumnType
 * @typedef {string} ViewerType
 * @typedef {string} ObjectType
 * @typedef {string} ViewerPropertyType
 * @typedef {string} Type
 * @typedef {string} SemType
 * @typedef {string} SimilarityMetric
 * @typedef {string} ColorType
 * @typedef {string} DockType
 **/
/** @enum {AggregationType} */
export var AGG;
(function (AGG) {
  AGG['KEY'] = 'key';
  AGG['PIVOT'] = 'pivot';
  AGG['FIRST'] = 'first';
  AGG['TOTAL_COUNT'] = 'count';
  AGG['VALUE_COUNT'] = 'values';
  AGG['UNIQUE_COUNT'] = 'unique';
  AGG['MISSING_VALUE_COUNT'] = 'nulls';
  AGG['MIN'] = 'min';
  AGG['MAX'] = 'max';
  AGG['SUM'] = 'sum';
  AGG['MED'] = 'med';
  AGG['AVG'] = 'avg';
  AGG['STDEV'] = 'stdev';
  AGG['VARIANCE'] = 'variance';
  AGG['SKEW'] = 'skew';
  AGG['KURT'] = 'kurt';
  AGG['Q1'] = 'q1';
  AGG['Q2'] = 'q2';
  AGG['Q3'] = 'q3';
  AGG['SELECTED_ROWS_COUNT'] = '#selected';
})(AGG || (AGG = {}));
/** @enum {SyncType} */
export var SYNC_TYPE;
(function (SYNC_TYPE) {
  SYNC_TYPE['CURRENT_ROW_TO_ROW'] = 'row to row';
  SYNC_TYPE['CURRENT_ROW_TO_SELECTION'] = 'row to selection';
  SYNC_TYPE['CURRENT_ROW_TO_FILTER'] = 'row to filter';
  SYNC_TYPE['MOUSE_OVER_ROW_TO_SELECTION'] = 'mouse-over to selection';
  SYNC_TYPE['MOUSE_OVER_ROW_TO_FILTER'] = 'mouse-over to filter';
  SYNC_TYPE['FILTER_TO_FILTER'] = 'filter to filter';
  SYNC_TYPE['FILTER_TO_SELECTION'] = 'filter to selection';
  SYNC_TYPE['SELECTION_TO_FILTER'] = 'selection to filter';
  SYNC_TYPE['SELECTION_TO_SELECTION'] = 'selection to selection';
})(SYNC_TYPE || (SYNC_TYPE = {}));
export const INT_NULL = -2147483648;
export const FLOAT_NULL = 2.6789344063684636e-34;
/** @enum {JoinType} */
export var JOIN_TYPE;
(function (JOIN_TYPE) {
  JOIN_TYPE['INNER'] = 'inner';
  JOIN_TYPE['OUTER'] = 'outer';
  JOIN_TYPE['LEFT'] = 'left';
  JOIN_TYPE['RIGHT'] = 'right';
})(JOIN_TYPE || (JOIN_TYPE = {}));
/** @enum {ColumnType} */
export var COLUMN_TYPE;
(function (COLUMN_TYPE) {
  COLUMN_TYPE['STRING'] = 'string';
  COLUMN_TYPE['INT'] = 'int';
  COLUMN_TYPE['FLOAT'] = 'double';
  COLUMN_TYPE['BOOL'] = 'bool';
  COLUMN_TYPE['DATE_TIME'] = 'datetime';
  COLUMN_TYPE['BIG_INT'] = 'bigint';
  COLUMN_TYPE['QNUM'] = 'qnum';
})(COLUMN_TYPE || (COLUMN_TYPE = {}));
/** @enum {Type} */
export var TYPE;
(function (TYPE) {
  TYPE['INT'] = 'int';
  TYPE['BIG_INT'] = 'bigint';
  TYPE['FLOAT'] = 'double';
  TYPE['NUM'] = 'num';
  TYPE['QNUM'] = 'qnum';
  TYPE['BOOL'] = 'bool';
  TYPE['STRING'] = 'string';
  TYPE['STRING_LIST'] = 'string_list';
  TYPE['DATE_TIME'] = 'datetime';
  TYPE['OBJECT'] = 'object';
  TYPE['DATA_FRAME'] = 'dataframe';
  TYPE['DATA_FRAME_LIST'] = 'dataframe_list';
  TYPE['CELL'] = 'cell';
  TYPE['COLUMN'] = 'column';
  TYPE['COLUMN_LIST'] = 'column_list';
  TYPE['GRAPHICS'] = 'graphics';
  TYPE['ROW_FILTER'] = 'tablerowfiltercall';
  TYPE['COLUMN_FILTER'] = 'colfiltercall';
  TYPE['BIT_SET'] = 'bitset';
  TYPE['MAP'] = 'map';
  TYPE['DYNAMIC'] = 'dynamic';
  TYPE['VIEWER'] = 'viewer';
  TYPE['LIST'] = 'list';
  TYPE['SEM_VALUE'] = 'semantic_value';
  TYPE['FUNC_CALL'] = 'funccall';
  TYPE['PROPERTY'] = 'property';
  TYPE['CATEGORICAL'] = 'categorical';
  TYPE['NUMERICAL'] = 'numerical';
  TYPE['GRID_CELL_RENDER_ARGS'] = 'GridCellRenderArgs';
  TYPE['VIEW'] = 'view';
  TYPE['TABLE_VIEW'] = 'TableView';
  TYPE['USER'] = 'User';
  TYPE['MENU'] = 'Menu';
  TYPE['PROJECT'] = 'Project';
  TYPE['SEMANTIC_VALUE'] = 'semantic_value';
  TYPE['EVENT_DATA'] = 'event_data';
  TYPE['PROGRESS_INDICATOR'] = 'progressindicator';
  TYPE['CREDENTIALS'] = 'Credentials';
  TYPE['SCRIPT_ENVIRONMENT'] = 'ScriptEnvironment';
  TYPE['NOTEBOOK'] = 'Notebook';
})(TYPE || (TYPE = {}));
export const TYPES_SCALAR = new Set([TYPE.INT, TYPE.BIG_INT, TYPE.FLOAT, TYPE.NUM, TYPE.BOOL, TYPE.STRING]);
/** @enum {ViewerPropertyType} */
export var VIEWER_PROPERTY_TYPE;
(function (VIEWER_PROPERTY_TYPE) {
  VIEWER_PROPERTY_TYPE['STRING'] = 'string';
  VIEWER_PROPERTY_TYPE['INT'] = 'int';
  VIEWER_PROPERTY_TYPE['FLOAT'] = 'double';
  VIEWER_PROPERTY_TYPE['BOOL'] = 'bool';
  VIEWER_PROPERTY_TYPE['DATE_TIME'] = 'datetime';
  VIEWER_PROPERTY_TYPE['BIG_INT'] = 'bigint';
  VIEWER_PROPERTY_TYPE['COLUMN'] = 'column';
  VIEWER_PROPERTY_TYPE['COLUMN_LIST'] = 'column_list';
  VIEWER_PROPERTY_TYPE['DATA_FRAME'] = 'dataframe';
})(VIEWER_PROPERTY_TYPE || (VIEWER_PROPERTY_TYPE = {}));
export var VIEW_TYPE;
(function (VIEW_TYPE) {
  VIEW_TYPE['TABLE_VIEW'] = 'TableView';
})(VIEW_TYPE || (VIEW_TYPE = {}));
;
///////
/** @enum {SemType} */
export const SEMTYPE = {
  EMAIL: 'Email Address',
  URL: 'URL',
  PHONE_NUMBER: 'Phone Number',
  CITY: 'City',
  COUNTRY: 'Country',
  GENDER: 'Gender',
  STATE: 'State',
  COUNTY: 'County',
  PLACE_NAME: 'Place Name',
  ZIP_CODE: 'Zip Code',
  AREA_CODE: 'Area Code',
  STREET_ADDRESS: 'Street Address',
  TEXT: 'Text',
  DURATION: 'Duration',
  LATITUDE: 'Latitude',
  LONGITUDE: 'Longitude',
  IP_ADDRESS: 'IP Address',
  MOLECULE: 'Molecule',
  SUBSTRUCTURE: 'Substructure',
  MONEY: 'Money',
  IMAGE: 'Image',
  FILE: 'File',
};
/////// Stats
/** @enum {STATS} */
export const STATS = {
  TOTAL_COUNT: 'count',
  VALUE_COUNT: 'values',
  UNIQUE_COUNT: 'unique',
  MISSING_VALUE_COUNT: 'nulls',
  MIN: 'min',
  MAX: 'max',
  SUM: 'sum',
  MED: 'med',
  AVG: 'avg',
  STDEV: 'stdev',
  VARIANCE: 'variance',
  SKEW: 'skew',
  KURT: 'kurt',
  Q1: 'q1',
  Q2: 'q2',
  Q3: 'q3'
};
/////// Tags
/** @enum {TAGS} */
export const TAGS = {
  LAYOUT_ID: 'layout-id',
  DESCRIPTION: 'description',
  TOOLTIP: '.tooltip',
  /** JSON-encoded list of strings to be used in a cell editor. Applicable for string columns only. */
  CHOICES: '.choices',
  /** When set to 'true', switches the cell editor to a combo box that only allows to choose values
     from a list of already existing values in the column.
     Applicable for string columns only.
     See also [TAGS_CHOICES]. */
  AUTO_CHOICES: '.auto-choices',
  ID: '.id',
  COLOR_CODING_CONDITIONAL: '.color-coding-conditional'
};
////// Viewers
/** @enum {VIEWER} */
export var VIEWER;
(function (VIEWER) {
  VIEWER['HISTOGRAM'] = 'Histogram';
  VIEWER['BAR_CHART'] = 'Bar chart';
  VIEWER['BOX_PLOT'] = 'Box plot';
  VIEWER['CALENDAR'] = 'Calendar';
  VIEWER['CORR_PLOT'] = 'Correlation plot';
  VIEWER['DENSITY_PLOT'] = 'Density plot';
  VIEWER['FILTERS'] = 'Filters';
  VIEWER['FORM'] = 'Form';
  VIEWER['GLOBE'] = 'Globe';
  VIEWER['GRID'] = 'Grid';
  VIEWER['GOOGLE_MAP'] = 'Google map';
  VIEWER['HEAT_MAP'] = 'Heat map';
  VIEWER['LINE_CHART'] = 'Line chart';
  VIEWER['SHAPE_MAP'] = 'Shape map';
  VIEWER['MARKUP'] = 'Markup';
  VIEWER['MATRIX_PLOT'] = 'Matrix plot';
  VIEWER['NETWORK_DIAGRAM'] = 'Network diagram';
  VIEWER['PC_PLOT'] = 'PC Plot';
  VIEWER['PIE_CHART'] = 'Pie chart';
  VIEWER['SCATTER_PLOT'] = 'scatter plot';
  VIEWER['SCATTER_PLOT_3D'] = '3d scatter plot';
  VIEWER['STATISTICS'] = 'Statistics';
  VIEWER['TILE_VIEWER'] = 'Tile viewer';
  VIEWER['TREE_MAP'] = 'Tree map';
  VIEWER['TRELLIS_PLOT'] = 'Trellis plot';
  VIEWER['WORD_CLOUD'] = 'Word cloud';
})(VIEWER || (VIEWER = {}));
/** @enum {SIMILARITY_METRIC} */
export var SIMILARITY_METRIC;
(function (SIMILARITY_METRIC) {
  SIMILARITY_METRIC['TANIMOTO'] = 'tanimoto';
  SIMILARITY_METRIC['DICE'] = 'dice';
  SIMILARITY_METRIC['COSINE'] = 'cosine';
  SIMILARITY_METRIC['SOKAL'] = 'sokal';
  SIMILARITY_METRIC['RUSSEL'] = 'russel';
  SIMILARITY_METRIC['ROGOT_GOLDBERG'] = 'rogot-goldberg';
  SIMILARITY_METRIC['KULCZYNSKI'] = 'kulczynski';
  SIMILARITY_METRIC['MC_CONNAUGHEY'] = 'mc-connaughey';
  SIMILARITY_METRIC['ASYMMETRIC'] = 'asymmetric';
  SIMILARITY_METRIC['BRAUN_BLANQUET'] = 'braun-blanquet';
})(SIMILARITY_METRIC || (SIMILARITY_METRIC = {}));
/** @enum {DemoDatasetName} */
export var DEMO_DATASET;
(function (DEMO_DATASET) {
  DEMO_DATASET['WELLS'] = 'wells';
  DEMO_DATASET['DEMOG'] = 'demog';
  DEMO_DATASET['BIOSENSOR'] = 'biosensor';
  DEMO_DATASET['RANDOM_WALK'] = 'random walk';
  DEMO_DATASET['GEO'] = 'geo';
})(DEMO_DATASET || (DEMO_DATASET = {}));
/** @enum {DockType} */
export var DOCK_TYPE;
(function (DOCK_TYPE) {
  DOCK_TYPE['LEFT'] = 'left';
  DOCK_TYPE['RIGHT'] = 'right';
  DOCK_TYPE['TOP'] = 'top';
  DOCK_TYPE['DOWN'] = 'down';
  DOCK_TYPE['FILL'] = 'fill';
})(DOCK_TYPE || (DOCK_TYPE = {}));
