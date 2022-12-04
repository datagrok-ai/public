
/** @enum {AGG} */
export enum AGG {
  KEY = "key",      // Special case= to be used in the 'group by' statement.
  PIVOT = "pivot",  // Special case= to be used as a pivot.
  FIRST = "first",
  TOTAL_COUNT = "count",
  VALUE_COUNT = "values",
  UNIQUE_COUNT = "unique",
  MISSING_VALUE_COUNT = "nulls",
  MIN = "min",
  MAX = "max",
  SUM = "sum",
  MED = "med",
  AVG = "avg",
  STDEV = "stdev",
  VARIANCE = "variance",
  SKEW = "skew",
  KURT = "kurt",
  Q1 = "q1",
  Q2 = "q2",
  Q3 = "q3",
  SELECTED_ROWS_COUNT = "#selected"
}

export enum STR_AGG {
  CONCAT_ALL = 'concat all',
  CONCAT_UNIQUE = 'concat unique',
  LONGEST = 'longest',
  SHORTEST = 'shortest',
  MOST_FREQUENT = 'most frequent',
  CONCAT_COUNTS = 'concat counts',
}

export enum STAT_COUNTS {
  TOTAL_COUNT = 'count',
  VALUE_COUNT = 'values',
  UNIQUE_COUNT = 'unique',
  MISSING_VALUE_COUNT = 'nulls',
}

/** @enum {SYNC_TYPE} */
export enum SYNC_TYPE {
  CURRENT_ROW_TO_ROW = 'row to row',
  CURRENT_ROW_TO_SELECTION = 'row to selection',
  CURRENT_ROW_TO_FILTER = 'row to filter',
  MOUSE_OVER_ROW_TO_SELECTION = 'mouse-over to selection',
  MOUSE_OVER_ROW_TO_FILTER = 'mouse-over to filter',
  FILTER_TO_FILTER = 'filter to filter',
  FILTER_TO_SELECTION = 'filter to selection',
  SELECTION_TO_FILTER = 'selection to filter',
  SELECTION_TO_SELECTION = 'selection to selection',
}

export const INT_NULL = -2147483648;
export const FLOAT_NULL = 2.6789344063684636e-34;

/** @enum {JoinType} */
export enum JOIN_TYPE {
  INNER = 'inner',
  OUTER = 'outer',
  LEFT = 'left',
  RIGHT = 'right'
}

/** @enum {COLUMN_TYPE} */
export enum COLUMN_TYPE {
  STRING = 'string',
  INT = 'int',
  FLOAT = 'double',
  BOOL = 'bool',
  BYTE_ARRAY = 'byte_array',
  DATE_TIME = 'datetime',
  BIG_INT = 'bigint',
  QNUM = 'qnum',
  DATA_FRAME = 'dataframe',
  OBJECT = 'object',
}


/** @enum {TYPE} */
export enum TYPE {
  INT = 'int',
  BIG_INT = 'bigint',
  FLOAT = 'double',
  NUM = 'num',
  QNUM = 'qnum',
  BOOL = 'bool',
  STRING = 'string',
  STRING_LIST = 'string_list',
  DATE_TIME = 'datetime',
  OBJECT = 'object',
  BYTE_ARRAY = 'byte_array',
  DATA_FRAME = 'dataframe',
  DATA_FRAME_LIST = 'dataframe_list',
  CELL = 'cell',
  COLUMN = 'column',
  COLUMN_LIST = 'column_list',
  GRAPHICS = 'graphics',
  ROW_FILTER = 'tablerowfiltercall',
  COLUMN_FILTER = 'colfiltercall',
  BIT_SET = 'bitset',
  MAP = 'map',
  DYNAMIC = 'dynamic',
  VIEWER = 'viewer',  // [ViewerBase] subclasses
  LIST = 'list',
  SEM_VALUE = 'semantic_value',
  FUNC = 'func',
  FUNC_CALL = 'funccall',
  PROPERTY = 'property',
  CATEGORICAL = 'categorical',
  NUMERICAL = 'numerical',
  GRID_CELL_RENDER_ARGS = 'GridCellRenderArgs',

  VIEW = 'view',
  TABLE_VIEW = 'TableView',
  USER = 'User',
  MENU = 'Menu',
  PROJECT = 'Project',
  SEMANTIC_VALUE = 'semantic_value',
  EVENT_DATA = 'event_data',
  PROGRESS_INDICATOR = 'progressindicator',
  CREDENTIALS = 'Credentials',
  SCRIPT_ENVIRONMENT = 'ScriptEnvironment',
  NOTEBOOK = 'Notebook'
}

// export type FILTER_TYPE =
//   'histogram' | 'categorical' | 'multi-value' | 'bool-columns' |
//   'free-text' | 'column-free-text' | 'Chem:substructureFilter';

export enum FILTER_TYPE {
  HISTOGRAM = 'histogram',
  CATEGORICAL = 'categorical',
  MULTI_VALUE = 'multi-value',
  BOOL_COLUMNS = 'bool-columns',
  FREE_TEXT = 'free-text',
  COLUMN_FREE_TEXT = 'column-free-text',
  SUBSTRUCTURE = 'Chem:substructureFilter'
}

export const TYPES_SCALAR = new Set([TYPE.INT, TYPE.BIG_INT, TYPE.FLOAT, TYPE.NUM, TYPE.BOOL, TYPE.STRING]);

/** @enum {VIEWER_PROPERTY_TYPE} */
export enum VIEWER_PROPERTY_TYPE {
  STRING = 'string',
  INT = 'int',
  FLOAT = 'double',
  BOOL = 'bool',
  DATE_TIME = 'datetime',
  BIG_INT = 'bigint',

  COLUMN = 'column',
  COLUMN_LIST = 'column_list',
  DATA_FRAME = 'dataframe'
}

export enum VIEW_TYPE {
  TABLE_VIEW = 'TableView'
}

///////
/** @enum {SEMTYPE} */
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
  MACROMOLECULE: 'Macromolecule',
  HELM: 'HELM',
  SUBSTRUCTURE: 'Substructure',
  MONEY: 'Money',
  IMAGE: 'Image',
  FILE: 'File',
}

export const UNITS = {
  Molecule : {
    SMILES: 'smiles',
    MOLBLOCK: 'molblock',
    V3K_MOLBLOCK: 'v3Kmolblock',
    INCHI: 'inchi'
  }
}

/////// Stats

/** @enum {STATS} */
export const STATS = {
  TOTAL_COUNT: "count",
  VALUE_COUNT: "values",
  UNIQUE_COUNT: "unique",
  MISSING_VALUE_COUNT: "nulls",
  MIN: "min",
  MAX: "max",
  SUM: "sum",
  MED: "med",
  AVG: "avg",
  STDEV: "stdev",
  VARIANCE: "variance",
  SKEW: "skew",
  KURT: "kurt",
  Q1: "q1",
  Q2: "q2",
  Q3: "q3"
}

/////// Tags

/** @enum {TAGS} */
export const TAGS = {
  LAYOUT_ID: 'layout-id',  // when set in a column tag, it gets used for layout column matching
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

  COLOR_CODING_TYPE: '.color-coding-type',
  COLOR_CODING_CONDITIONAL: '.color-coding-conditional',
  COLOR_CODING_CATEGORICAL: '.color-coding-categorical',
  COLOR_CODING_LINEAR: '.color-coding-linear',

  MARKER_CODING: '.marker-coding',
  FORMULA_LINES: '.formula-lines',

  MULTI_VALUE_SEPARATOR: '.multi-value-separator',
  /** When a dataframe is loaded from a CSV, the maximum number of significant digits
   in the fractional part for each numeric column is determined  */
  SOURCE_PRECISION: '.source-precision',
  /** Set on a dataframe column; used to format column contents in grids, CSV export, passing to scripts */
  FORMAT: 'format',
  FORMULA: 'formula',
  SEMTYPE: 'quality',

  IGNORE_CUSTOM_FILTER: '.ignore-custom-filter',
  STRUCTURE_FILTER_TYPE: '.structure-filter-type',

  CELL_RENDERER: 'cell.renderer',
  UNITS: 'units',  // see DG.UNITS

  FRIENDLY_NAME: 'friendlyName',
  ALLOW_RENAME: '.allow-rename',

  CHEM: {
    SCAFFOLD: 'chem-scaffold'
  }
}

export const FUNC_TYPES = {
  /** An application that gets shown in the app store.
    * Signature: app() */
  APP: 'app',

  /** Context-specific widget that appears on the property panel
    * Signature: panel(x: any): Widget */
  PANEL: 'panel',

  /** Gets invoked when the containing package is initialized
    * Signature: init() */
  INIT: 'init',

  /** Gets invoked at platform startup. Use it wisely as the whole package will get initialized.
    * Signature: autostart() */
  AUTOSTART: 'autostart',

  /** Semantic type detector for a column. Gets invoked when a new dataframe is imported into the platform.
   *  Implementation should either set column.semType directly, or return the semantic type that will get assigned.
   *  Signature: semTypeDetector(Column): string */
  SEM_TYPE_DETECTOR: 'semTypeDetector',

  /** Creates a viewer (or editor) for a file with the specified extension.
   *  The extension is derived from the `fileViewer-<extension>` tag.
   *  Used in the file system browser.
   *  Signature: fileViewer(FileInfo): View */
  FILE_VIEWER: 'fileViewer',

  /** Exports a file. Gets added to the "export" menu at startup.
   *  Signature: fileExporter() */
  FILE_EXPORTER: 'fileExporter',

  /** Handles custom file formats.
   * The `meta.ext` parameter should contain a comma-separated list of extensions.
   * Signature: fileImporter(x: string | TypedArray): DataFrame[] */
  FILE_IMPORTER: 'file-handler',

  /** Creates a cell renderer that is used for rendering cells for specific semantic types.
   *  Semantic type is derived from the `cellRenderer-<semType>` tag.
   *  Signature: cellRenderer(): GridCellRenderer */
  CELL_RENDERER: 'cellRenderer',

  /** Edits package settings.
   *  Signature: packageSettingsEditor(): Widget */
  PACKAGE_SETTINGS_EDITOR: 'packageSettingsEditor',

  /** Makes a widget appear on the welcome screen
   *  Signature: dashboard(): DG.Widget */
  DASHBOARD: 'dashboard',

  /**
   * Function analysis. Examples: sensitivity analysis, parameter editor
   * Func => View */
  FUNCTION_ANALYSIS: 'functionAnalysis',

  /** Converts values. Has one input and one output */
  CONVERTER: 'converter'
}

export enum LOG_LEVEL {
  DEBUG = 'debug',
  INFO = 'info',
  WARNING = 'warning',
  ERROR = 'error',
  AUDIT = 'audit',
  USAGE=  'usage'
}

////// Viewers
/** @enum {VIEWER} */
export enum VIEWER {
  HISTOGRAM = 'Histogram',
  BAR_CHART = 'Bar chart',
  BOX_PLOT = 'Box plot',
  CALENDAR = 'Calendar',
  CORR_PLOT = 'Correlation plot',
  DENSITY_PLOT = 'Density plot',
  FILTERS = 'Filters',
  FORM = 'Form',
  GLOBE = 'Globe',
  GRID = 'Grid',
  GOOGLE_MAP = 'Google map',
  HEAT_MAP = 'Heat map',
  LINE_CHART = 'Line chart',
  SHAPE_MAP = 'Shape Map',
  MARKUP = 'Markup',
  MATRIX_PLOT = 'Matrix plot',
  NETWORK_DIAGRAM = 'Network diagram',
  PC_PLOT = 'PC Plot',
  PIE_CHART = 'Pie chart',
  SCATTER_PLOT = 'Scatter plot',
  SCATTER_PLOT_3D = '3d scatter plot',
  STATISTICS = 'Statistics',
  TILE_VIEWER = 'Tile Viewer',
  TREE_MAP = 'Tree map',
  TRELLIS_PLOT = 'Trellis plot',
  WORD_CLOUD = 'Word cloud',
  TIMELINES = 'TimelinesViewer',
  SURFACE_PLOT = 'SurfacePlot'
}

/** @enum {LINE_CHART_SERIES_TYPE} */
export enum LINE_CHART_SERIES_TYPE {
  LINE = 'Line Chart',
  AREA = 'Area Chart',
  BAR = 'Bar Chart',
  STACKED_BAR = 'Stacked Bar Chart',
  STACKED_AREA = 'Stacked Area Chart'
}

/** @enum {SIMILARITY_METRIC} */
export enum SIMILARITY_METRIC {
  TANIMOTO = 'tanimoto',
  DICE = 'dice',
  COSINE = 'cosine',
  SOKAL = 'sokal',
  RUSSEL = 'russel',
  ROGOT_GOLDBERG = 'rogot-goldberg',
  KULCZYNSKI = 'kulczynski',
  MC_CONNAUGHEY = 'mc-connaughey',
  ASYMMETRIC = 'asymmetric',
  BRAUN_BLANQUET = 'braun-blanquet'
}

/** @enum {DEMO_DATASET} */
export enum DEMO_DATASET {
  WELLS = 'wells',
  DEMOG = 'demog',
  BIOSENSOR = 'biosensor',
  RANDOM_WALK = 'random walk',
  GEO = 'geo',
  MOLECULES = 'molecules',
  DOSE_RESPONSE = 'dose-response',
}

/** @enum {DOCK_TYPE} */
export enum DOCK_TYPE {
  LEFT = "left",
  RIGHT = "right",
  TOP = "up",
  DOWN = "down",
  FILL = "fill",
}

export enum LEGEND_POSITION {
  LEFT = "left",
  RIGHT = "right",
  TOP = "top",
  BOTTOM = "bottom",
}

export enum COLOR_CODING_TYPE {
  CATEGORICAL = 'Categorical',
  CONDITIONAL = 'Conditional',
  LINEAR = 'Linear',
  OFF = 'Off',
}

export enum SCRIPT_LANGUAGE {
  JAVASCRIPT = 'javascript',
  JULIA = 'julia',
  OCTAVE = 'octave',
  PYTHON = 'python',
  R = 'r',
  NODE = 'nodejs',
  GROK = 'grok',
}

export enum NAMED_VALIDATORS {
  CONTAINS_MISSING_VALUES = 'containsMissingValues',
  COLUMN_NAME = 'columnName',
  COLUMN_IS_NUMERICAL = 'columnIsNumerical',
  COLUMN_IS_CATEGORICAL = 'columnIsCategorical',
  NOT_EMPTY = 'notEmpty',
}

export enum MARKER_TYPE {
  CIRCLE = "circle",
  CIRCLE_BORDER = "circle border",
  SQUARE = "square",
  SQUARE_BORDER = "square border",
  CROSS_BORDER = "cross border",
  CROSS_X_BORDER = "cross x border",
  OUTLIER = "outlier",
  DIAMOND = "diamond",
  DIAMOND_BORDER = "diamond border",
  TRIANGLE_TOP = "triangle top",
  TRIANGLE_RIGHT = "triangle right",
  TRIANGLE_BOTTOM = "triangle bottom",
  TRIANGLE_LEFT = "triangle left",
  TRIANGLE_TOP_BORDER = "triangle top border",
  TRIANGLE_RIGHT_BORDER = "triangle right border",
  TRIANGLE_BOTTOM_BORDER = "triangle bottom border",
  TRIANGLE_LEFT_BORDER = "triangle left border",
  ASTERISK = "asterisk",
  STAR = "star",
  DOT = "dot",
  GRADIENT = "gradient",
}

export enum USER_STATUS {
  STATUS_NEW = "new",
  STATUS_ACTIVE = "active",
  STATUS_BLOCKED = "blocked",
  STATUS_GUEST = "guest"
}

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
 * 
 * @typedef {Object} ElementOptions
 * @property {string} id
 * @property {string} classes
 * @property {Object} style
 * 
 * @typedef {Object} CsvImportOptions
 * @property {string} delimiter
 * @property {string} decimalSeparator
 * @property {string} thousandSeparator
 * 
 * @typedef {function(number): boolean} IndexPredicate
 * @typedef {function(String): boolean} StringPredicate
 **/

export type AggregationType = `${AGG}`;
export type ColumnAggregationType = `${AGG}` | `${STR_AGG}` | string;
export type SyncType = `${SYNC_TYPE}`;
export type JoinType = `${JOIN_TYPE}`;
export type ColumnType = `${COLUMN_TYPE}`;
export type ViewerType = `${VIEWER}` | string;
export type ObjectType = string;
export type ViewerPropertyType = string;
export type Type = `${TYPE}`;
export type SemType = string;
export type SimilarityMetric = `${SIMILARITY_METRIC}`;
export type ColorType = number;
export type ColorCodingType = `${COLOR_CODING_TYPE}`;
export type MarkerCodingType = `${MARKER_TYPE}`;
export type DemoDatasetName = `${DEMO_DATASET}`;
export type DockType = `${DOCK_TYPE}`;
export type LegendPosition = `${LEGEND_POSITION}`;
export type ColumnInfo = {name: string, type?: string, semType?: string};
export type CsvImportOptions = {
  delimiter?: string, decimalSeparator?: string, thousandSeparator?: string, headerRow?: boolean,
  columnFilterNames?: string[], columnFilterRegexp?: string, mergeDelimiters?: boolean, maxRows?: number,
  rowFilterTop?: number, rowFilterProb?: number, nullStrings?: string[], columnImportOptions?: ColumnInfo[]};
export type IndexPredicate = (ind: number) => boolean;
export type StringPredicate = (str: string) => boolean;
export type ScriptLanguage = `${SCRIPT_LANGUAGE}`;

export type ElementOptions = {
  id?: string;
  classes?: string;
  style?: object;
  processNode?: (node: HTMLElement) => void;
  onClick?: (node: HTMLElement) => void;
};

/** Metadata associated with the semantic type. */
export interface SemTypeInfo {

  /** Semantic type id */
  name: string;

  /** Semantic type description (shown in tooltips, etc) */
  description: string;

  /** Specifies the value data type.
   * For example, the `itemType` for semantic type `Molecule` is `String`.
   * Used for the automatic detection of semantic data types. */
  itemType?: ColumnType;

  /** Regular expression to check against the column name
   * Used for the automatic detection of semantic data types. */
  columnNameRegexp?: string;

  /** Regular expression to check against the values (only applies to strings)
   * Used for the automatic detection of semantic data types. */
  valueRegexp?: string;
}
