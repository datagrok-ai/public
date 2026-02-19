
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
  FILE = 'file',
  BLOB = 'blob',
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

  ELEMENT = 'element',
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

export enum GRID_COLUMN_TAGS {

}

/** Commonly used options on the function level */
export enum FUNC_OPTIONS {
  DEMO_PATH = 'demoPath',    // Demo path, such as 'Viewers | Radar'
  IS_DEMO_SCRIPT = 'isDemoScript'
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
  TABLE_VIEW = 'TableView',
  APPS = 'apps',
  SETTINGS = 'settings',
  WELCOME = 'welcome',
  SCRIPT = 'script',
  SKETCH = 'sketch',
  FORUM = 'forum',
  PROJECTS = 'projects',
  NOTEBOOKS = 'notebooks',
  HELP = 'help',
  OPEN_TEXT = 'text',
  DATABASES = 'databases',
  WEB_SERVICES = 'webservices',
  VIEW_LAYOUTS = 'layouts',
  FUNCTIONS = 'functions',
  DATA_CONNECTIONS = 'connections',
  DATA_JOB_RUNS = 'jobs',
  FILES = 'files',
  DATA_QUERY_RUNS = 'queryruns',
  EMAILS = 'emails',
  GROUPS = 'groups',
  MODELS = 'models',
  QUERIES = 'queries',
  SCRIPTS = 'scripts',
  USERS = 'users',
  PACKAGES = 'packages',
  PACKAGE_REPOSITORIES = 'repositories',
  JS_EDITOR = 'js',
  BROWSE = 'browse',
  HOME = 'datagrok',
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
  MOLECULE3D: 'Molecule3D',

  CONCENTRATION: 'Concentration',
  VOLUME: 'Volume',

  PDB_ID: 'PDB_ID',
  NEWICK: 'Newick',
  HELM: 'HELM',
  SUBSTRUCTURE: 'Substructure',
  MONEY: 'Money',
  IMAGE: 'Image',
  FILE: 'File',
  CHEMICAL_REACTION: 'ChemicalReaction',

  IC50: 'IC50',  //	[nM, µM] Half-maximal inhibitory concentration (lower = more potent)
  EC50: 'EC50',  // [nM, µM] Half-maximal effective concentration
  Ki: 'Ki',      // [nM, µM] Inhibition constant (binding affinity to target)
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
  COLOR_CODING_LINEAR_BELOW_MIN_COLOR: '.%color-coding-linear-below-min-color',
  COLOR_CODING_LINEAR_ABOVE_MAX_COLOR: '.%color-coding-linear-above-max-color',
  COLOR_CODING_LINEAR_ABSOLUTE: '.%color-coding-linear-absolute',
  COLOR_CODING_LINEAR_IS_ABSOLUTE: '.%color-coding-linear-is-absolute',
  COLOR_CODING_SCHEME_MAX: '.color-coding-scheme-max',
  COLOR_CODING_SCHEME_MIN: '.color-coding-scheme-min',
  COLOR_CODING_MATCH_TYPE: '.color-coding-match-type',
  COLOR_CODING_FALLBACK_COLOR: '.color-coding-fallback-color',

  MARKER_CODING: '.marker-coding',
  FORMULA_LINES: '.formula-lines',
  ANNOTATION_REGIONS: '.annotation-regions',

  /** When a dataframe is loaded from a CSV, the maximum number of significant digits
   in the fractional part for each numeric column is determined  */
  SOURCE_PRECISION: '.source-precision',
  /** Set on a dataframe column; used to format column contents in grids, CSV export, passing to scripts */
  FORMAT: 'format',
  FORMULA: 'formula',
  SEMTYPE: 'quality',

  /** Separator used to parse a cell value into multiple values for filter categories. */
  MULTI_VALUE_SEPARATOR: '.multi-value-separator',
  /** Boolean flag to control custom filters visibility. */
  IGNORE_CUSTOM_FILTER: '.ignore-custom-filter',
  /** Filter type for molecular columns: "Sketch" | "Categorical". See [DG.STRUCTURE_FILTER_TYPE] */
  STRUCTURE_FILTER_TYPE: '.structure-filter-type',
  /** Custom filter type to be used by default for a column: "<PackageName\>:<FilterType\>".
   * Takes precedence over [IGNORE_CUSTOM_FILTER] */
  CUSTOM_FILTER_TYPE: '.custom-filter-type',

  CELL_RENDERER: 'cell.renderer',
  UNITS: 'units',  // see DG.UNITS

  FRIENDLY_NAME: 'friendlyName',
  ALLOW_RENAME: '.allow-rename',

  CHEM: {
    SCAFFOLD: 'chem-scaffold'
  },

  LINK_CLICK_BEHAVIOR: '.linkClickBehavior',
  GROUP: 'group',
}

export enum LINK_CLICK_BEHAVIOR {
  OPEN_IN_NEW_TAB = 'Open in new tab',
  OPEN_IN_CONTEXT_PANEL = 'Open in context panel',
  CUSTOM = 'Custom',
}

export const FUNC_TYPES = {
  /** An application that gets shown in the app store.
    * Signature: app() */
  APP: 'app',

  /** Context-specific widget that appears on the context panel
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
  CONVERTER: 'converter',

  WIDGET: 'widget',
  WIDGETS: 'widgets',
  EDITOR: 'editor',
  TRANSFORM: 'Transform',
  FILTER: 'filter',
  VIEWER: 'viewer',
  VALUE_EDITOR: 'valueEditor',
  CELL_EDITOR: 'cellEditor',
  UNIT_CONVERTER: 'unitConverter',
  MOLECULE_SKETCHER: 'moleculeSketcher',
  TOOLTIP: 'tooltip',
  FOLDER_VIEWER: 'folderViewer',
  SCRIPT_HANDLER: 'scriptHandler',

  HIT_TRIAGE_FUNCTION: 'HitTriageFunction',
  HIT_TRIAGE_DATA_SOURCE: 'HitTriageDataSource',
  HIT_TRIAGE_SUBMIT_FUNCTION: 'HitTriageSubmitFunction',
  HIT_DESIGNER_FUNCTION: 'HitDesignerFunction',

  DIM_RED_PREPROCESS: 'dim-red-preprocessing-function',
  DIM_RED_POSTPROCESS: 'dim-red-postprocessing-function',

  MONOMER_LIB_PROVIDER: 'monomer-lib-provider',

  SEARCH_PROVIDER: 'searchProvider',
  NOTATION_REFINER: 'notationRefiner',
}


export interface FuncRoleDescription {
  role: string;
  description: string;
  header: string;
  signature?: string;
}


export const functionRoles: FuncRoleDescription[] = [
  {
    role: FUNC_TYPES.APP,
    description: 'An application that gets shown in the app store.',
    header: 'tags',
    signature: 'app(): void | View'
  },
  {
    role: FUNC_TYPES.PANEL,
    description: 'Context-specific widget that appears on the context panel.',
    header: 'tags',
    signature: 'panel(...args): Widget | Viewer | graphics | void'
  },
  {
    role: FUNC_TYPES.INIT,
    description: 'Gets invoked when the containing package is initialized.',
    header: 'tags',
    signature: 'init(): void'
  },
  {
    role: FUNC_TYPES.AUTOSTART,
    description: 'Gets invoked at platform startup. Use it wisely as the whole package will get initialized.',
    header: 'tags',
    signature: 'autostart(): void'
  },
  {
    role: FUNC_TYPES.SEM_TYPE_DETECTOR,
    description: 'Semantic type detector for a column. Gets invoked when a new dataframe is imported into the platform.\n   *  Implementation should either set column.semType directly, or return the semantic type that will get assigned.',
    header: 'tags',
    signature: 'semTypeDetector(col: Column): string'
  },
  {
    role: FUNC_TYPES.FILE_VIEWER,
    header: 'tags',
    description: 'Creates a viewer (or editor) for a file with the specified extension.\n   *  The extension is derived from the `fileViewer-[extension]` tag.\n   *  Used in the file system browser.',
    signature: 'fileViewer(file: FileInfo): View'
  },
  {
    role: FUNC_TYPES.FILE_EXPORTER,
    header: 'tags',
    description: 'Exports a file. Gets added to the "export" menu at startup.',
    signature: 'fileExporter(): void'
  },
  {
    role: FUNC_TYPES.FILE_IMPORTER,
    header: 'tags',
    description: 'Handles custom file formats.\n   * The `meta.ext` parameter should contain a comma-separated list of extensions',
    signature: 'fileImporter(x: string | TypedArray): DataFrame[]'
  },
  {
    role: FUNC_TYPES.CELL_RENDERER,
    header: 'tags',
    description: 'Creates a cell renderer that is used for rendering cells for specific semantic types.\n   *  Semantic type is derived from the `cellRenderer-[semType]` tag.',
    signature: 'cellRenderer(): GridCellRenderer'
  },
  {
    role: FUNC_TYPES.PACKAGE_SETTINGS_EDITOR,
    header: 'tags',
    description: 'Edits package settings.',
    signature: 'packageSettingsEditor(): Widget'
  },
  {
    role: FUNC_TYPES.DASHBOARD,
    description: 'Makes a widget appear on the welcome screen.',
    header: 'tags',
    signature: 'dashboard(): Widget'
  },
  {
    role: FUNC_TYPES.FUNCTION_ANALYSIS,
    description: 'Function analysis that gets added to the function view. Examples: sensitivity analysis, parameter editor',
    header: 'tags',
    signature: 'functionAnalysis(x: Function): View'
  },
  {
    role: FUNC_TYPES.CONVERTER,
    description: 'Converts values. Has one input and one output',
    header: 'role',
    signature: 'converter(x: any): any'
  },
  {
    role: FUNC_TYPES.WIDGET,
    description: 'Creates a custom widget.',
    header: 'tags',
    signature: 'widget(...args): Widget',
  },
  {
    role: FUNC_TYPES.WIDGETS,
    description: 'Creates a custom widget.',
    header: 'tags',
    signature: 'widgets(...args): Widget',
  },
  {
    role: FUNC_TYPES.EDITOR,
    description: 'Creates a custom editor for a function call',
    header: 'tags',
    signature: 'editor(call: FuncCall): Widget | View | void',
  },
  {
    role: FUNC_TYPES.TRANSFORM,
    description: '',
    header: 'tags',
    signature: 'transform(table: DataFrame, ...args): any',
  },
  {
    role: FUNC_TYPES.FILTER,
    description: 'Creates a custom table filter',
    header: 'tags',
    signature: 'filter(): Filter',
  },
  {
    role: FUNC_TYPES.VIEWER,
    description: 'Creates a custom viewer',
    header: 'tags',
    signature: 'viewer(): Viewer',
  },
  {
    role: FUNC_TYPES.VALUE_EDITOR,
    description: 'Custom editor for specific input types',
    header: 'tags',
    signature: 'valueEditor(...args): any',
  },
  {
    role: FUNC_TYPES.CELL_EDITOR,
    description: 'CCustom editor for grid cells',
    header: 'tags',
    signature: 'cellEditor(cell: GridCell): void',
  },
  {
    role: FUNC_TYPES.UNIT_CONVERTER,
    description: 'Converts between measurement units',
    header: 'tags',
    signature: 'unitConverter(value: string, source: string, target: string): string',
  },
  {
    role: FUNC_TYPES.MOLECULE_SKETCHER,
    description: 'Creates a molecule sketcher widget',
    header: 'tags',
    signature: 'moleculeSketcher(): Widget',
  },
  {
    role: FUNC_TYPES.TOOLTIP,
    description: 'Provides a custom tooltip for columns',
    header: 'tags',
    signature: 'tooltip(col: Column): Widget',
  },
  {
    role: FUNC_TYPES.FOLDER_VIEWER,
    description: 'Provides a custom folder content preview.',
    header: 'tags',
    signature: 'folderViewer(folder: File, files: list<file>): Widget',
  },
  {
    role: FUNC_TYPES.HIT_TRIAGE_FUNCTION,
    description: 'Compute function for Hit Triage campaigns.',
    header: 'tags',
    signature: 'hitTriageFunction(table: DataFrame, moleculeCol: Column): DataFrame'
  },
  {
    role: FUNC_TYPES.HIT_TRIAGE_DATA_SOURCE,
    description: 'Provides a datasource for Hit Triage campaigns. Must return a dataframe containing molecules.',
    header: 'tags',
    signature: 'hitTriageDataSource(...args): DataFrame'
  },
  {
    role: FUNC_TYPES.HIT_TRIAGE_SUBMIT_FUNCTION,
    description: 'Processes or saves the computed dataset.',
    header: 'tags',
    signature: 'hitTriageSubmitFunction(df: DG.DataFrame, moleculesCol: string): void'
  },
  {
    role: FUNC_TYPES.HIT_DESIGNER_FUNCTION,
    description: 'Compute function for Hit Design campaigns.',
    header: 'tags',
    signature: 'hitDesignerFunction(molecule: string, ...args): DataFrame'
  },
  {
    role: FUNC_TYPES.DIM_RED_PREPROCESS,
    description: 'Preprocessing function for dimensionality reduction.',
    header: 'tags',
    signature: 'preprocess(col: Column, metric: string, ...args): any'
  },
  {
    role: FUNC_TYPES.DIM_RED_POSTPROCESS,
    description: 'Postprocessing function for dimensionality reduction.',
    header: 'tags',
    signature: 'postprocess(xCol: Column, yCol: Column, ...args): void'
  },
  {
    role: FUNC_TYPES.MONOMER_LIB_PROVIDER,
    description: 'Provides a monomer library provider.',
    header: 'tags',
    signature: 'getMonomerLibProvider(): IMonomerLibProvider'
  },
  {
    role: FUNC_TYPES.SEARCH_PROVIDER,
    description: 'Marks a function to be used as a search provider in the global search.',
    header: 'tags',
    signature: 'searchProvider(): SearchProvider'
  },
  {
    role: FUNC_TYPES.NOTATION_REFINER,
    description: 'Refines the biological sequence notation based on company specific rules',
    header: 'tags',
    signature: 'notationRefiner(column: Column, stats: any, separator: string): bool'
  }
]

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
  TIMELINES = 'Timelines',
  RADAR_VIEWER = 'Radar',
  SURFACE_PLOT = 'Surface plot',
  SCAFFOLD_TREE = 'Scaffold Tree',
  PIVOT_TABLE = 'Pivot table',
  CONFUSION_MATRIX = 'Confusion matrix'
}

/** @enum {CORE_VIEWER} */
export enum CORE_VIEWER {
  HISTOGRAM = 'Histogram',
  BAR_CHART = 'Bar chart',
  BOX_PLOT = 'Box plot',
  CALENDAR = 'Calendar',
  CORR_PLOT = 'Correlation plot',
  DENSITY_PLOT = 'Density plot',
  FILTERS = 'Filters',
  FORM = 'Form',
  GRID = 'Grid',
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
  PIVOT_TABLE = 'Pivot table',
  CONFUSION_MATRIX = 'Confusion matrix'
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

/** @enum {STRUCTURE_FILTER_TYPE} */
export enum STRUCTURE_FILTER_TYPE {
  Sketch = 'Sketch',
  Categorical = 'Categorical'
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
  LEFT = 'left',
  RIGHT = 'right',
  TOP = 'up',
  DOWN = 'down',
  FILL = 'fill',
}

export enum LEGEND_POSITION {
  AUTO = 'auto',
  LEFT = 'left',
  RIGHT = 'right',
  TOP = 'top',
  BOTTOM = 'bottom',
  // TODO: commented temporarily
  // LEFT_TOP = 'left top',
  // LEFT_BOTTOM = 'left bottom',
  // RIGHT_TOP = 'right top',
  // RIGHT_BOTTOM = 'right bottom'
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
  PYODIDE = 'pyodide'
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

export enum PERMISSION {
  EDIT = 'Edit',
  VIEW = 'View',
  SHARE = 'Share',
  DELETE = 'Delete',
}

/**
 * Event type identifiers used with `grok.events.onEvent()` and `__obs()`.
 *
 * Events are grouped by prefix:
 * - `d4-*`: Platform/UI events
 * - `grok-*`: Grok-specific events (views, projects)
 * - `ddt-*`: DataFrame/table events (used internally via DataFrame event getters)
 *
 * @see Events class in events.ts for typed event getters
 * @see {@link ./dart-interop.md} for event system documentation
 */
export enum EVENT_TYPE {
  // Context menu events
  CONTEXT_MENU = 'd4-context-menu',
  CONTEXT_MENU_CLOSED = 'd4-menu-closed',

  // View events
  CURRENT_VIEW_CHANGED = 'd4-current-view-changed',
  CURRENT_VIEW_CHANGING = 'd4-current-view-changing',
  VIEW_ADDED = 'grok-view-added',
  VIEW_ADDING = 'grok-view-adding',
  VIEW_REMOVED = 'grok-view-removed',
  VIEW_REMOVING = 'grok-view-removing',
  VIEW_RENAMED = 'grok-view-renamed',
  VIEW_CHANGED = 'grok-view-changed',
  VIEW_CHANGING = 'grok-view-changing',

  // Object events
  CURRENT_OBJECT_CHANGED = 'd4-current-object-changed',
  CURRENT_CELL_CHANGED = 'd4-current-cell-changed',

  // Table events
  TABLE_ADDED = 'd4-table-added',
  TABLE_REMOVED = 'd4-table-removed',

  // Query events
  QUERY_STARTED = 'd4-query-started',
  QUERY_FINISHED = 'd4-query-finished',

  // Project events
  CURRENT_PROJECT_CHANGED = 'grok-current-project-changed',
  PROJECT_SAVED = 'grok-project-saved',
  PROJECT_SAVING = 'grok-project-saving',
  PROJECT_OPENED = 'grok-project-opened',
  PROJECT_CLOSING = 'grok-project-closing',
  PROJECT_CLOSED = 'grok-project-closed',
  PROJECT_MODIFIED = 'grok-project-modified',

  // Viewer events
  VIEWER_ADDED = 'd4-viewer-added',
  VIEWER_CLOSED = 'd4-viewer-closed',

  // Layout events
  VIEW_LAYOUT_GENERATED = 'd4-view-layout-generated',
  VIEW_LAYOUT_APPLYING = 'd4-view-layout-applying',
  VIEW_LAYOUT_APPLIED = 'd4-view-layout-applied',

  // Tooltip events
  TOOLTIP_REQUEST = 'd4-tooltip-request',
  TOOLTIP_SHOWN = 'd4-tooltip-shown',
  TOOLTIP_CLOSED = 'd4-tooltip-closed',

  // Dialog and input events
  DIALOG_SHOWN = 'd4-dialog-showed',
  INPUT_CREATED = 'd4-input-created',
  ACCORDION_CONSTRUCTED = 'd4-accordion-constructed',

  // Form events
  FORM_CREATING = 'd4-form-creating',

  // Filter events
  RESET_FILTER_REQUEST = 'd4-reset-filter-request',

  // File events
  FILE_EDITED = 'grok-file-edited',
  FILE_IMPORT_REQUEST = 'd4-file-import-request',

  // Grid events
  GRID_CELL_LINK_CLICKED = 'd4-grid-cell-link-clicked-global',

  // Package events
  PACKAGE_LOADED = 'd4-package-loaded',

  // AI events
  AI_GENERATION_ABORT = 'd4-ai-generation-abort',
  AI_PANEL_TOGGLE = 'd4-ai-panel-toggle',

  // Tree view events
  TREE_VIEW_NODE_ADDED = 'd4-tree-view-node-added',

  // Server events
  SERVER_MESSAGE = 'server-message',
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
export type ViewType = `${VIEW_TYPE}` | string;
export type ObjectType = string;
export type ViewerPropertyType = string;
export type Type = `${TYPE}`;
export type SemType = string;
export type SimilarityMetric = `${SIMILARITY_METRIC}`;
export type StructureFilterType = `${STRUCTURE_FILTER_TYPE}`;
export type ColorCodingType = `${COLOR_CODING_TYPE}`;
export type MarkerCodingType = `${MARKER_TYPE}`;
export type DemoDatasetName = `${DEMO_DATASET}`;
export type DockType = `${DOCK_TYPE}`;
export type LegendPosition = `${LEGEND_POSITION}`;
export type CsvImportColumnOptions = {name: string, type?: string, semType?: string};
export type CsvImportOptions = {
  delimiter?: string, decimalSeparator?: string, thousandSeparator?: string, headerRow?: boolean,
  columnFilterNames?: string[], columnFilterRegexp?: string, mergeDelimiters?: boolean, maxRows?: number,
  doublePrecision?: boolean, rowFilterTop?: number, rowFilterProb?: number, nullStrings?: string[], columnImportOptions?: CsvImportColumnOptions[]};
export type IndexPredicate = (ind: number) => boolean;
export type StringPredicate = (str: string) => boolean;
export type ScriptingLanguage = `${SCRIPT_LANGUAGE}`;
type CSSProperties = Partial<Record<keyof CSSStyleDeclaration, string>>;

export type ElementOptions = {
  id?: string;
  classes?: string;
  style?: CSSProperties;
  //tooltip?: string;
  processNode?: (node: HTMLElement) => void;
  onClick?: (event: PointerEvent) => void;
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
