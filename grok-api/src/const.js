/**
 * @typedef {string} AggregationType
 * @typedef {string} SyncType 
 **/

/** @enum {AggregationType} */
export const AGG = {
    KEY: "key",      // Special case: to be used in the 'group by' statement.
    PIVOT: "pivot",  // Special case: to be used as a pivot.
    FIRST: "first",
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
    Q3: "q3",
    SELECTED_ROWS_COUNT: "#selected"
};

/** @enum {SyncType} */
export const SYNC_TYPE = {
    CURRENT_ROW_TO_ROW: 'row to row',
    CURRENT_ROW_TO_SELECTION: 'row to selection',
    CURRENT_ROW_TO_FILTER: 'row to filter',
    MOUSE_OVER_ROW_TO_SELECTION: 'mouse-over to selection',
    MOUSE_OVER_ROW_TO_FILTER: 'mouse-over to filter',
    FILTER_TO_FILTER: 'filter to filter',
    FILTER_TO_SELECTION: 'filter to selection',
    SELECTION_TO_FILTER: 'selection to filter',
    SELECTION_TO_SELECTION: 'selection to selection',
};

export const INT_NULL = -2147483648;
export const FLOAT_NULL = 2.6789344063684636e-34;
export const DATE_TIME_NULL = -62135578800000000.0;

export const JOIN_TYPE_INNER = 'inner';
export const JOIN_TYPE_OUTER = 'outer';
export const JOIN_TYPE_LEFT = 'left';
export const JOIN_TYPE_RIGHT = 'right';

export const TYPE_INT = 'int';
export const TYPE_BIG_INT = 'bigint';
export const TYPE_FLOAT = 'double';
export const TYPE_NUM = 'num';
export const TYPE_BOOL = 'bool';
export const TYPE_STRING = 'string';
export const TYPE_STRING_LIST = 'string_list';
export const TYPE_DATE_TIME = 'datetime';
export const TYPE_OBJECT = 'object';
export const TYPE_PROJECT = 'project';
export const TYPE_DATA_FRAME = 'dataframe';
export const TYPE_DATA_FRAME_LIST = 'dataframe_list';
export const TYPE_CELL = 'cell';
export const TYPE_COLUMN = 'column';
export const TYPE_COLUMN_LIST = 'column_list';
export const TYPE_GRAPHICS = 'graphics';
export const TYPE_ROW_FILTER = 'tablerowfiltercall';
export const TYPE_COLUMN_FILTER = 'colfiltercall';
export const TYPE_BIT_SET = 'bitset';
export const TYPE_MAP = 'map';
export const TYPE_DYNAMIC = 'dynamic';
export const TYPE_VIEWER = 'viewer';  // [ViewerBase] subclasses
export const TYPE_LIST = 'list';
export const TYPE_SEM_VALUE = 'semantic_value';
export const TYPE_FUNC_CALL = 'funccall';
export const TYPE_PROPERTY = 'property';
export const TYPE_CATEGORICAL = 'categorical';
export const TYPE_NUMERICAL = 'numerical';
export const TYPE_GRID_CELL_RENDER_ARGS = 'GridCellRenderArgs';

export const TYPE_VIEW = 'view';
export const TYPE_TABLE_VIEW = 'TableView';
export const TYPE_USER = 'User';
export const TYPE_MENU = 'Menu';
export const TYPE_SEMANTIC_VALUE = 'semantic_value';
export const TYPE_EVENT_DATA = 'event_data';

export const TYPES_SCALAR = new Set([TYPE_INT, TYPE_BIG_INT, TYPE_FLOAT, TYPE_NUM, TYPE_BOOL, TYPE_STRING]);

export const VIEW_TYPE_TABLE_VIEW = 'TableView';

/////// Semantic typesSEMTYPE_

export const SEMTYPE_EMAIL = 'Email Address';
export const SEMTYPE_URL = 'URL';
export const SEMTYPE_PHONE_NUMBER = 'Phone Number';
export const SEMTYPE_CITY = 'City';
export const SEMTYPE_COUNTRY = 'Country';
export const SEMTYPE_GENDER = 'Gender';
export const SEMTYPE_STATE = 'State';
export const SEMTYPE_COUNTY = 'County';
export const SEMTYPE_PLACE_NAME = 'Place Name';
export const SEMTYPE_ZIP_CODE = 'Zip Code';
export const SEMTYPE_AREA_CODE = 'Area Code';
export const SEMTYPE_STREET_ADDRESS = 'Street Address';
export const SEMTYPE_TEXT = 'Text';
export const SEMTYPE_DURATION = 'Duration';
export const SEMTYPE_LATITUDE = 'Latitude';
export const SEMTYPE_LONGITUDE = 'Longitude';
export const SEMTYPE_IP_ADDRESS = 'IP Address';
export const SEMTYPE_MOLECULE = 'Molecule';
export const SEMTYPE_SUBSTRUCTURE = 'Substructure';
export const SEMTYPE_MONEY = 'Money';
export const SEMTYPE_IMAGE = 'Image';
export const SEMTYPE_FILE = 'File';

/////// Stats

export const STATS_TOTAL_COUNT = "count";
export const STATS_VALUE_COUNT = "values";
export const STATS_UNIQUE_COUNT = "unique";
export const STATS_MISSING_VALUE_COUNT = "nulls";
export const STATS_MIN = "min";
export const STATS_MAX = "max";
export const STATS_SUM = "sum";
export const STATS_MED = "med";
export const STATS_AVG = "avg";
export const STATS_STDEV = "stdev";
export const STATS_VARIANCE = "variance";
export const STATS_SKEW = "skew";
export const STATS_KURT = "kurt";
export const STATS_Q1 = "q1";
export const STATS_Q2 = "q2";
export const STATS_Q3 = "q3";

/////// Tags

export const TAGS_DESCRIPTION = 'description';

export const TAGS_TOOLTIP = '.tooltip';

/** JSON-encoded list of strings to be used in a cell editor. Applicable for string columns only. */
export const TAGS_CHOICES = '.choices';

/** When set to 'true', switches the cell editor to a combo box that only allows to choose values
 from a list of already existing values in the column.
 Applicable for string columns only.
 See also [TAGS_CHOICES]. */
export const TAGS_AUTO_CHOICES = '.auto-choices';

////// Viewers
export const VIEWER_HISTOGRAM = 'Histogram';
export const VIEWER_BAR_CHART = 'Bar chart';
export const VIEWER_BOX_PLOT = 'Box plot';
export const VIEWER_CALENDAR = 'Calendar';
export const VIEWER_CORR_PLOT = 'Correlation plot';
export const VIEWER_DENSITY_PLOT = 'Density plot';
export const VIEWER_FILTERS = 'Filters';
export const VIEWER_FORM = 'Form';
export const VIEWER_GLOBE = 'Globe';
export const VIEWER_GRID = 'Grid';
export const VIEWER_GOOGLE_MAP = 'Google map';
export const VIEWER_HEAT_MAP = 'Heat map';
export const VIEWER_LINE_CHART = 'Line chart';
export const VIEWER_SHAPE_MAP = 'Shape map';
export const VIEWER_MARKUP = 'Markup';
export const VIEWER_MATRIX_PLOT = 'Matrix plot';
export const VIEWER_NETWORK_DIAGRAM = 'Network diagram';
export const VIEWER_PC_PLOT = 'PC Plot';
export const VIEWER_PIE_CHART = 'Pie chart';
export const VIEWER_SCATTER_PLOT = 'scatter plot';
export const VIEWER_SCATTER_PLOT_3D = '3d scatter plot';
export const VIEWER_STATISTICS = 'Statistics';
export const VIEWER_TILE_VIEWER = 'Tile viewer';
export const VIEWER_TREE_MAP = 'Tree map';
export const VIEWER_TRELLIS_PLOT = 'Trellis plot';
export const VIEWER_WORD_CLOUD = 'Word cloud';

export const METRIC_TANIMOTO = 'tanimoto';
export const METRIC_DICE = 'dice';
export const METRIC_COSINE = 'cosine';
export const METRIC_SOKAL = 'sokal';
export const METRIC_RUSSEL = 'russel';
export const METRIC_ROGOT_GOLDBERG = 'rogot-goldberg';
export const METRIC_KULCZYNSKI = 'kulczynski';
export const METRIC_MC_CONNAUGHEY = 'mc-connaughey';
export const METRIC_ASYMMETRIC = 'asymmetric';
export const METRIC_BRAUN_BLANQUET = 'braun-blanquet';
