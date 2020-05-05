
const AGG_KEY = "key";      // Special case: to be used in the 'group by' statement.
const AGG_PIVOT = "pivot";  // Special case: to be used as a pivot.
const AGG_FIRST = "first";
const AGG_TOTAL_COUNT = "count";
const AGG_VALUE_COUNT = "values";
const AGG_UNIQUE_COUNT = "unique";
const AGG_MISSING_VALUE_COUNT = "nulls";
const AGG_MIN = "min";
const AGG_MAX = "max";
const AGG_SUM = "sum";
const AGG_MED = "med";
const AGG_AVG = "avg";
const AGG_STDEV = "stdev";
const AGG_VARIANCE = "variance";
const AGG_SKEW = "skew";
const AGG_KURT = "kurt";
const AGG_Q1 = "q1";
const AGG_Q2 = "q2";
const AGG_Q3 = "q3";
const AGG_SELECTED_ROWS_COUNT = "#selected";

const CURRENT_ROW_TO_ROW = 'row to row';
const CURRENT_ROW_TO_SELECTION = 'row to selection';
const CURRENT_ROW_TO_FILTER = 'row to filter';
const MOUSE_OVER_ROW_TO_SELECTION = 'mouse-over to selection';
const MOUSE_OVER_ROW_TO_FILTER = 'mouse-over to filter';
const FILTER_TO_FILTER = 'filter to filter';
const FILTER_TO_SELECTION = 'filter to selection';
const SELECTION_TO_FILTER = 'selection to filter';
const SELECTION_TO_SELECTION = 'selection to selection';

const JOIN_TYPE_INNER = 'inner';
const JOIN_TYPE_OUTER = 'outer';
const JOIN_TYPE_LEFT = 'left';
const JOIN_TYPE_RIGHT = 'right';

const TYPE_INT = 'int';
const TYPE_BIG_INT = 'bigint';
const TYPE_FLOAT = 'double';
const TYPE_NUM = 'num';
const TYPE_BOOL = 'bool';
const TYPE_STRING = 'string';
const TYPE_STRING_LIST = 'string_list';
const TYPE_DATE_TIME = 'datetime';
const TYPE_OBJECT = 'object';
const TYPE_PROJECT = 'project';
const TYPE_DATA_FRAME = 'dataframe';
const TYPE_DATA_FRAME_LIST = 'dataframe_list';
const TYPE_CELL = 'cell';
const TYPE_COLUMN = 'column';
const TYPE_COLUMN_LIST = 'column_list';
const TYPE_GRAPHICS = 'graphics';
const TYPE_ROW_FILTER = 'tablerowfiltercall';
const TYPE_COLUMN_FILTER = 'colfiltercall';
const TYPE_BIT_SET = 'bitset';
const TYPE_MAP = 'map';
const TYPE_DYNAMIC = 'dynamic';
const TYPE_VIEWER = 'viewer';  // [ViewerBase] subclasses
const TYPE_LIST = 'list';
const TYPE_SEM_VALUE = 'semantic_value';
const TYPE_FUNC_CALL = 'funccall';
const TYPE_PROPERTY = 'property';
const TYPE_CATEGORICAL = 'categorical';
const TYPE_NUMERICAL = 'numerical';

const TYPE_VIEW = 'view';
const TYPE_TABLE_VIEW = 'TableView';
const TYPE_USER = 'User';
const TYPE_MENU = 'Menu';
const TYPE_SEMANTIC_VALUE = 'semantic_value';
const TYPE_EVENT_DATA = 'event_data';

const TYPES_SCALAR = new Set([TYPE_INT, TYPE_BIG_INT, TYPE_FLOAT, TYPE_NUM, TYPE_BOOL, TYPE_STRING]);

const VIEW_TYPE_TABLE_VIEW = 'TableView';

/////// Semantic typesSEMTYPE_

const SEMTYPE_EMAIL = 'Email Address';
const SEMTYPE_URL = 'URL';
const SEMTYPE_PHONE_NUMBER = 'Phone Number';
const SEMTYPE_CITY = 'City';
const SEMTYPE_COUNTRY = 'Country';
const SEMTYPE_GENDER = 'Gender';
const SEMTYPE_STATE = 'State';
const SEMTYPE_COUNTY = 'County';
const SEMTYPE_PLACE_NAME = 'Place Name';
const SEMTYPE_ZIP_CODE = 'Zip Code';
const SEMTYPE_AREA_CODE = 'Area Code';
const SEMTYPE_STREET_ADDRESS = 'Street Address';
const SEMTYPE_TEXT = 'Text';
const SEMTYPE_DURATION = 'Duration';
const SEMTYPE_LATITUDE = 'Latitude';
const SEMTYPE_LONGITUDE = 'Longitude';
const SEMTYPE_IP_ADDRESS = 'IP Address';
const SEMTYPE_MOLECULE = 'Molecule';
const SEMTYPE_SUBSTRUCTURE = 'Substructure';
const SEMTYPE_MONEY = 'Money';
const SEMTYPE_IMAGE = 'Image';
const SEMTYPE_FILE = 'File';

/////// Stats

const STATS_TOTAL_COUNT = "count";
const STATS_VALUE_COUNT = "values";
const STATS_UNIQUE_COUNT = "unique";
const STATS_MISSING_VALUE_COUNT = "nulls";
const STATS_MIN = "min";
const STATS_MAX = "max";
const STATS_SUM = "sum";
const STATS_MED = "med";
const STATS_AVG = "avg";
const STATS_STDEV = "stdev";
const STATS_VARIANCE = "variance";
const STATS_SKEW = "skew";
const STATS_KURT = "kurt";
const STATS_Q1 = "q1";
const STATS_Q2 = "q2";
const STATS_Q3 = "q3";

/////// Tags

const TAGS_DESCRIPTION = 'description';

const TAGS_TOOLTIP = '.tooltip';

/** JSON-encoded list of strings to be used in a cell editor. Applicable for string columns only. */
const TAGS_CHOICES = '.choices';

/** When set to 'true', switches the cell editor to a combo box that only allows to choose values
 from a list of already existing values in the column.
 Applicable for string columns only.
 See also [TAGS_CHOICES]. */
const TAGS_AUTO_CHOICES = '.auto-choices';

////// Viewers
const VIEWER_HISTOGRAM = 'Histogram';
const VIEWER_BAR_CHART = 'Bar chart';
const VIEWER_BOX_PLOT = 'Box plot';
const VIEWER_CALENDAR = 'Calendar';
const VIEWER_CORR_PLOT = 'Correlation plot';
const VIEWER_DENSITY_PLOT = 'Density plot';
const VIEWER_FILTERS = 'Filters';
const VIEWER_FORM = 'Form';
const VIEWER_GLOBE = 'Globe';
const VIEWER_GRID = 'Grid';
const VIEWER_GOOGLE_MAP = 'Google map';
const VIEWER_HEAT_MAP = 'Heat map';
const VIEWER_LINE_CHART = 'Line chart';
const VIEWER_SHAPE_MAP = 'Shape map';
const VIEWER_MARKUP = 'Markup';
const VIEWER_MATRIX_PLOT = 'Matrix plot';
const VIEWER_NETWORK_DIAGRAM = 'Network diagram';
const VIEWER_PC_PLOT = 'PC Plot';
const VIEWER_PIE_CHART = 'Pie chart';
const VIEWER_SCATTER_PLOT = 'scatter plot';
const VIEWER_SCATTER_PLOT_3D = '3d scatter plot';
const VIEWER_STATISTICS = 'Statistics';
const VIEWER_TILE_VIEWER = 'Tile viewer';
const VIEWER_TREE_MAP = 'Tree map';
const VIEWER_TRELLIS_PLOT = 'Trellis plot';
const VIEWER_WORD_CLOUD = 'Word cloud';

const METRIC_TANIMOTO = 'tanimoto';
const METRIC_DICE = 'dice';
const METRIC_COSINE = 'cosine';
const METRIC_SOKAL = 'sokal';
const METRIC_RUSSEL = 'russel';
const METRIC_ROGOT_GOLDBERG = 'rogot-goldberg';
const METRIC_KULCZYNSKI = 'kulczynski';
const METRIC_MC_CONNAUGHEY = 'mc-connaughey';
const METRIC_ASYMMETRIC = 'asymmetric';
const METRIC_BRAUN_BLANQUET = 'braun-blanquet';
