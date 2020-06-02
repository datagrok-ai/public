/**
 * @typedef {string} AggregationType
 * @typedef {string} SyncType
 * @typedef {string} JoinType
 * @typedef {string} ColumnType
 * @typedef {string} ViewerType
 * @typedef {string} ObjectType
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

/** @enum {JoinType} */
export const JOIN_TYPE = {
    INNER: 'inner',
    OUTER: 'outer',
    LEFT: 'left',
    RIGHT: 'right'
}

/** @enum {ColumnType} */
export const COLUMN_TYPE = {
    STRING : 'string',
    INT : 'int',
    FLOAT : 'double',
    BOOL : 'bool',
    DATE_TIME : 'datetime',
    BIG_INT : 'bigint'
}

/** @enum {TYPE} */
export const TYPE = {
    INT : 'int',
    BIG_INT : 'bigint',
    FLOAT : 'double',
    NUM : 'num',
    BOOL : 'bool',
    STRING : 'string',
    STRING_LIST : 'string_list',
    DATE_TIME : 'datetime',
    OBJECT : 'object',
    DATA_FRAME : 'dataframe',
    DATA_FRAME_LIST : 'dataframe_list',
    CELL : 'cell',
    COLUMN : 'column',
    COLUMN_LIST : 'column_list',
    GRAPHICS : 'graphics',
    ROW_FILTER : 'tablerowfiltercall',
    COLUMN_FILTER : 'colfiltercall',
    BIT_SET : 'bitset',
    MAP : 'map',
    DYNAMIC : 'dynamic',
    VIEWER : 'viewer',  // [ViewerBase] subclasses
    LIST : 'list',
    SEM_VALUE : 'semantic_value',
    FUNC_CALL : 'funccall',
    PROPERTY : 'property',
    CATEGORICAL : 'categorical',
    NUMERICAL : 'numerical',
    GRID_CELL_RENDER_ARGS : 'GridCellRenderArgs',

    VIEW : 'view',
    TABLE_VIEW : 'TableView',
    USER : 'User',
    MENU : 'Menu',
    PROJECT : 'Project',
    SEMANTIC_VALUE : 'semantic_value',
    EVENT_DATA : 'event_data'
}

export const TYPES_SCALAR = new Set([TYPE.INT, TYPE.BIG_INT, TYPE.FLOAT, TYPE.NUM, TYPE.BOOL, TYPE.STRING]);

export const VIEW_TYPE = {
    TABLE_VIEW : 'TableView'
};

/////// Semantic typesSEMTYPE_
/** @enum {SEMTYPE} */
export const SEMTYPE = {
    EMAIL : 'Email Address',
    URL : 'URL',
    PHONE_NUMBER : 'Phone Number',
    CITY : 'City',
    COUNTRY : 'Country',
    GENDER : 'Gender',
    STATE : 'State',
    COUNTY : 'County',
    PLACE_NAME : 'Place Name',
    ZIP_CODE : 'Zip Code',
    AREA_CODE : 'Area Code',
    STREET_ADDRESS : 'Street Address',
    TEXT : 'Text',
    DURATION : 'Duration',
    LATITUDE : 'Latitude',
    LONGITUDE : 'Longitude',
    IP_ADDRESS : 'IP Address',
    MOLECULE : 'Molecule',
    SUBSTRUCTURE : 'Substructure',
    MONEY : 'Money',
    IMAGE : 'Image',
    FILE : 'File',
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
}

////// Viewers
/** @enum {VIEWER} */
export const VIEWER = {
    HISTOGRAM : 'Histogram',
    BAR_CHART : 'Bar chart',
    BOX_PLOT : 'Box plot',
    CALENDAR : 'Calendar',
    CORR_PLOT : 'Correlation plot',
    DENSITY_PLOT : 'Density plot',
    FILTERS : 'Filters',
    FORM : 'Form',
    GLOBE : 'Globe',
    GRID : 'Grid',
    GOOGLE_MAP : 'Google map',
    HEAT_MAP : 'Heat map',
    LINE_CHART : 'Line chart',
    SHAPE_MAP : 'Shape map',
    MARKUP : 'Markup',
    MATRIX_PLOT : 'Matrix plot',
    NETWORK_DIAGRAM : 'Network diagram',
    PC_PLOT : 'PC Plot',
    PIE_CHART : 'Pie chart',
    SCATTER_PLOT : 'scatter plot',
    SCATTER_PLOT_3D : '3d scatter plot',
    STATISTICS : 'Statistics',
    TILE_VIEWER : 'Tile viewer',
    TREE_MAP : 'Tree map',
    TRELLIS_PLOT : 'Trellis plot',
    WORD_CLOUD : 'Word cloud'
}

/** @enum {SIMILARITY_METRIC} */
export const SIMILARITY_METRIC = {
    TANIMOTO : 'tanimoto',
    DICE : 'dice',
    COSINE : 'cosine',
    SOKAL : 'sokal',
    RUSSEL : 'russel',
    ROGOT_GOLDBERG : 'rogot-goldberg',
    KULCZYNSKI : 'kulczynski',
    MC_CONNAUGHEY : 'mc-connaughey',
    ASYMMETRIC : 'asymmetric',
    BRAUN_BLANQUET : 'braun-blanquet'
};
