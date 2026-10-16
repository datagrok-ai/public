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
  /** @deprecated Use {@link TYPE.SEMANTIC_VALUE}. Removed in 1.29. */
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

export const TYPES_SCALAR = new Set([TYPE.INT, TYPE.BIG_INT, TYPE.FLOAT, TYPE.NUM, TYPE.BOOL, TYPE.STRING]);

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

  FUNCTION_NAME: 'FunctionName',  // namespace-qualified function name, such as 'Chem:SmilesToMw'
  COLUMN_NAME: 'ColumnName',      // a string parameter that names a column; `options.table` names the parameter that supplies its table
  TABLE_NAME: 'TableName',        // a string parameter that names an open table
}
