import * as DG from 'datagrok-api/dg';

//help with examples

export const DEFAULT_CATEGORY = 'Misc';

export enum FUNC_PROPS_FIELD {
  NAME = 'name',
  DESCRIPTION = 'description',
  LANGUAGE = 'language',
  HELP_URL = 'helpUrl',
  REFERENCE = 'reference',
  LOGIN = 'author?.login',
  SAMPLE = 'sample',
  ENVIRONMENT = 'environment',
  TAGS = 'tags',
  CONNECTION = "connection"
}

export const tooltipMessage = {
  'caption': 'Custom field caption',
  'postfix': 'Field postfix',
  'units': 'Value unit name',
  'editor': 'Editor',
  'semType': 'Semantic type',
  'columns': 'Numerical or categorical columns will be loaded',
  'type': 'Column type',
  'format': 'Datetime format',
  'allowNulls': 'Validation of the missing values presence',
  'action': 'For output parameters only',
  'choices': 'List of choices for string parameter',
  'suggestions': 'List of suggestions for string parameter',
  'min': 'Minimum value',
  'max': 'Maximum value',
}

export const obligatoryFuncProps = ['name', 'description', 'helpUrl', 'language'];

export const functionPropsLabels = (key: FUNC_PROPS_FIELD) => {
  switch (key) {
    case FUNC_PROPS_FIELD.NAME: return 'Name';
    case FUNC_PROPS_FIELD.CONNECTION: return 'Connection';
    case FUNC_PROPS_FIELD.DESCRIPTION: return 'Description';
    case FUNC_PROPS_FIELD.LOGIN: return 'Author';
    case FUNC_PROPS_FIELD.HELP_URL: return 'Help URL';
    case FUNC_PROPS_FIELD.LANGUAGE: return 'Language';
    case FUNC_PROPS_FIELD.REFERENCE: return 'Reference';
    case FUNC_PROPS_FIELD.SAMPLE: return 'Sample';
    case FUNC_PROPS_FIELD.ENVIRONMENT: return 'Environment';
    case FUNC_PROPS_FIELD.TAGS: return 'Tags';
  }
};


export const getChoicesByName = (name) => {
    switch (name) {
        case FUNC_PROPS_FIELD.LANGUAGE: return {choices: languages};
        case FUNC_PROPS_FIELD.DESCRIPTION: return {editor: 'textarea'};
        case FUNC_PARAM_FIELDS.DIRECTION: return {choices: [...Object.values(DIRECTION)]};
        case FUNC_PARAM_FIELDS.TYPE: return {choices: funcParamTypes};
        case OPTIONAL_TAG_NAME.COLUMNS: return {choices: ['numerical', 'categorical']};
        case OPTIONAL_TAG_NAME.TYPE: return {choices: ['numerical', 'categorical', 'dateTime']}; 
    }
}
export const highlightModeByLang = (key: LANGUAGE) => {
  switch (key) {
    case LANGUAGE.JS:
    case LANGUAGE.NODEJS: return 'javascript';
    case LANGUAGE.PYTHON:
    case LANGUAGE.GROK: return 'python';
    case LANGUAGE.OCTAVE: return 'octave';
    case LANGUAGE.JULIA: return 'julia';
    case LANGUAGE.R: return 'r';
    case LANGUAGE.SQL: return 'sql';
  }
};

export const functionPropsCode = (key: string) => {
  switch (key) {
    case FUNC_PROPS_FIELD.NAME: return 'name';
    case FUNC_PROPS_FIELD.CONNECTION: return 'connection';
    case FUNC_PROPS_FIELD.DESCRIPTION: return 'description';
    case FUNC_PROPS_FIELD.LOGIN: return 'author';
    case FUNC_PROPS_FIELD.HELP_URL: return 'helpUrl';
    case FUNC_PROPS_FIELD.LANGUAGE: return 'language';
    case FUNC_PROPS_FIELD.REFERENCE: return 'reference';
    case FUNC_PROPS_FIELD.SAMPLE: return 'sample';
    case FUNC_PROPS_FIELD.ENVIRONMENT: return 'environment';
    case FUNC_PROPS_FIELD.TAGS: return 'tags';
  }
};

export const funcParamTypes = [
  DG.TYPE.BOOL,
  DG.TYPE.COLUMN_LIST,
  DG.TYPE.COLUMN,
  DG.TYPE.DATA_FRAME,
  DG.TYPE.DATE_TIME,
  DG.TYPE.FLOAT,
  DG.TYPE.GRAPHICS,
  DG.TYPE.INT,
  DG.TYPE.STRING,
];

export const enum FUNC_PARAM_FIELDS {
  DIRECTION = 'direction',
  TYPE = 'propertyType',
  NAME = 'name',
  DEFAULT_VALUE = 'defaultValue',
  DESCRIPTION = 'description',
  CATEGORY = 'category',
}

export enum DIRECTION {
  INPUT = 'input',
  OUTPUT = 'output',
}



export const functionParamsMapping = {
  [FUNC_PARAM_FIELDS.DIRECTION]: 'Direction',
  [FUNC_PARAM_FIELDS.TYPE]: 'Type',
  [FUNC_PARAM_FIELDS.NAME]: 'Name',
  [FUNC_PARAM_FIELDS.DEFAULT_VALUE]: 'Default value',
  [FUNC_PARAM_FIELDS.DESCRIPTION]: 'Description',
  [FUNC_PARAM_FIELDS.CATEGORY]: 'Category',
  'Direction': FUNC_PARAM_FIELDS.DIRECTION,
  'Type': FUNC_PARAM_FIELDS.TYPE,
  'Name': FUNC_PARAM_FIELDS.NAME,
  'Default value': FUNC_PARAM_FIELDS.DEFAULT_VALUE,
  'Description': FUNC_PARAM_FIELDS.DESCRIPTION,
  'Category': FUNC_PARAM_FIELDS.CATEGORY,
};

export enum COMMON_TAG_NAME {
  VALIDATORS = 'validators',
  CAPTION = 'caption',
  POSTFIX = 'postfix',
  UNITS = 'units',
  EDITOR = 'editor',
  SEM_TYPE = 'semType',
}

export enum OPTIONAL_TAG_NAME {
  COLUMNS = 'columns',
  TYPE = 'type',
  FORMAT = 'format',
  ALLOW_NULLS = 'allowNulls',
  ACTION = 'action',
  CHOICES = 'choices',
  SUGGESTIONS = 'suggestions',
  MIN = 'min',
  MAX = 'max',
}

export const COMMON_TAG_NAMES = [...Object.values(COMMON_TAG_NAME)];
export const OPTIONAL_TAG_NAMES = [...Object.values(OPTIONAL_TAG_NAME)];
export const FUNC_PROPS_FIELDS = [...Object.values(FUNC_PROPS_FIELD)];
export const DF_TAG_NAMES = [
  ...Object.values(COMMON_TAG_NAME),
  OPTIONAL_TAG_NAME.COLUMNS,
  OPTIONAL_TAG_NAME.ACTION];
export const COLUMN_TAG_NAMES = [...Object.values(COMMON_TAG_NAME), OPTIONAL_TAG_NAME.TYPE, OPTIONAL_TAG_NAME.FORMAT,
  OPTIONAL_TAG_NAME.ALLOW_NULLS, OPTIONAL_TAG_NAME.ACTION];
export const STRING_TAG_NAMES = [...Object.values(COMMON_TAG_NAME),
  OPTIONAL_TAG_NAME.CHOICES , OPTIONAL_TAG_NAME.SUGGESTIONS];
export const INT_TAG_NAMES = [...Object.values(COMMON_TAG_NAME), OPTIONAL_TAG_NAME.MIN, OPTIONAL_TAG_NAME.MAX];

export const optionTags = ((param: DG.Property) => {
  switch (param.propertyType) {
    case DG.TYPE.INT:
      return INT_TAG_NAMES;
    case DG.TYPE.DATA_FRAME:
      return DF_TAG_NAMES;
    case DG.TYPE.COLUMN_LIST:
    case DG.TYPE.COLUMN:
      return COLUMN_TAG_NAMES;
    case DG.TYPE.STRING:
      return STRING_TAG_NAMES;
    default:
      return COMMON_TAG_NAMES;
  }
});

export enum LANGUAGE {
  JS = 'javascript',
  PYTHON = 'python',
  R = 'r',
  JULIA = 'julia',
  OCTAVE = 'octave',
  NODEJS = 'nodejs',
  GROK = 'grok',
  SQL = 'sql'
};

export const languages = ['javascript', 'python', 'r', 'julia', 'octave', 'nodejs', 'grok', 'sql'];

export const headerSign = (lang: LANGUAGE) => {
  switch (lang) {
    case LANGUAGE.JS:
    case LANGUAGE.NODEJS:
      return '//';
    case LANGUAGE.R:
    case LANGUAGE.GROK:
    case LANGUAGE.JULIA:
    case LANGUAGE.PYTHON:
    case LANGUAGE.OCTAVE:
      return '#';
    case LANGUAGE.SQL:
      return '--';
  }
};

export const DATA_QUERY_VIEW = 'DataQueryView';