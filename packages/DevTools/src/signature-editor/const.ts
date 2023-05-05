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


export const functionPropsLabels: {[_: string]: string} = {
  [FUNC_PROPS_FIELD.NAME]: 'Name',
  [FUNC_PROPS_FIELD.CONNECTION]: 'Connection',
  [FUNC_PROPS_FIELD.DESCRIPTION]: 'Description',
  [FUNC_PROPS_FIELD.LOGIN]: 'Author',
  [FUNC_PROPS_FIELD.HELP_URL]: 'Help URL',
  [FUNC_PROPS_FIELD.LANGUAGE]: 'Language',
  [FUNC_PROPS_FIELD.REFERENCE]: 'Reference',
  [FUNC_PROPS_FIELD.SAMPLE]: 'Sample',
  [FUNC_PROPS_FIELD.ENVIRONMENT]: 'Environment',
  [FUNC_PROPS_FIELD.TAGS]: 'Tags'
};

export const languages = ['javascript', 'python', 'r', 'julia', 'octave', 'nodejs', 'grok', 'sql'];

export enum DIRECTION {
  INPUT = 'input',
  OUTPUT = 'output',
}

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

export const getChoicesByName: {[_: string]: {}} = {
  [FUNC_PROPS_FIELD.LANGUAGE]: {choices: languages},
  [FUNC_PROPS_FIELD.DESCRIPTION]: {editor: 'textarea'},
  [FUNC_PARAM_FIELDS.DIRECTION]: {choices: [...Object.values(DIRECTION)]},
  [FUNC_PARAM_FIELDS.TYPE]: {choices: funcParamTypes},
  [OPTIONAL_TAG_NAME.COLUMNS]: {choices: ['numerical', 'categorical']},
  [OPTIONAL_TAG_NAME.TYPE]: {choices: ['numerical', 'categorical', 'dateTime']}     
}

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

export const highlightModeByLang: {[_: string]: string}  = {
  [LANGUAGE.JS]: 'javascript',
  [LANGUAGE.NODEJS]: 'javascript',
  [LANGUAGE.PYTHON]: 'python',
  [LANGUAGE.GROK]: 'python',
  [LANGUAGE.OCTAVE]: 'octave',
  [LANGUAGE.JULIA]: 'julia',
  [LANGUAGE.R]: 'r',
  [LANGUAGE.SQL]: 'sql'
};

export const functionPropsCode: {[_: string]: string} = {
  [FUNC_PROPS_FIELD.NAME]: 'name',
  [FUNC_PROPS_FIELD.CONNECTION]: 'connection',
  [FUNC_PROPS_FIELD.DESCRIPTION]: 'description',
  [FUNC_PROPS_FIELD.LOGIN]: 'author',
  [FUNC_PROPS_FIELD.HELP_URL]: 'helpUrl',
  [FUNC_PROPS_FIELD.LANGUAGE]: 'language',
  [FUNC_PROPS_FIELD.REFERENCE]: 'reference',
  [FUNC_PROPS_FIELD.SAMPLE]: 'sample',
  [FUNC_PROPS_FIELD.ENVIRONMENT]: 'environment',
  [FUNC_PROPS_FIELD.TAGS]: 'tags'
}

export const enum FUNC_PARAM_FIELDS {
  DIRECTION = 'direction',
  TYPE = 'propertyType',
  NAME = 'name',
  DEFAULT_VALUE = 'defaultValue',
  DESCRIPTION = 'description',
  CATEGORY = 'category',
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

export const headerSign: {[_: string]: string} = {
  [LANGUAGE.JS]: '//',
  [LANGUAGE.NODEJS]: '//',
  [LANGUAGE.R]: '#',
  [LANGUAGE.GROK]: '#',
  [LANGUAGE.JULIA]: '#',
  [LANGUAGE.PYTHON]: '#',
  [LANGUAGE.OCTAVE]: '#',
  [LANGUAGE.SQL]: '--'
};

export const DATA_QUERY_VIEW = 'DataQueryView';