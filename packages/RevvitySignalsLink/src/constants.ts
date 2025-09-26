export type RevvitySignalsCompoundTypeMapping = {
    compoundType: string,
    viewName: string,
}

export const compoundTypeAndViewNameMapping = [
    {compoundType: 'asset', viewName: 'assets'},
    {compoundType: 'batch', viewName: 'batches'}
]

export const MOL_COL_NAME = 'molecule';
export const ID_COL_NAME = 'id';
export const HIDDEN_ID_COL_NAME = '~id';
export const NAME = 'name';

export const STORAGE_NAME = 'RevvitySignalsSearch';
export const QUERY_KEY = 'query';
export const PARAMS_KEY = 'params';

export const FIELDS_TO_EXCLUDE_FROM_CORPORATE_ID_WIDGET = ['eid', 'system.ExID', 'system.Keywords'];
export const ENTITY_FIELDS_TO_EXCLUDE = ['type', 'isTemplate', 'tags', 'eid', 'id'];
export const TAGS_TO_EXCLUDE = ['system.Keywords'];
export const USER_FIELDS = ['owner', 'createdBy', 'editedBy'];
export const FIRST_COL_NAMES = [MOL_COL_NAME, 'name', 'owner'];
export const LAST_COL_NAMES = ['createdAt', 'createdBy', 'modifiedAt', 'editedBy', 'description', 'ExID', HIDDEN_ID_COL_NAME]
export const FIELDS_SECTION_NAME = 'fields';
export const TABS_TO_EXCLUDE_FROM_WIDGET = ['chemicalDrawing'];
export const FIELDS_TO_EXCLUDE_FROM_WIDGET = ['assetTypeId', 'assetId', 'eid', 'type', 'library'];
export const MOLECULAR_FORMULA_FIELD_NAME = 'Molecular Formula';
export const SUBMITTER_FIELD_NAME = 'Submitter';

export const API_URL_PARAM_NAME = 'apiUrl';
export const API_KEY_PARAM_NAME = 'apiKey';