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

export const STORAGE_NAME = 'RevvitySignalsSearch';
export const QUERY_KEY = 'query';
export const PARAMS_KEY = 'params';

export const ENTITY_FIELDS_TO_EXCLUDE = ['type', 'isTemplate', 'tags', 'eid'];
export const TAGS_TO_EXCLUDE = ['system.Keywords'];
export const USER_FIELDS = ['owner', 'createdBy', 'editedBy'];
export const FIRST_COL_NAMES = [MOL_COL_NAME, 'id', 'name'];
export const FIELDS_SECTION_NAME = 'fields';
export const TABS_TO_EXCLUDE_FROM_WIDGET = ['chemicalDrawing'];
export const FIELDS_TO_EXCLUDE_FROM_WIDGET = ['assetTypeId', 'assetId', 'eid', 'type', 'library'];
export const MOLECULAR_FORMULA_FIELD_NAME = 'Molecular Formula';
export const SUBMITTER_FIELD_NAME = 'Submitter';

export const API_URL_PARAM_NAME = 'apiUrl';
export const API_KEY_PARAM_NAME = 'apiKey';