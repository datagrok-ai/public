// Signals /entities/search Query TypeScript Interfaces

// Base for all field-based operators
interface FieldOperatorBase {
  field: string;
  in?: 'tags'; // Extend as needed
}

// $in operator
export interface QueryInOperator {
  $in: {
    field: string;
    in?: 'tags';
    values: (string | number | boolean)[];
  };
}

// $gt, $lt, $gte, $lte
export interface QueryGtOperator {
  $gt: { field: string; in?: 'tags'; value: number };
}
export interface QueryLtOperator {
  $lt: { field: string; in?: 'tags'; value: number };
}
export interface QueryGteOperator {
  $gte: { field: string; in?: 'tags'; value: number };
}
export interface QueryLteOperator {
  $lte: { field: string; in?: 'tags'; value: number };
}

// $range
export interface QueryRangeOperator {
  $range: { field: string; in?: 'tags'; from: number; to: number };
}

// $match, $prefix
export interface QueryMatchOperator {
  $match: { field: string; in?: 'tags'; value: string, mode?: string };
}
export interface QueryPrefixOperator {
  $prefix: { field: string; in?: 'tags'; value: string, mode: string };
}

// $exists
export interface QueryExistsOperator {
  $exists: { field: string; in?: 'tags' };
}

// $simple
export interface QuerySimpleOperator {
  $simple: { query: string };
}

// $intersect
export interface QueryIntersectOperator {
  $intersect: {
    field: string;
    in?: 'tags';
    values: (string | number | boolean)[];
    minimum?: number;
  };
}

// $chemsearch
export interface QueryChemSearchOperator {
  $chemsearch: {
    molecule: string;
    mime: string;
    options?: string;
  };
}

// Logical operators (recursive)
export interface QueryOrOperator {
  $or: QueryOperator[];
}
export interface QueryAndOperator {
  $and: QueryOperator[];
}
export interface QueryNotOperator {
  $not: QueryOperator[];
}

// Union of all possible operators
export type QueryOperator =
  | QueryInOperator
  | QueryGtOperator
  | QueryLtOperator
  | QueryGteOperator
  | QueryLteOperator
  | QueryRangeOperator
  | QueryMatchOperator
  | QueryPrefixOperator
  | QueryExistsOperator
  | QuerySimpleOperator
  | QueryIntersectOperator
  | QueryChemSearchOperator
  | QueryOrOperator
  | QueryAndOperator
  | QueryNotOperator;

// The full query object
export interface SignalsSearchQuery {
  query: QueryOperator;
  options?: {
    sort?: { [field: string]: 'asc' | 'desc' };
    // ...other options
  };
} 

export interface SignalsSearchParams {
  'page[offset]'?: number;
  'page[limit]'?: number;
  include?: string[];
  sort?: string;
  source?: string;
  'fields[entity]'?: string;
  stopAfterItems?: number;
} 

export enum SignalsSearchField {
  EID = 'eid',
  IS_TEMPLATE = 'isTemplate',
  NAME = 'name',
  CREATED_AT = 'createdAt',
  CREATED_BY = 'createdBy',
  MODIFIED_AT = 'modifiedAt',
  EDITED_BY = 'editedBy',
  OWNER = 'owner',
  STATE = 'state',
  TYPE = 'type',
  NOTEBOOK_EID = 'notebookEid',
  EXPERIMENT_EID = 'experimentEid',
  REACTANT = 'reactant',
  PRODUCT = 'product',
  FORMULA = 'formula',
  SOLVENT = 'solvent',
  FM = 'fm',
  MW = 'mw',
  EM = 'em',
  YIELD = 'yield',
  PURITY = 'purity',
  DURATION = 'duration',
  PRESSURE = 'pressure',
  ACTUAL = 'actual',
  SAMPLE = 'sample',
  MOLARITY = 'molarity',
  VOLUME = 'volume',
  TEMPERATURE = 'temperature',
  SUPPLIER = 'supplier',
  LOT_NUMBER = 'lot_number',
  LOAD = 'load',
  STATE_AT = 'stateAt',
  REVIEWERS = 'reviewers',
  SIGNED_BY = 'signedBy',
  SIGNED_AT = 'signedAt',
  ASSET_TYPE_EID = 'assetTypeEid',
}

export const SignalsSearchFieldDescriptions: Record<SignalsSearchField, string> = {
  [SignalsSearchField.EID]: 'ID',
  [SignalsSearchField.IS_TEMPLATE]: 'Is Template',
  [SignalsSearchField.NAME]: 'Name',
  [SignalsSearchField.CREATED_AT]: 'Creation time',
  [SignalsSearchField.CREATED_BY]: 'Creator ID',
  [SignalsSearchField.MODIFIED_AT]: 'Latest Editing time',
  [SignalsSearchField.EDITED_BY]: 'Latest Editor ID',
  [SignalsSearchField.OWNER]: 'Owner ID',
  [SignalsSearchField.STATE]: 'Experiment state (open: Active, closed: Closed, in_review: In Review, rejected: Rejected, voided: Voided)',
  [SignalsSearchField.TYPE]: 'Entity type (experiment, request, task, journal, chemicalDrawing, grid, materialsTable, acronym, text, pdf, viewonly, presentation, excel, imageResource, uploadedResource, spotfiredxp, sample, assetType, asset, batch)',
  [SignalsSearchField.NOTEBOOK_EID]: 'Parent Notebook ID',
  [SignalsSearchField.EXPERIMENT_EID]: 'Parent Experiment ID',
  [SignalsSearchField.REACTANT]: 'Stoic: Reactant',
  [SignalsSearchField.PRODUCT]: 'Stoic: Product',
  [SignalsSearchField.FORMULA]: 'Stoic: MF',
  [SignalsSearchField.SOLVENT]: 'Stoic: Solvent',
  [SignalsSearchField.FM]: 'Stoic: FM',
  [SignalsSearchField.MW]: 'Stoic: MW',
  [SignalsSearchField.EM]: 'Stoic: EM',
  [SignalsSearchField.YIELD]: 'Stoic: Yield',
  [SignalsSearchField.PURITY]: 'Stoic: Purity',
  [SignalsSearchField.DURATION]: 'Stoic: Reaction Time',
  [SignalsSearchField.PRESSURE]: 'Stoic: Pressure',
  [SignalsSearchField.ACTUAL]: 'Stoic: Actual Mass',
  [SignalsSearchField.SAMPLE]: 'Stoic: Sample Mass',
  [SignalsSearchField.MOLARITY]: 'Stoic: Reaction Molarity',
  [SignalsSearchField.VOLUME]: 'Stoic: Total Volume',
  [SignalsSearchField.TEMPERATURE]: 'Stoic: Temperature',
  [SignalsSearchField.SUPPLIER]: 'Stoic: Supplier',
  [SignalsSearchField.LOT_NUMBER]: 'Stoic: Lot Number',
  [SignalsSearchField.LOAD]: 'Stoic: Number',
  [SignalsSearchField.STATE_AT]: 'Experiment: State Modification Time',
  [SignalsSearchField.REVIEWERS]: 'Experiment: Reviewer IDs',
  [SignalsSearchField.SIGNED_BY]: 'Experiment: Signer IDs',
  [SignalsSearchField.SIGNED_AT]: 'Experiment: Latest Signing Time',
  [SignalsSearchField.ASSET_TYPE_EID]: 'Material: Library Name',
}; 

// Entity type values for the 'type' field
// experiment: Experiment
// request: Request
// task: Task
// journal: Notebook
// chemicalDrawing: Chemical Drawing
// grid: Admin Defined Table
// materialsTable: Materials Table
// acronym: Chemical Acronym
// text: Text
// pdf: PDF Document
// viewonly: Word Document
// presentation: PowerPoint Document
// excel: Excel Document
// imageResource: Image
// uploadedResource: Uploaded File
// spotfiredxp: Spotfire DXP Document
// sample: Sample
// assetType: Material Library
// asset: Material Asset
// batch: Material Batch
export enum SignalsEntityType {
  EXPERIMENT = 'experiment',
  REQUEST = 'request',
  TASK = 'task',
  JOURNAL = 'journal',
  CHEMICAL_DRAWING = 'chemicalDrawing',
  GRID = 'grid',
  MATERIALS_TABLE = 'materialsTable',
  ACRONYM = 'acronym',
  TEXT = 'text',
  PDF = 'pdf',
  VIEWONLY = 'viewonly',
  PRESENTATION = 'presentation',
  EXCEL = 'excel',
  IMAGE_RESOURCE = 'imageResource',
  UPLOADED_RESOURCE = 'uploadedResource',
  SPOTFIREDXP = 'spotfiredxp',
  SAMPLE = 'sample',
  ASSET_TYPE = 'assetType',
  ASSET = 'asset',
  BATCH = 'batch',
} 