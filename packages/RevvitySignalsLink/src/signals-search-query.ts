// Signals /entities/search Query TypeScript Interfaces

import { ComplexCondition, Operators, SimpleCondition } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { NOT_IN_TAGS } from './properties';
import { RevvityLibrary } from './libraries';

export enum OPERATORS {
  MATCH = '$match',
  IN = '$in',
  GT = '$gt',
  LT = '$lt',
  GTE = '$gte',
  LTE = '$lte',
  RANGE = '$range',
  PREFIX = '$prefix',
  EXISTS = '$exists',
  SIMPLE = '$simple',
  INTERSECT = '$intersect',
  CHEMSEARCH = '$chemsearch',
  AND = '$and',
  OR = '$or',
  NOT = '$not',
}

const REVVITY_OPERATORS_MAPPING = {
  [Operators.EQ]: OPERATORS.MATCH,
  [Operators.IN]: OPERATORS.IN,
  [Operators.GT]: OPERATORS.GT,
  [Operators.LT]: OPERATORS.LT,
  [Operators.GTE]: OPERATORS.GTE,
  [Operators.LTE]: OPERATORS.LTE,
  [Operators.BETWEEN]: OPERATORS.RANGE,
  [Operators.STARTS_WITH]: OPERATORS.PREFIX,
  [Operators.Logical.and]: OPERATORS.AND,
  [Operators.Logical.or]: OPERATORS.OR,
  [Operators.NOT_EQ]: OPERATORS.NOT,
  [Operators.BEFORE]: OPERATORS.LT,
  [Operators.AFTER]: OPERATORS.GT,
  [Operators.CONTAINS]: OPERATORS.CHEMSEARCH,
  [Operators.IS_SIMILAR]: OPERATORS.CHEMSEARCH,
}


/**
 * Checks if a field should include the "in": "tags" property
 * Returns true if the field name is NOT in REVVITY_DEFAULT_FILTERS_NAMES
 */
function shouldIncludeTagsProperty(fieldName: string): boolean {
  return !NOT_IN_TAGS.includes(fieldName);
}

type AsType = 'text' | 'date' | 'double'

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
  $gt: { field: string; in?: 'tags'; as: AsType; value: number | string };
}
export interface QueryLtOperator {
  $lt: { field: string; in?: 'tags'; as: AsType; value: number | string };
}
export interface QueryGteOperator {
  $gte: { field: string; in?: 'tags'; as: AsType; value: number | string };
}
export interface QueryLteOperator {
  $lte: { field: string; in?: 'tags'; as: AsType; value: number | string };
}

// $range
export interface QueryRangeOperator {
  $range: { field: string; in?: 'tags'; as: AsType; from: number | string; to: number | string };
}

// $match, $prefix
export interface QueryMatchOperator {
  $match: { field: string; in?: 'tags'; value: string | boolean, mode?: string };
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
  field?: string;
  in?: string;
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
  PARAEXP = 'paraexp',
  PARASUBEXP = 'parasubexp',
  SAMPLE_CONTAINER = 'sampleContainer',
  WORKSHEET = 'worksheet',
}

/**
 * Converts a SimpleCondition to the appropriate QueryOperator based on the operator mapping
 */
export function convertSimpleConditionToQueryOperator(condition: SimpleCondition): QueryOperator | QueryNotOperator {
  const { field, operator, value } = condition;
  
  
  // Helper function to convert dayjs values to string
  const convertDateValue = (val: any): string => {
    if (val && typeof val === 'object' && val.format) {
      // It's a dayjs object
      return val.format('YYYY-MM-DDTHH:mm:ss.SSS[Z]');
    }
    return val;
  };
  
  // Helper function to check if a value is a dayjs object
  const isDateValue = (val: any): boolean => {
    return val && typeof val === 'object' && val.format;
  };
  
  // Helper function to determine the 'as' type based on value
  const getAsType = (val: any): AsType => {
    if (isDateValue(val)) {
      return 'date';
    } else if (typeof val === 'number') {
      return 'double';
    } else if (typeof val === 'string') {
      return 'text';
    } else {
      // Default to text for other types
      return 'text';
    }
  };
  
  // Find the corresponding Revvity operator
  const revvityOperator = REVVITY_OPERATORS_MAPPING[operator as keyof typeof REVVITY_OPERATORS_MAPPING];
  
  if (!revvityOperator) {
    throw new Error(`Unsupported operator: ${operator}`);
  }
  
  // Create the appropriate QueryOperator based on the Revvity operator
  let result: any;
  
  switch (revvityOperator) {
    case OPERATORS.IN:
      if (!Array.isArray(value)) {
        throw new Error(`Operator ${operator} requires an array value`);
      }
      result = {
        field: field,
        values: value
      };
      break;
      
    case OPERATORS.GT:
      result = {
        field: field,
        value: convertDateValue(value),
        as: getAsType(value)
      };
      break;
       
    case OPERATORS.LT:
      result = {
        field: field,
        value: convertDateValue(value),
        as: getAsType(value)
      };
      break;
      
    case OPERATORS.GTE:
      result = {
        field: field,
        value: convertDateValue(value),
        as: getAsType(value)
      };
      break;
       
    case OPERATORS.LTE:
      result = {
        field: field,
        value: convertDateValue(value),
        as: getAsType(value)
      };
      break;
      
    case OPERATORS.RANGE:
      if (!Array.isArray(value)) {
        throw new Error(`Operator ${operator} requires an array`);
      }
      result = {
        field: field,
        from: convertDateValue(value[0]),
        to: convertDateValue(value[1]),
        as: getAsType(value[0]) // Use the first value to determine type
      };
      break;
      
    case OPERATORS.MATCH:
      result = {
        field: field,
        value
      };
      if (typeof field === 'string') {
        result.mode = 'keyword'
      }
      break;
      
    case OPERATORS.PREFIX:
      result = {
        field: field,
        value,
        mode: 'keyword'
      };
      break;
      
    case OPERATORS.EXISTS:
      result = {
        field: field
      };
      break;

    case OPERATORS.INTERSECT:
      if (!Array.isArray(value)) {
        throw new Error(`Operator ${operator} requires an array value`);
      }
      result = {
        field: field,
        values: value
      };
      break;

    case OPERATORS.NOT:
      const notMatchResult = {
        field: field,
        value
      };
      if (shouldIncludeTagsProperty(field)) {
        (notMatchResult as any).in = 'tags';
      }
      return {
        $not:
          [{
            $match: notMatchResult
          }]
      } as QueryNotOperator;

    case OPERATORS.CHEMSEARCH:
      const chemsearchResult: any = {
        molecule: value.molecule ?? value,
        mime: 'chemical/x-daylight-smiles'
      };
      if (value.threshold)
        chemsearchResult.options = `full=no,similar=yes,simthreshold=${value.threshold * 100}`;
      return { $chemsearch: chemsearchResult } as QueryChemSearchOperator;

    default:
      throw new Error(`Unsupported Revvity operator: ${revvityOperator}`);
  }
  
  // Add tags property if needed (for all operators except NOT and CHEMSEARCH)
  if (shouldIncludeTagsProperty(field)) {
    (result as any).in = 'tags';
  }
  
  // Return the appropriate operator with the result
  switch (revvityOperator) {
    case OPERATORS.IN:
      return { $in: result } as QueryInOperator;
    case OPERATORS.GT:
      return { $gt: result } as QueryGtOperator;
    case OPERATORS.LT:
      return { $lt: result } as QueryLtOperator;
    case OPERATORS.GTE:
      return { $gte: result } as QueryGteOperator;
    case OPERATORS.LTE:
      return { $lte: result } as QueryLteOperator;
    case OPERATORS.RANGE:
      return { $range: result } as QueryRangeOperator;
    case OPERATORS.MATCH:
      return { $match: result } as QueryMatchOperator;
    case OPERATORS.PREFIX:
      return { $prefix: result } as QueryPrefixOperator;
    case OPERATORS.EXISTS:
      return { $exists: result } as QueryExistsOperator;
    case OPERATORS.INTERSECT:
      return { $intersect: result } as QueryIntersectOperator;
    default:
      throw new Error(`Unsupported Revvity operator: ${revvityOperator}`);
  }
}


/**
 * Converts a ComplexCondition to QueryOperator recursively
 */
export function convertComplexConditionToQueryOperator(condition: ComplexCondition): QueryOperator {
  if (!condition.conditions || condition.conditions.length === 0) {
    throw new Error('ComplexCondition must have at least one condition');
  }

  // Convert all nested conditions to QueryOperators
  const queryOperators: QueryOperator[] = [];
  
  for (const nestedCondition of condition.conditions) {
    if ('field' in nestedCondition) {
      // This is a SimpleCondition
      const queryOperator = convertSimpleConditionToQueryOperator(nestedCondition);
      queryOperators.push(queryOperator);
    } else {
      // This is a nested ComplexCondition
      const nestedQueryOperator = convertComplexConditionToQueryOperator(nestedCondition);
      queryOperators.push(nestedQueryOperator);
    }
  }

  // If there's only one condition, return it directly (no need for logical wrapper)
  if (queryOperators.length === 1) {
    return queryOperators[0];
  }

  // Create logical operator based on the condition's logical operator
  const logicalOp = condition.logicalOperator === Operators.Logical.and ? OPERATORS.AND : OPERATORS.OR;
  
  switch (logicalOp) {
    case OPERATORS.AND:
      return { $and: queryOperators } as QueryAndOperator;
    case OPERATORS.OR:
      return { $or: queryOperators } as QueryOrOperator;
    default:
      throw new Error(`Unsupported logical operator: ${condition.logicalOperator}`);
  }
}

/**
 * Converts a ComplexCondition to a complete SignalsSearchQuery
 */
export function convertComplexConditionToSignalsSearchQuery(
  condition: ComplexCondition,
  options?: {
    sort?: { [field: string]: 'asc' | 'desc' };
    [key: string]: any;
  }
): SignalsSearchQuery {
  const queryOperator = convertComplexConditionToQueryOperator(condition);
  
  return {
    query: queryOperator,
    options
  };
}

export function createMaterialsSearch(library: RevvityLibrary, condition: ComplexCondition) {

}

