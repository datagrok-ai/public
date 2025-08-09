// Signals /entities/search Query TypeScript Interfaces

import { RevvityLibrary } from "./package";
import { ComplexCondition, Operators } from "./query-builder";
import { SimpleCondition } from './query-builder';
import { PROPERTY_NAMES_TO_QUERY_MAPPING } from './properties';

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
 * Maps internal field names to Revvity API field names using the property mapping
 */
function mapFieldName(fieldName: string, explicitMapping?: string): string {
  return explicitMapping || PROPERTY_NAMES_TO_QUERY_MAPPING[fieldName as keyof typeof PROPERTY_NAMES_TO_QUERY_MAPPING] || fieldName;
}

type AsType = 'text' | 'date'

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
export function convertSimpleConditionToQueryOperator(condition: SimpleCondition, fieldName?: string): QueryOperator | QueryNotOperator {
  const { field, operator, value } = condition;
  
  // Apply field name mapping if applicable
  const targetField = mapFieldName(field, fieldName);
  
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
  
  // Find the corresponding Revvity operator
  const revvityOperator = REVVITY_OPERATORS_MAPPING[operator as keyof typeof REVVITY_OPERATORS_MAPPING];
  
  if (!revvityOperator) {
    throw new Error(`Unsupported operator: ${operator}`);
  }
  
  // Create the appropriate QueryOperator based on the Revvity operator
  switch (revvityOperator) {
    case OPERATORS.IN:
      if (!Array.isArray(value)) {
        throw new Error(`Operator ${operator} requires an array value`);
      }
             return {
         $in: {
           field: targetField,
           values: value
         }
       } as QueryInOperator;
      
         case OPERATORS.GT:
       const gtResult: any = {
         field: targetField,
         value: convertDateValue(value)
       };
       if (isDateValue(value))
         gtResult.as = 'date';
       return { $gt: gtResult } as QueryGtOperator;
       
     case OPERATORS.LT:
       const ltResult: any = {
         field: targetField,
         value: convertDateValue(value)
       };
       if (isDateValue(value))
         ltResult.as = 'date';
       return { $lt: ltResult } as QueryLtOperator;
      
         case OPERATORS.GTE:
       const gteResult: any = {
         field: targetField,
         value: convertDateValue(value)
       };
       if (isDateValue(value))
         gteResult.as = 'date';
       return { $gte: gteResult } as QueryGteOperator;
       
     case OPERATORS.LTE:
       const lteResult: any = {
         field: targetField,
         value: convertDateValue(value)
       };
       if (isDateValue(value))
         lteResult.as = 'date';
       return { $lte: lteResult } as QueryLteOperator;
      
         case OPERATORS.RANGE:
       if (!Array.isArray(value)) {
         throw new Error(`Operator ${operator} requires an array`);
       }
       const rangeResult: any = {
         field: targetField,
         from: convertDateValue(value[0]),
         to: convertDateValue(value[1])
       };
       if (isDateValue(value[0]) || isDateValue(value[1]))
         rangeResult.as = 'date';
       return { $range: rangeResult } as QueryRangeOperator;
      
    case OPERATORS.MATCH:
             return {
         $match: {
           field: targetField,
           value
         }
       } as QueryMatchOperator;
      
    case OPERATORS.PREFIX:
             return {
         $prefix: {
           field: targetField,
           value,
           mode: 'keyword'
         }
       } as QueryPrefixOperator;
      
    case OPERATORS.EXISTS:
             return {
         $exists: {
           field: targetField
         }
       } as QueryExistsOperator;

    case OPERATORS.NOT:
      return {
        $not:
          [{
            $match: {
              field: targetField,
              value
            }
          }]
      } as QueryNotOperator;

    case OPERATORS.CHEMSEARCH:
      const chemsearchResult: any = {
        molecule: value.molecule ?? value,
        mime: 'chemical/x-daylight-smiles'
      }
        ;
      if (value.threshold)
        chemsearchResult.options = `full=no,similar=yes,simthreshold=${value.threshold * 100}`
      return { $chemsearch: chemsearchResult } as QueryChemSearchOperator;

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

