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
  $match: { field: string; in?: 'tags'; value: string };
}
export interface QueryPrefixOperator {
  $prefix: { field: string; in?: 'tags'; value: string };
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