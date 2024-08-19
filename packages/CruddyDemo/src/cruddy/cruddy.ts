import * as DG from 'datagrok-api/dg';

export type TFilter = {[name: string]: string | null};



export interface IFilterDescription {
  type: 'distinct' | 'combo' | 'radio' | 'search' | 'range' | 'expression';
  column: string;
}

//type GridColumnNames = { gridColumnsNames?: string[]; }


export namespace Exp {
  export const CONTAINS = 'contains';
  export const STARTS_WITH = 'starts with';
  export const ENDS_WITH = 'ends with';
  export const GT = '>';
  export const LT = '<';

  export const typeOperators: {[key: string]: string[]} = {
    [DG.TYPE.STRING]: [ STARTS_WITH, ENDS_WITH ],
    [DG.TYPE.INT]: [ GT, LT ],
    [DG.TYPE.FLOAT]: [ GT, LT ],
  };

  export function getCondition(column: string, op: string, value: any): TFilter {
    if (!value || value === '')
      return {};

    switch (op) {
      case STARTS_WITH: return {[`starts_with(${column}, '${value}')`]: null};
      case ENDS_WITH: return {[`${column} like '%${value}'`]: null};
      case LT: return {[`${column} < ${value}`]: null};
      case GT: return {[`${column} > ${value}`]: null};
    }
    return {};
  }
}