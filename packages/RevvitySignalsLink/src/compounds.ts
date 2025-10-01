import * as DG from 'datagrok-api/dg';
import { queryStructureById } from './revvity-api';
import { ComplexCondition, Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';

export const assetsQuery = {
    "query": {
        "$match": {
            "field": "type",
            "value": "asset",
            "mode": "keyword"
        }
    },
    "options": {
        "offset": 0,
        "limit": 20,
        "stop-after-items": 1000000
    }
};

export const batchesQuery = {
    "query": {
        "$match": {
            "field": "type",
            "value": "batch",
            "mode": "keyword"
        }
    },
    "options": {
        "offset": 0,
        "limit": 20,
        "stop-after-items": 1000000
    }
};

export const retrieveQueriesMap: {[key: string]: any} = {
    'asset': assetsQuery,
    'batch': batchesQuery,
}

export const materialsCondition: ComplexCondition = {
    logicalOperator: Operators.Logical.and,
    conditions: [
        {
            field: "isMaterial",
            operator: Operators.EQ,
            value: true
        },
        {
            field: "type",
            operator: Operators.NOT_EQ,
            value: "assetType"
        }
    ]
}

export async function addMoleculeStructures(moleculeIds: string[], molCol: DG.Column): Promise<void> {  
    const promises = moleculeIds.map((id, index) => 
    queryStructureById(id)
      .then(molecule => molCol.set(index, molecule))
      .catch(() => {})
  );

  await Promise.allSettled(promises);
}