import * as DG from 'datagrok-api/dg';
import { queryStructureById } from './revvity-api';
import { ComplexCondition, Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { SignalsSearchQuery } from './signals-search-query';

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

export function getConditionForLibAndType(type: string, assetTypeId: string, isMaterial: boolean): any[] {
    const innerAndConditions: any[] = [
        {
            "$match": {
                "field": "assetTypeEid",
                "value": assetTypeId,
            }
        },
        {
            "$match": {
                "field": "type",
                "value": type,
                "mode": "keyword"
            }
        },
    ];
    if (isMaterial) {
        innerAndConditions.push({
            "$and": [
                {
                    "$match": {
                        "field": "isMaterial",
                        "value": true
                    }
                },
                {
                    "$not": [
                        {
                            "$match": {
                                "field": "type",
                                "value": "assetType"
                            }
                        }
                    ]
                }
            ]
        })
    }
    return innerAndConditions;
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

    let counter = 0;
    const pb = DG.TaskBarProgressIndicator.create('Loading molecule structures...');
    const promises = moleculeIds.map((id, index) =>
        queryStructureById(id)
            .then(molecule => {
                counter++;
                molCol.set(index, molecule);
                pb.update((counter / molCol.length) * 100, `Loaded ${counter} of ${molCol.length} molecules`);
            })
            .catch(() => { })
    );

    try {
        await Promise.allSettled(promises);
    } finally {
        pb.close();
    }
}