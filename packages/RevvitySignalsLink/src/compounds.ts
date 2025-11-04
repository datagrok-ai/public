import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { queryStructureById } from './revvity-api';
import { ComplexCondition, Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { SignalsSearchQuery } from './signals-search-query';
import { delay } from '@datagrok-libraries/utils/src/test';

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
    const BATCH_SIZE = 100;
    let counter = 0;
    const pb = DG.TaskBarProgressIndicator.create('Loading molecule structures...');
    
    // Split IDs into batches
    const batches: string[][] = [];
    for (let i = 0; i < moleculeIds.length; i += BATCH_SIZE) {
        batches.push(moleculeIds.slice(i, i + BATCH_SIZE));
    }

    const failedIdxs: number[] = [];

    try {
        // Process each batch sequentially
        for (let batchIndex = 0; batchIndex < batches.length; batchIndex++) {
            const batch = batches[batchIndex];
            const promises = batch.map(async (id, batchItemIndex) => {
                const globalIndex = batchIndex * BATCH_SIZE + batchItemIndex;
                try {
                    const molecule = await queryStructureById(id);
                    counter++;
                    molCol.set(globalIndex, molecule);
                    pb.update((counter / molCol.length) * 100, `Loaded ${counter} of ${molCol.length} molecules`);
                } catch (e) {
                    failedIdxs.push(globalIndex);
                }
            });
            
            await Promise.allSettled(promises);
            //make a delay after each batch to avoid 429 (Too many requests) error
            await delay(1000);
        }
    } catch (e: any) {
        throw e;
    } finally {
        pb.close();
        if (failedIdxs.length)
            grok.shell.warning(`Smiles for molecules ${failedIdxs.join(',')} failed to load`);
    }
}