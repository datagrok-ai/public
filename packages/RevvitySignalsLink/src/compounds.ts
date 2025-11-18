import * as DG from 'datagrok-api/dg';
import { queryStructureById } from './revvity-api';
import { ComplexCondition, Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { SignalsSearchQuery } from './signals-search-query';
import { delay } from '@datagrok-libraries/utils/src/test';
import { funcs } from './package-api';

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
    const pb = DG.TaskBarProgressIndicator.create('Loading molecule structures...', {cancelable: true});
    
    
    // Split IDs into batches
    const batches: string[][] = [];
    for (let i = 0; i < moleculeIds.length; i += BATCH_SIZE) {
        batches.push(moleculeIds.slice(i, i + BATCH_SIZE));
    }

    try {
        // Process each batch sequentially
        for (let batchIndex = 0; batchIndex < batches.length; batchIndex++) {
            const batch = batches[batchIndex];
            
            // Process all promises in this batch and collect results with their indices
            // Update progress bar as each promise completes, but maintain correct indexing
            const promises = batch.map(async (id, batchItemIndex) => {
                const globalIndex = batchIndex * BATCH_SIZE + batchItemIndex;
                try {
                    const molecule = await funcs.getStructureById(id);
                    // Update counter and progress bar as each promise completes
                    // Using post-increment to get the current count before incrementing
                    const currentCount = ++counter;
                    pb.update((currentCount / molCol.length) * 100, `Loaded ${currentCount} of ${molCol.length} molecules`);
                    return { index: globalIndex, molecule, success: true };
                } catch (e) {
                    console.log(`***************** 10 attemps failed for molecule ${globalIndex}`);
                    // Update counter even for failures
                    const currentCount = ++counter;
                    pb.update((currentCount / molCol.length) * 100, `Loaded ${currentCount} of ${molCol.length} molecules`);
                    return { index: globalIndex, molecule: null, success: false };
                }
            });
            
            // Wait for all promises in this batch to complete
            const results = await Promise.allSettled(promises);

            if (pb.canceled)
                return;
            
            // Process results and set molecules at correct indices
            // This ensures deterministic assignment regardless of promise completion order
            for (let resultIndex = 0; resultIndex < results.length; resultIndex++) {
                const result = results[resultIndex];
                if (result.status === 'fulfilled') {
                    const { index, molecule, success } = result.value;
                    if (success && molecule !== null) {
                        // Set molecule at the correct index - this is now deterministic
                        molCol.set(index, molecule);
                    }
                } else {
                    // Handle rejected promise (shouldn't happen with allSettled, but handle it)
                    const globalIndex = batchIndex * BATCH_SIZE + resultIndex;
                    console.log(`Failed to load molecule at index ${globalIndex}`);
                }
            }
            
            // Make a delay after each batch to avoid 429 (Too many requests) error
            if (batchIndex < batches.length - 1) {
                await delay(1000);
            }
        }
    } catch (e: any) {
        throw e;
    } finally {
        pb.close();
    }
}