import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { queryStructureById } from './revvity-api';

export const assetsQuery = {
    "query": {
        "$and": [
            {
                "$match": {
                    "field": "assetTypeEid",
                    "value": "assetType:686ecf60e3c7095c954bd94f",
                    "mode": "keyword"
                }
            },
            {
                "$match": {
                    "field": "type",
                    "value": "asset",
                    "mode": "keyword"
                }
            },
            {
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
            }
        ]
    },
    "options": {
        "offset": 0,
        "limit": 20,
        "stop-after-items": 1000000
    },
    "meta": {
        "reason": "Advanced Search"
    }
}

export const batchesQuery = {
    "query": {
        "$and": [
            {
                "$match": {
                    "field": "assetTypeEid",
                    "value": "assetType:686ecf60e3c7095c954bd94f",
                    "mode": "keyword"
                }
            },
            {
                "$match": {
                    "field": "type",
                    "value": "batch",
                    "mode": "keyword"
                }
            },
            {
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
            }
        ]
    },
    "options": {
        "offset": 0,
        "limit": 20,
        "stop-after-items": 1000000
    },
    "meta": {
        "reason": "Advanced Search"
    }
}

export const MOL_COL_NAME = 'molecule';

export async function addMoleculeStructures(moleculeIds: string[], molCol: DG.Column): Promise<void> {  
    const promises = moleculeIds.map((id, index) => 
    queryStructureById(id)
      .then(molecule => molCol.set(index, molecule))
      .catch(() => {})
  );

  await Promise.allSettled(promises);
}