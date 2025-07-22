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