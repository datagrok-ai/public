export const testFilterCondition = {

    "operator": "and",
    "conditions": [
        {
            "field": "Created",
            "operator": "after",
            "values": [
                "2023-10-27T10:00:00.000Z",
                null
            ]
        },
        {
            "operator": "or",
            "conditions": [
                {
                    "field": "Id",
                    "operator": "=",
                    "values": [
                        "id123"
                    ]
                },
                {
                    "field": "Created by",
                    "operator": "=",
                    "values": [
                        "User1"
                    ]
                },
                {
                    "field": "Created",
                    "operator": "between",
                    "values": [
                        null,
                        null,
                    ]
                },
            ]
        },
        {
            "field": "MW",
            "operator": "<",
            "values": [
                1234567
            ]
        }
    ]
}