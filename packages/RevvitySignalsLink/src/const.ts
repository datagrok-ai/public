export const testFilterCondition = {

    "logicalOperator": "and",
    "conditions": [
        {
            "field": "MW",
            "operator": ">",
            "value": 1.2345
        },
        {
            "logicalOperator": "or",
            "conditions": [
                {
                    "field": "Editor",
                    "operator": "=",
                    "value": "id123"
                },
                {
                    "field": "Creator",
                    "operator": "=",
                    "value": "User1"
                },
               
            ]
        },
        {
            "field": "MW",
            "operator": "<",
            "value": 1234567
        }
    ]
}