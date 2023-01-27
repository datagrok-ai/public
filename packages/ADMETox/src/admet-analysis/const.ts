export const properties: any = {
    "Absorption": {
        "name": "Absorption",
        "description": "Models to predict absorbtion",
        "models": [
        /*{
            "name": "Caco-2",
        },*/
        {
            "name": "Pgp-Inhibitor",
        },
        {
            "name": "Pgp-Substrate",
        },
        {
            "name": "HIA",
        },
        {
            "name": "F(20%)",
        },
        {
            "name": "F(30%)",
        }]
    },
    "Distribution": {
        "name": "Distribution",
        "description": "Models to predict distribution",
        "models": [
        /*{
            "name": "PPB",
        },
        {
            "name": "VD",
        },*/
        {
            "name": "BBB",
        }]
    },
    "Metabolism": {
        "name": "Metabolism",
        "description": "Models to predict metabolism",
        "models": [
        {
            "name": "CYP1A2-Inhibitor",
        },
        {
            "name": "CYP1A2-Substrate",
        },
        {
            "name": "CYP3A4-Inhibitor",
        },
        {
            "name": "CYP3A4-Substrate",
        },
        {
            "name": "CYP2C19-Inhibitor",
        },
        {
            "name": "CYP2C19-Substrate",
        },
        {
            "name": "CYP2C9-Inhibitor",
        },
        {
            "name": "CYP2C9-Substrate",
        },
        {
            "name": "CYP2D6-Inhibitor",
        },
        {
            "name": "CYP2D6-Substrate",
        }]
    },
    /*"Excretion": {
        "name": "Excretion",
        "description": "Models to predict excretion",
        "models": [
        {
            "name": "Clearance",
        },
        {
            "name": "T1/2",
        }]
    },*/
    "Toxicity": {
        "name": "Toxicity",
        "description": "Models to predict toxicity",
        "models": [
        /*{
            "name": "hERG",
        },
        {
            "name": "LD50 of acute toxicity",
        }
        {
            "name": "H-HT",
        },*/
        {
            "name": "Ames",
        },
        {
            "name": "SkinSen",
        }]
    }
}

export const models: any = {
    "Pgp-Inhibitor": "<= 0.5: Non-inhibitor, > 0.5: Inhibitor",
    "Pgp-Substrate": "<= 0.5: Non-Substrate, > 0.5: Substrate",
    "HIA": "<= 0.5: HIA-, > 0.5: HIA+",
    "F(20%)": "<= 0.5: F20-, > 0.5: F20+",
    "F(30%)": "<= 0.5: F30-, > 0.5: F30+",
    "BBB": "<= 0.5: BB ratio <0.1: BBB-, > 0.5: BB ratio >=0.1: BBB+",
    "CYP1A2-Inhibitor": "<= 0.5: Non-inhibitor, > 0.5: Inhibitor", 
    "CYP1A2-Substrate": "<= 0.5: Non-Substrate, > 0.5: Substrate",
    "CYP3A4-Inhibitor": "<= 0.5: Non-inhibitor, > 0.5: Inhibitor",
    "CYP3A4-Substrate": "<= 0.5: Non-Substrate, > 0.5: Substrate", 
    "CYP2C19-Inhibitor": "<= 0.5: Non-inhibitor, > 0.5: Inhibitor",
    "CYP2C19-Substrate": "<= 0.5: Non-Substrate, > 0.5: Substrate",
    "CYP2C9-Inhibitor": "<= 0.5: Non-inhibitor, > 0.5: Inhibitor",
    "CYP2C9-Substrate": "<= 0.5: Non-Substrate, > 0.5: Substrate",
    "CYP2D6-Inhibitor": "<= 0.5: Non-inhibitor, > 0.5: Inhibitor",
    "CYP2D6-Substrate": "<= 0.5: Non-Substrate, > 0.5: Substrate",
    "Ames": "<= 0.5: Ames negative(-), > 0.5: Ames positive(+)",
    "SkinSen": "<= 0.5:  Non-sensitizer, > 0.5: Sensitizer",
}