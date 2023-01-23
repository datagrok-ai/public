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