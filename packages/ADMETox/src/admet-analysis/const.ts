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
    "Pgp-Inhibitor": "Red: Non-inhibitor, Green: Inhibitor",
    "Pgp-Substrate": "Red: Non-Substrate, Green: Substrate",
    "HIA": "Red: HIA-, Green: HIA+",
    "F(20%)": "Red: F20-, Green: F20+",
    "F(30%)": "Red: F30-, Green: F30+",
    "BBB": "Red: BB ratio <0.1: BBB-, Green: BB ratio >=0.1: BBB+",
    "CYP1A2-Inhibitor": "Red: Non-inhibitor, Green: Inhibitor", 
    "CYP1A2-Substrate": "Red: Non-Substrate, Green: Substrate",
    "CYP3A4-Inhibitor": "Red: Non-inhibitor, Green: Inhibitor",
    "CYP3A4-Substrate": "Red: Non-Substrate, Green: Substrate", 
    "CYP2C19-Inhibitor": "Red: Non-inhibitor, Green: Inhibitor",
    "CYP2C19-Substrate": "Red: Non-Substrate, Green: Substrate",
    "CYP2C9-Inhibitor": "Red: Non-inhibitor, Green: Inhibitor",
    "CYP2C9-Substrate": "Red: Non-Substrate, Green: Substrate",
    "CYP2D6-Inhibitor": "Red: Non-inhibitor, Green: Inhibitor",
    "CYP2D6-Substrate": "Red: Non-Substrate, Green: Substrate",
    "Ames": "Red: Ames negative(-), Green: Ames positive(+)",
    "SkinSen": "Red:  Non-sensitizer, Green: Sensitizer",
}