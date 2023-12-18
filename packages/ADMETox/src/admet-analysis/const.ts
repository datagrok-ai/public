export const properties: any = {
    "Absorption": {
        "name": "Absorption",
        "description": "Models to predict absorbtion",
        "models": [
        {
            "skip": false,
            "name": "Pgp-Inhibitor"
        },
        {
            "skip": false,
            "name": "Pgp-Substrate"
        },
        {
            "skip": true,
            "name": "HIA"
        },
        {
            "skip": true,
            "name": "F(20%)"
        },
        {
            "skip": true,
            "name": "F(30%)"
        },
        {
            "skip": false,
            "name": "Caco2",
        },
        {
            "skip": false,
            "name": "Lipophilicity",
        },
        {
            "skip": false,
            "name": "Solubility",
        }]
    },
    "Distribution": {
        "name": "Distribution",
        "description": "Models to predict distribution",
        "models": [
        {
            "skip": false,
            "name": "PPBR"
        },
        {
            "skip": false,
            "name": "VDss"
        },
        {
            "skip": true,
            "name": "BBB"
        }]
    },
    "Metabolism": {
        "name": "Metabolism",
        "description": "Models to predict metabolism",
        "models": [
        {
            "skip": false,
            "name": "CYP1A2-Inhibitor"
        },
        {
            "skip": true,
            "name": "CYP1A2-Substrate"
        },
        {
            "skip": false,
            "name": "CYP3A4-Inhibitor"
        },
        {
            "skip": false,
            "name": "CYP3A4-Substrate"
        },
        {
            "skip": false,
            "name": "CYP2C19-Inhibitor"
        },
        {
            "skip": true,
            "name": "CYP2C19-Substrate"
        },
        {
            "skip": false,
            "name": "CYP2C9-Inhibitor"
        },
        {
            "skip": false,
            "name": "CYP2C9-Substrate"
        },
        {
            "skip": false,
            "name": "CYP2D6-Inhibitor"
        },
        {
            "skip": false,
            "name": "CYP2D6-Substrate"
        }]
    },
    "Excretion": {
        "name": "Excretion",
        "description": "Models to predict excretion",
        "models": [
        {
            "skip": false,
            "name": "CL-Hepa"
        },
        {
            "skip": false,
            "name": "CL-Micro"
        },
        {
            "skip": false,
            "name": "Half-Life"
        }]
    },
    "Toxicity": {
        "name": "Toxicity",
        "description": "Models to predict toxicity",
        "models": [
        {
            "skip": true,
            "name": "hERG"
        },
        {
            "skip": true,
            "name": "LD50"
        },
        {
            "skip": true,
            "name": "H-HT"
        },
        {
            "skip": true,
            "name": "Ames"
        },
        {
            "skip": true,
            "name": "SkinSen"
        }]
    },
    "Lipophilicity": {
        "name": "Lipophilicity",
        "description": "Lipophilicity",
        "models": [
        {
            "skip": true,
            "name": "logP"
        },
        {
            "skip": true,
            "name": "logD"
        }]
    },
    "Solubility": {
        "name": "Solubility",
        "description": "Solubility",
        "models": [
        {
            "skip": true,
            "name": "logS"
        }]
    }
}

export const models: any = {
    "Pgp-Inhibitor": {
        "0": "Non-inhibitor",
        "0.5": "Inhibitor"
    },
    "Pgp-Substrate": {
        "0": "Non-Substrate",
        "0.5": "Substrate"
    },
    "HIA": {
        "0": "HIA-", 
        "0.5": "HIA+"
    },
    "Caco2": {

    },
    "Lipophilicity": {

    },
    "Solubility": {
        
    },
    "PPBR": {
        "0": "Significant with drugs that are highly protein-bound and have a low therapeutic index",
        "90": "Significant with drugs that are highly protein-bound and have a low therapeutic index"
    },
    "VDss": {
        "0.07": "Evenly distributed",
        "0.7": "Bound to tissue components, highly lipophilic"
    },
    "CYP1A2-Inhibitor": {
        "0": "Non-inhibitor",
        "0.5": "Inhibitor"
    },
    "CYP3A4-Inhibitor": {
        "0": "Non-inhibitor", 
        "0.5": "Inhibitor"
    },
    "CYP3A4-Substrate": {
        "0": "Non-Substrate", 
        "0.5": "Substrate"
    }, 
    "CYP2C19-Inhibitor": {
        "0": "Non-inhibitor", 
        "0.5": "Inhibitor"
    },
    "CYP2C9-Inhibitor": {
        "0": "Non-inhibitor", 
        "0.5": "Inhibitor"
    },
    "CYP2C9-Substrate": {
        "0": "Non-Substrate", 
        "0.5": "Substrate"
    },
    "CYP2D6-Inhibitor": {
        "0": "Non-inhibitor", 
        "0.5": "Inhibitor"
    },
    "CYP2D6-Substrate": {
        "0": "Non-Substrate", 
        "0.5": "Substrate"
    },
    "CL-Hepa": {

    },
    "CL-Micro": {

    },
    "Half-Life": {
        "-100": "Low",
        "3": "Moderate",
        "8": "High"
    }
}