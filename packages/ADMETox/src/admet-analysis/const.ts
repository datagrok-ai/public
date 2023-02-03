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
        {
            "name": "PPB",
        },
        {
            "name": "VD",
        },
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
    "Excretion": {
        "name": "Excretion",
        "description": "Models to predict excretion",
        "models": [
        {
            "name": "Clearance",
        },
        {
            "name": "T",
        }]
    },
    "Toxicity": {
        "name": "Toxicity",
        "description": "Models to predict toxicity",
        "models": [
        {
            "name": "hERG",
        },
        {
            "name": "LD50",
        },
        {
            "name": "H-HT",
        },
        {
            "name": "Ames",
        },
        {
            "name": "SkinSen",
        }]
    },
    "Lipophilicity": {
        "name": "Lipophilicity",
        "description": "Lipophilicity",
        "models": [
        {
            "name": "logP",
        },
        {
            "name": "logD",
        }]
    },
    "Solubility": {
        "name": "Solubility",
        "description": "Solubility",
        "models": [
        {
            "mame": "logS",
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
    "F(20%)": {
        "0": "F20-", 
        "0.5": "F20+"
    },
    "F(30%)": {
        "0": "F30-", 
        "0.5": "F30+"
    },
    "PPB": {
        "0": "Significant with drugs that are highly protein-bound and have a low therapeutic index",
        "90": "Significant with drugs that are highly protein-bound and have a low therapeutic index"
    },
    "VD": {
        "0.07": "Evenly distributed",
        "0.7": "Bound to tissue components, highly lipophilic"
    },
    "BBB": {
        "0": "BB ratio <0.1: BBB-",
        "0.5": "BB ratio >=0.1: BBB+"
    },
    "CYP1A2-Inhibitor": {
        "0": "Non-inhibitor",
        "0.5": "Inhibitor"
    }, 
    "CYP1A2-Substrate": {
        "0": "Non-Substrate",
        "0.5": "Substrate"
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
    "CYP2C19-Substrate": {
        "0": "Non-Substrate", 
        "0.5": "Substrate"
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
    "Clearance": {
        "-100": "Low",
        "5": "Moderate",
        "15": "High"
    },
    "T": {
        "-100": "Low",
        "3": "Moderate",
        "8": "High"
    },
    "hERG": {
        "0": "Non-blockers",
        "0.5": "Blockers"
    },
    "LD50": {
        "1": "High-toxicity",
        "50": "Toxicity",
        "500": "Low-toxicity"
    },
    "H-HT": {
        "0": "H-HT negative(-)",
        "0.5": "H-HT positive(+)"
    },
    "Ames": {
        "0": "Ames negative(-)", 
        "0.5": "Ames positive(+)"
    },
    "SkinSen": {
        "0": "Non-sensitizer", 
        "0.5": "Sensitizer"
    },
    "logP": {
        "-100": "Poor lipid bilayer permeability",
        "0": "Optimal",
        "3": "Poor aqueous solubility"
    },
    "logD": {
        "-100": "Solubility high; Permeability low by passive transcellular diffusion; Permeability possible via paracellular if MW < 200; Metabolism low",
        "1": "Solubility moderate; Permeability moderate; Metabolism low",
        "3": "Solubility low; Permeability high; Metabolism moderate to high",
        "5": "Solubility low; Permeability high; Metabolism high"
    },
    "logS": {
        "-100": "Low solubility",
        "10": "Moderate solubility",
        "60": "High solubility"
    }
}