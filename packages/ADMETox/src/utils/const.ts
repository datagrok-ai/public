import * as DG from 'datagrok-api/dg';

export const properties: any = {
  "Absorption": {
  },
  "Distribution": {
    "name": "Distribution",
    "description": "Models to predict distribution",
    "models": [
    {
      "skip": false,
      "specific": false,
      "units": "%",
      "interpretation": {
        "range": {
          "below_40": "Low",
          "to_40_90": "Moderate",
          "above_90": "High",
        }
      },
      "coloring": {
        "type": "Conditional",
        "> 90": `${DG.Color.toHtml(DG.Color.cyan)}`,
        "40-90": `${DG.Color.toHtml(DG.Color.maroon)}`,
        "< 40": `${DG.Color.toHtml(DG.Color.yellow)}`,
      },
      "preference": "higher",
      "name": "PPBR"
    },
    {
      "skip": false,
      "specific": false,
      "units": "L/kg",
      "interpretation": {
        "range": {
          "below_1": "Low",
          "to_1_2": "Moderate",
          "above_2": "High",
        }
      },
      "coloring": {
        "type": "Conditional",
        "> 2": `${DG.Color.toHtml(DG.Color.cyan)}`,
        "1-2": `${DG.Color.toHtml(DG.Color.maroon)}`,
        "< 1": `${DG.Color.toHtml(DG.Color.yellow)}`,
      },
      "preference": "higher",
      "name": "VDss"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "",
      "coloring": "",
      "preference": "",
      "name": "BBB"
    }]
  },
  "Metabolism": {
    "name": "Metabolism",
    "description": "Models to predict metabolism",
    "models": [
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP1A2-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP1A2-Substrate"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP3A4-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP3A4-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP2C19-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP2C19-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP2C9-Inhibitor"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP2C9-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP2D6-Inhibitor"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "CYP2D6-Substrate"
    }]
  },
  "Excretion": {
    "name": "Excretion",
    "description": "Models to predict excretion",
    "models": [
    {
      "skip": false,
      "specific": true,
      "units": "mL/min",
      "interpretation": {
        "range": {
          "below_10": "Low",
          "to_10_100": "Moderate",
          "above_100": "Bad",
        }
      },
      "coloring": {
        "type": "Linear",
        "min": "10",
        "max": "100",
        "colors": `[${DG.Color.green}, ${DG.Color.red}]`,
      },
      "preference": "higher",
      "name": "CL-Hepa"
    },
    {
      "skip": false,
      "specific": true,
      "units": "μL/min/million cells",
      "interpretation": {
        "range": {
          "below_10": "Low",
          "to_10_100": "Moderate",
          "above_100": "Bad",
        }
      },
      "coloring": {
        "type": "Linear",
        "min": "10",
        "max": "100",
        "colors": `[${DG.Color.green}, ${DG.Color.red}]`,
      },
      "preference": "higher",
      "name": "CL-Micro"
    },
    {
      "skip": false,
      "specific": true,
      "units": "hours",
      "interpretation": "No good or bad values",
      "coloring": "",
      "preference": "higher",
      "name": "Half-Life"
    }]
  },
  "Toxicity": {
    "name": "Toxicity",
    "description": "Models to predict toxicity",
    "models": [
    {
      "skip": true,
      "specific": true,
      "name": "hERG"
    },
    {
      "skip": true,
      "specific": true,
      "name": "LD50"
    },
    {
      "skip": true,
      "specific": true,
      "name": "H-HT"
    },
    {
      "skip": true,
      "specific": true,
      "name": "Ames"
    },
    {
      "skip": true,
      "specific": true,
      "name": "SkinSen"
    }]
  },
  "Lipophilicity": {
    "name": "Lipophilicity",
    "description": "Lipophilicity",
    "models": [
    {
      "skip": true,
      "specific": true,
      "name": "logP"
    },
    {
      "skip": true,
      "specific": true,
      "name": "logD"
    }]
  },
  "Solubility": {
    "name": "Solubility",
    "description": "Solubility",
    "models": [
    {
      "skip": true,
      "specific": true,
      "name": "logS"
    }]
  }
}

export const models: any = {
  "CL-Micro": "Estimates intrinsic clearance which is related ro compound's metabolic stability. This characterestic is specific to analyzed problem.",
  "Half-Life": "This estimates metabolic half-time of the drug or the time when concentration drops by half. This characterestic is specific to analyzed problem."
}


export const template = {
  "name": "Demo template",
  "subgroup": [{
    "name": "Absorption",
		"models": [
			{
				"name": "Pgp-Substrate",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "A probability of being a substrate of P-glycoprotein which is responsible for cell membrane permeability. Compounds with high molecular mass and a large number of polar atoms are the most probable substrates. Binding the substrate leads to low cell permeability of substance.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "-",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Lower is better",
          "enable": false, 
        }],
        "coloring": {
          "type": "Linear",
          "min": "0",
          "max": "1",
          "colors": `[${DG.Color.red}, ${DG.Color.green}]`,
        },
        "preference": "lower",
			},
      {
				"name": "Caco2",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "It is a rate of compound flux across polarized human colon carcinoma Caco-2 cells layers, cm/s (VALIDATE). High values indicate good permeability characteristics.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "cm/s",
          "enable": false,
        }],
				"coloring": {
          "type": "Linear",
          "min": "-6",
          "max": "-5.5",
          "colors": `[${DG.Color.red}, ${DG.Color.green}]`,
        },
        "preference": "higher",
			},
      {
				"name": "Lipophilicity",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "This is a charachteristic of compound quatifying the distribution between lipophilic and hydrophiliv environments. High lipophilicity values  are good for the most of memranes permeability, though low values are related to achievment of therapeutical plasma concentrations. This characterestic is specific to analyzed problem.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "log-ratio",
          "enable": false,
        }],
				"coloring": {
          "type": "Linear",
          "min": "0",
          "max": "5",
          "colors": `[${DG.Color.red}, ${DG.Color.green}]`,
        },
        "preference": "",
			},
      {
				"name": "Solubility",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "Indicates the fraction of compound that can be solved in water solution, e.g. plasma. This characterestic is specific to analyzed problem.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "log mol/L",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": false,
        }],
				"coloring": {
          "type": "Linear",
          "min": "-4.5",
          "max": "-1",
          "colors": `[${DG.Color.red}, ${DG.Color.green}]`,
        },
        "preference": "higher",
			},
    ]
  },
  {
    "name": "Distribution",
    "models": [
      {
        "name": "PPBR",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "This value stands for plasma protein binding. PPB strongly influences drug behavior as only free drug is able to penetrate mebranes. Low values indicate good permeability characteristics.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "%",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      },
      {
        "name": "VDss",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "Volume of distribution (L) is a pharmacokinetic parameter representing the ability of drug to remain in plasma or redistribute to other tissues. Hih VD has a propensity to leave the plasma to other tissues, low VD indicates the ability of drug to achieve desired concentration at lower doses.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "L/kg",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      }
    ]
  },
  {
    "name": "Metabolism",
    "models": [
      {
        "name": "CYP1A2-Inhibitor",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "A probability of being an inhibitor of cytochrome CYP1A2 which catalyzes drug metabolism. Inhibiting CYP1A2 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "-",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      },
      {
        "name": "CYP2C19-Inhibitor",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "A probability of being an inhibitor of cytochrome CYP2C19 which catalyzes drug metabolism. Inhibiting CYP2C19 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "-",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      },
      {
        "name": "CYP2C9-Inhibitor",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "A probability of being an inhibitor of cytochrome CYP2C9 which catalyzes drug metabolism. Inhibiting CYP2C9 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "-",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      },
      {
        "name": "CYP2C9-Substrate",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "A probability of being an inhibitor of substrate CYP2C9 which catalyzes drug metabolism. Binding CYP2C9 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "-",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true, 
        }]
      },
      {
        "name": "CYP2D6-Inhibitor",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "A probability of being an inhibitor of cytochrome CYP2D6 which catalyzes drug metabolism. Inhibiting CYP2D6 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "-",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      },
      {
        "name": "CYP2D6-Substrate",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "A probability of being an inhibitor of substrate CYP2D6 which catalyzes drug metabolism. Binding CYP2D6 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "-",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      }
    ]
  },
  {
    "name": "Excretion",
    "models": [
      {
        "name": "CL-Hepa",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "Estimate for a hepatocyte clearance is excretion characteristic which contributes to  pojection of dose and drug exposure. This characterestic is specific to analyzed problem.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "mL/min",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      },
      {
        "name": "CL-Micro",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "Estimates intrinsic clearance which is related ro compound's metabolic stability. This characterestic is specific to analyzed problem.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "μL/min/million",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      },
      {
        "name": "Half-Life",
        "properties": [{
          "name": "description",
          "inputType": "TextArea",
          "defaultValue": "This estimates metabolic half-time of the drug or the time when concentration drops by half. This characterestic is specific to analyzed problem.",
          "enable": false,
        },
        {
          "name": "units",
          "inputType": "Text",
          "defaultValue": "hours",
          "enable": false,
        },
        {
          "name": "radio",
          "inputType": "Radio",
          "choices": ["Lower is better", "Higher is better"],
          "defaultValue": "Higher is better",
          "enable": true,
        }]
      }
    ]
  }
  ]
}