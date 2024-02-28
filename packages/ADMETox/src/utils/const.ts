import * as DG from 'datagrok-api/dg';

export const properties: any = {
  "Absorption": {
    "name": "Absorption",
    "description": "Models to predict absorbtion",
    "models": [
    {
      "skip": true,
      "specific": false,
      "units": "-",
      "interpretation": "",
      "coloring": "",
      "name": "Pgp-Inhibitor"
    },
		{
      "skip": false,
      "specific": false,
      "units": "-",
      "interpretation": "Lower is better",
      "coloring": {
        "type": "Linear",
        "min": "0",
        "max": "1",
        "colors": `[${DG.Color.red}, ${DG.Color.green}]`,
      },
      "name": "Pgp-Substrate"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "",
      "coloring": "",
      "name": "HIA"
    },
    {
      "skip": false,
      "specific": false,
      "units": "cm/s",
      "interpretation": {
        "range": {
          "below_-6": "Bad",
          "to_-6_-5.5": "Intermediate",
          "above_-5.5": "Good",
        }
      },
      "coloring": {
        "type": "Linear",
        "min": "-6",
        "max": "-5.5",
        "colors": `[${DG.Color.red}, ${DG.Color.green}]`,
      },
      "name": "Caco2"
    },
    {
      "skip": false,
      "specific": false,
      "units": "log-ratio",
      "interpretation": {
        "range": {
          "below_0": "Bad",
          "to_2_3": "Ideal",
          "above_5": "Bad",
        }
      },
      "coloring": "",
      "name": "Lipophilicity"
    },
    {
      "skip": false,
      "specific": true,
      "units": "log mol/L",
      "interpretation": "The higher the better",
      "coloring": {
        "type": "Linear",
        "min": "-4.5",
        "max": "-1",
        "colors": `[${DG.Color.red}, ${DG.Color.green}]`,
      },
      "name": "Solubility"
    }]
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
      "name": "VDss"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "",
      "coloring": "",
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
      "name": "CYP1A2-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP1A2-Substrate"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP3A4-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP3A4-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP2C19-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP2C19-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP2C9-Inhibitor"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP2C9-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
      "name": "CYP2D6-Inhibitor"
    },
    {
      "skip": false,
      "specific": true,
      "units": "-",
      "interpretation": "No good or bad values",
      "coloring": "",
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
      "name": "CL-Micro"
    },
    {
      "skip": false,
      "specific": true,
      "units": "hours",
      "interpretation": "No good or bad values",
      "coloring": "",
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
  "Pgp-Inhibitor": "A probability of being an inhibitor of P-glycoprotein which is responsible for cell membrane transport. Inhibiting Pgp leads to low cell permeability of substance.",
  "Pgp-Substrate": "A probability of being a substrate of P-glycoprotein which is responsible for cell membrane permeability. Compounds with high molecular mass and a large number of polar atoms are the most probable substrates. Binding the substrate leads to low cell permeability of substance.",
  "Caco2": "It is a rate of compound flux across polarized human colon carcinoma Caco-2 cells layers, cm/s (VALIDATE). High values indicate good permeability characteristics.",
  "Lipophilicity": "This is a charachteristic of compound quatifying the distribution between lipophilic and hydrophiliv environments. High lipophilicity values  are good for the most of memranes permeability, though low values are related to achievment of therapeutical plasma concentrations. This characterestic is specific to analyzed problem.",
  "Solubility": "Indicates the fraction of compound that can be solved in water solution, e.g. plasma. This characterestic is specific to analyzed problem.",
  "PPBR": "This value stands for plasma protein binding. PPB strongly influences drug behavior as only free drug is able to penetrate mebranes. Low values indicate good permeability characteristics.",
  "VDss": "Volume of distribution (L) is a pharmacokinetic parameter representing the ability of drug to remain in plasma or redistribute to other tissues. Hih VD has a propensity to leave the plasma to other tissues, low VD indicates the ability of drug to achieve desired concentration at lower doses.",
  "CYP1A2-Inhibitor": "A probability of being an inhibitor of cytochrome CYP1A2 which catalyzes drug metabolism. Inhibiting CYP1A2 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
  "CYP3A4-Inhibitor": "A probability of being an inhibitor of cytochrome CYP3A4 which catalyzes drug metabolism. Inhibiting CYP3A4 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
  "CYP3A4-Substrate": "A probability of being an inhibitor of substrate CYP3A4 which catalyzes drug metabolism. Binding CYP3A4 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.", 
  "CYP2C19-Inhibitor": "A probability of being an inhibitor of cytochrome CYP2C19 which catalyzes drug metabolism. Inhibiting CYP2C19 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
  "CYP2C9-Inhibitor": "A probability of being an inhibitor of cytochrome CYP2C9 which catalyzes drug metabolism. Inhibiting CYP2C9 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
  "CYP2C9-Substrate": "A probability of being an inhibitor of substrate CYP2C9 which catalyzes drug metabolism. Binding CYP2C9 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
  "CYP2D6-Inhibitor": "A probability of being an inhibitor of cytochrome CYP2D6 which catalyzes drug metabolism. Inhibiting CYP2D6 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
  "CYP2D6-Substrate": " A probability of being an inhibitor of substrate CYP2D6 which catalyzes drug metabolism. Binding CYP2D6 could be both a desired and undesired property depending on drug development goals: elongation/reduction of drug effect, drug-drug interactions etc.",
  "CL-Hepa": "Estimate for a hepatocyte clearance is excretion characteristic which contributes to  pojection of dose and drug exposure. This characterestic is specific to analyzed problem.",
  "CL-Micro": "Estimates intrinsic clearance which is related ro compound's metabolic stability. This characterestic is specific to analyzed problem.",
  "Half-Life": "This estimates metabolic half-time of the drug or the time when concentration drops by half. This characterestic is specific to analyzed problem."
}
