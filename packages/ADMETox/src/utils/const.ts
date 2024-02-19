export const properties: any = {
  "Absorption": {
    "name": "Absorption",
    "description": "Models to predict absorbtion",
    "models": [
    {
      "skip": true,
      "specific": false,
      "name": "Pgp-Inhibitor"
    },
		{
      "skip": false,
      "specific": false,
      "name": "Pgp-Substrate"
    },
    {
      "skip": true,
      "specific": true,
      "name": "HIA"
    },
    {
      "skip": true,
      "specific": true,
      "name": "F(20%)"
    },
    {
      "skip": true,
      "specific": true,
      "name": "F(30%)"
    },
    {
      "skip": false,
      "specific": false,
      "name": "Caco2",
    },
    {
      "skip": false,
      "specific": false,
      "name": "Lipophilicity",
    },
    {
      "skip": false,
      "specific": true,
      "name": "Solubility",
    }]
  },
  "Distribution": {
    "name": "Distribution",
    "description": "Models to predict distribution",
    "models": [
    {
      "skip": false,
      "specific": false,
      "name": "PPBR"
    },
    {
      "skip": false,
      "specific": false,
      "name": "VDss"
    },
    {
      "skip": true,
      "specific": true,
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
      "name": "CYP1A2-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "name": "CYP1A2-Substrate"
    },
    {
      "skip": true,
      "specific": true,
      "name": "CYP3A4-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "name": "CYP3A4-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "name": "CYP2C19-Inhibitor"
    },
    {
      "skip": true,
      "specific": true,
      "name": "CYP2C19-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "name": "CYP2C9-Inhibitor"
    },
    {
      "skip": false,
      "specific": true,
      "name": "CYP2C9-Substrate"
    },
    {
      "skip": false,
      "specific": true,
      "name": "CYP2D6-Inhibitor"
    },
    {
      "skip": false,
      "specific": true,
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
      "name": "CL-Hepa"
    },
    {
      "skip": false,
      "specific": true,
      "name": "CL-Micro"
    },
    {
      "skip": false,
      "specific": true,
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
