import {IChemPropertyGroupMap} from './types';

export const ChemPropertyGroupMap: IChemPropertyGroupMap = {
  toxRisks: {
    name: 'Toxicity Risks',
    values: [
      {name: 'Mutagenicity', description: 'Mutagenicity', propertyName: 'mutagenicity'},
      {name: 'Tumorigenicity', description: 'Tumorigenicity', propertyName: 'tumorigenicity'},
      {name: 'Irritating effects', description: 'Irritating effects', propertyName: 'irritatingEffects'},
      {name: 'Reproductive effects', description: 'Reproductive effects', propertyName: 'reproductiveEffects'},
    ],
  },
  structuralAlerts: {
    name: 'Structural Alerts',
    values: [
      {name: 'PAINS', description: 'Pan Assay Interference Compounds filters', propertyName: 'pains'},
      {name: 'BMS', description: 'Bristol-Myers Squibb HTS Deck filters', propertyName: 'bms'},
      {name: 'SureChEMBL', description: 'MedChem unfriendly compounds from SureChEMBL', propertyName: 'sureChembl'},
      {name: 'MLSMR', description: 'NIH MLSMR Excluded Functionality filters', propertyName: 'mlsmr'},
      {name: 'Dandee', description: 'University of Dundee NTD Screening Library filters', propertyName: 'dandee'},
      {name: 'Inpharmatica', description: 'Inpharmatica filters', propertyName: 'inpharmatica'},
      {name: 'LINT', description: 'Pfizer LINT filters', propertyName: 'lint'},
      {name: 'Glaxo', description: 'Glaxo Wellcome Hard filters', propertyName: 'glaxo'},
    ],
  },
  chemProperties: {
    name: 'Chemical Properties',
    values: [
      {name: 'Molecular weight', description: 'Molecular weight', propertyName: 'MW'},
      {name: 'HBA', description: 'HBA', propertyName: 'HBA'},
      {name: 'HBD', description: 'HBD', propertyName: 'HBD'},
      {name: 'LogP', description: 'LogP', propertyName: 'logP'},
      {name: 'LogS', description: 'logS', propertyName: 'logS'},
      {name: 'PSA', description: 'PSA', propertyName: 'PSA'},
      {name: 'Rotatable bonds', description: 'Rotatable bonds', propertyName: 'rotatableBonds'},
      {name: 'Stereo centers', description: 'Stereo centers', propertyName: 'stereoCenters'},
      {name: 'Molecule charge', description: 'Molecule charge', propertyName: 'moleculeCharge'},
    ],
  },
};
