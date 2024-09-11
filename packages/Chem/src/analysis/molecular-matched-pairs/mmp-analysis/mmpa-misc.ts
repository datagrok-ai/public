import BitArray from '@datagrok-libraries/utils/src/bit-array';

export enum MMP_CONSTRICTIONS {
    CPU = 1E4,
    GPU = 1E5
  }

export enum MMP_ERRORS {
    FRAGMENTS_CPU = 'No GPU found and 10,000 molecules is upper limit for MMPa with CPU',
    FRAGMENTS_GPU = 'Upper limit for MMPa with GPU is 100,000 molecules',
    GPU_ABORTED = 'GPU calculations were aborted - faling back to CPU',
    PAIRS = 'Unable to calculate pairs in MMPa analysis',
    GENERATIONS = 'Unable to calculate generations in MMPa analysis'
}

export type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

export type MolecularPair = {
  first: number,
  second: number,
  core: string,
  firstR: string,
  secondR: string
};

export type MmpInitData = {
    molName: string;
    molecules: string[];
    activities: Float32Array [];
    activitiesNames: string[];
    activitiesCount: number;
};

export type MmpRulesBasedData = {
    fromFrag: string[];
    toFrag: string[];
    occasions: Int32Array;
    meanDiffs: Float32Array[];
  };

export type MmpAllCasesBasedData = {
    maxActs: number [];
    molFrom: string[];
    molTo: string[];
    pairNum: Int32Array;
    molNumFrom: Int32Array;
    molNumTo: Int32Array;
    pairsFromSmiles: string[];
    pairsToSmiles: string[];
    ruleNum: Int32Array;
    diffs: Float32Array[];
    activityPairsIdxs: BitArray[];
};

export type MmpGeneration = {
  allStructures: string [];
  allInitActivities: Float32Array;
  cores: string [];
  from: string [];
  to: string [];
  prediction: Float32Array;
  activityName: string [];
};
