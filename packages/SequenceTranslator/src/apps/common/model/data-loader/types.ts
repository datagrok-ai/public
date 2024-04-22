import {MODIFICATION_DATA_FIELDS as F} from './const';

type KeyToValue = {[key: string]: string};

export type Edges = {
  [key: string]: KeyToValue
}

export type ModificationData = {
  [F.COLOR]: string,
  [F.SUBSTITUTION]: string,
}

export type ModificationEntry = {
  [modification: string]: ModificationData
}

export type PatternAppData = {
  [format: string]: ModificationEntry
};

export type CodesInfo = {
    [key: string]: { // nucleoside or phosphate
      [code: string]: string
    }
  }


export type FormatToHELMDict = {
  [sourceFormat: string]: CodesInfo
}

export type CodeToSymbol = {
  [format: string]: {
    [code: string]: string
  }
}
