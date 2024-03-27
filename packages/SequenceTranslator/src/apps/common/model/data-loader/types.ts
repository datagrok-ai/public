type KeyToValue = {[key: string]: string};

export type Edges = {
  [key: string]: KeyToValue
}

export type AxolabsStyle = {
  [index: string]: {
    fullName: string,
    symbols: string[],
    color: string,
  }
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
