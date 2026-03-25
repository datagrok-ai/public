export interface FcsHeader {
  version: string;
  textStart: number;
  textEnd: number;
  dataStart: number;
  dataEnd: number;
  analysisStart: number;
  analysisEnd: number;
}

export interface FcsParameter {
  index: number;
  shortName: string;
  longName: string;
  bits: number;
  range: number;
  amplification: string;
}

export interface FcsSpillover {
  n: number;
  paramNames: string[];
  matrix: number[][];
}

export interface FcsParsed {
  header: FcsHeader;
  keywords: Map<string, string>;
  parameters: FcsParameter[];
  eventCount: number;
  paramCount: number;
  dataType: 'F' | 'D' | 'I';
  littleEndian: boolean;
  data: Float32Array[];
  spillover: FcsSpillover | null;
}
