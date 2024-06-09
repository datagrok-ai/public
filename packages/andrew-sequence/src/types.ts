export type Nucleotide = 'A'| 'T' | 'C' | 'G' | 'a' | 't' | 'c' | 'g';
export type AminoAcid = Nucleotide | 'R' | 'N' | 'D' | 'r' | 'n' | 'd';

export interface EnaSequence {
  seqType: string;
  id: string;
  genBank: string;
  sequence: string;
  code: string;
  description: string;
  name: string;
  extra: string;
  raw: string;
}
