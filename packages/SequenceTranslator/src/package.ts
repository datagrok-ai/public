/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

function r2C(nucleotide: string) {
  switch (nucleotide) {
    case "5": return "T";
    case "6": return "A";
    case "7": return "C";
    case "8": return "G";
    case "9": return "C";
    default: return nucleotide;
  }
}

//name: robotToClassic
//input: string nucleotides {semType: dna_sequence/robot}
//output: string result {semType: dna_sequence/classic}
export function robotToClassic(nucleotides: string) {
  return nucleotides
      .split("*")
      .map(r2C)
      .join("");
}
