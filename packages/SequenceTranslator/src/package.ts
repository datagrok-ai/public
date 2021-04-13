/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

//name: classicToBioSpring {semType: nucleotides}
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result
export function classicToBioSpring(nucleotides: string) {
  return nucleotides;
}

//name: classicToGCRS
//input: string nucleotides {semType: nucleotides}
//output: string result {semType: GCRS / Gapmers}
export function classicToGCRS(nucleotides: string) {
  return nucleotides;
}

//name: bioSpringToClassic
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: nucleotides}
export function bioSpringToClassic(nucleotides: string) {
  const obj: {[index: string]: string} = {"*": "", "5": "T", "6": "A", "7": "C", "8": "G", "9": "C"};
  return nucleotides.replace(/[*56789]/g, function (x: string) {return obj[x];});
}

//name: bioSpringToGCRS
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: GCRS / Gapmers}
export function bioSpringToGCRS(nucleotides: string) {
  const obj: {[index: string]: string} = {"*": "nps", "5": "moeU", "6": "moeA", "7": "moe5mC", "8": "moeG", "9": "5mC"};
  return nucleotides.replace(/[*56789]/g, function (x: string) {return obj[x];});
}

//name: GCRSToBioSpring
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function GCRSToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "moeT": "5", "moeA": "6", "moe5mC": "7", "moeG": "8", "5mC": "9", "nps": "*", "ps": "*", "U": "T"
  };
  return nucleotides.replace(/(moeT|moeA|moe5mC|moeG|5mC|nps|ps|U)/g, function (x: string) {return obj[x];});
}

//name: GCRSToClassic
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: nucleotides}
export function GCRSToClassic(nucleotides: string) {
  const obj: {[index: string]: string} = {"moe": "", "5m": "", "n": "", "ps": "", "U": "T"};
  return nucleotides.replace(/(moe|5m|n|ps|U)/g, function (x: string) {return obj[x];});
}

//name: classicComplement
//input: string nucleotides {semType: nucleotides}
//output: string result {semType: nucleotides}
export function classicComplement(nucleotides: string) {
  const obj: {[index: string]: string} = {"A": "T", "T": "A", "G": "C", "C": "G"};
  return nucleotides.replace(/[ATGC]/g, function (x: string) {return obj[x];});
}

//name: bioSpringComplement
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function bioSpringComplement(nucleotides: string) {
  const obj: {[index: string]: string} = {"A": "T", "T": "A", "G": "C", "C": "G", "5": "6", "6": "5", "7": "8", "8": "9", "9": "8"};
  return nucleotides.replace(/[56789ATGC]/g, function (x: string) {return obj[x];});
}

//name: GCRSComplement
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: GCRS / Gapmers}
export function GCRSComplement(nucleotides: string) {
  const obj: {[index: string]: string} = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"};
  return nucleotides.replace(/[ACGTU]/g, function (x: string) {return obj[x];});
}