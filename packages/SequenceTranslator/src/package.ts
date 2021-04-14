/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

//name: classicToBioSpring
//input: string nucleotides {semType: nucleotides}
//output: string result {semType: BioSpring / Gapmers}
export function classicToBioSpring(nucleotides: string) {
  let count: number = -1;
  const objForEdges: {[index: string]: string} = {"T": "5*", "A": "6*", "C": "7*", "G": "8*"};
  const objForCenter: {[index: string]: string} = {"C": "9*", "A": "A*", "T": "T*", "G": "G*"};
  return nucleotides.replace(/[ATCG]/g, function (x: string) {
    count++;
    if (count < 5) return objForEdges[x];
    if (count < 15) return objForCenter[x];
    return objForEdges[x];
  }).slice(0, 2 * count + 1);
}

//name: classicToGCRS
//input: string nucleotides {semType: nucleotides}
//output: string result {semType: GCRS / Gapmers}
export function classicToGCRS(nucleotides: string) {
  let count: number = -1;
  const objForEdges: {[index: string]: string} = {"T": "moeUnps", "A": "moeAnps", "C": "moe5mCnps", "G": "moeGnps"};
  const objForCenter: {[index: string]: string} = {"C": "Cps", "A": "Aps", "T": "Tps", "G": "Gps"};
  return nucleotides.replace(/[ATCG]/g, function (x: string) {
    count++;
    if (count < 5) return objForEdges[x];
    if (count < 15) return objForCenter[x];
    return objForEdges[x];
  }).slice(0, -3);
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
  const obj: {[index: string]: string} = {
    "5*": "moeUnps", "6*": "moeAnps", "7*": "moe5mCnps", "8*": "moeGnps", "9*": "5mCps", "A*": "Aps", "T*": "Tps",
    "G*": "Gps", "C*": "Cps", "5": "moeU", "6": "moeA", "7": "moe5mC", "8": "moeG"
  };
  return nucleotides.replace(/(5\*|6\*|7\*|8\*|9\*|A\*|T\*|G\*|C\*|5|6|7|8])/g, function (x: string) {return obj[x];});
}

//name: GCRSToBioSpring
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function GCRSToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "moeT": "5", "moeA": "6", "moe5mC": "7", "moeG": "8", "moeU": "5", "5mC": "9", "nps": "*", "ps": "*", "U": "T"
  };
  return nucleotides.replace(/(moeT|moeA|moe5mC|moeG|moeU|5mC|nps|ps|U)/g, function (x: string) {return obj[x];});
}

//name: GCRSToClassic
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: nucleotides}
export function GCRSToClassic(nucleotides: string) {
  const obj: {[index: string]: string} = {"moe": "", "5m": "", "n": "", "ps": "", "U": "T"};
  return nucleotides.replace(/(moe|5m|n|ps|U)/g, function (x: string) {return obj[x];});
}

//name: GCRSToAxolabs
//input: string nucleotides {semType: GCRS / siRNA}
//output: string result {semType: Axolabs / siRNA}
export function GCRSToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "fU": "Uf", "fA": "Af", "fC": "Cf", "fG": "Gf", "mU": "u", "mA": "a", "mC": "c", "mG": "g", "ps": "s"
  };
  return nucleotides.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: axolabsToGCRS
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: GCRS / siRNA}
export function axolabsToGCRS(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "fU", "Af": "fA", "Cf": "fC", "Gf": "fG", "u": "mU", "a": "mA", "c": "mC", "g": "mG", "s": "ps"
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: axolabsToBioSpring
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: BioSpring / siRNA}
export function axolabsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "1", "Af": "2", "Cf": "3", "Gf": "4", "u": "5", "a": "6", "c": "7", "g": "8", "s": "*"
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: axolabsToClassic
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: nucleotides}
export function axolabsToClassic(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "U", "Af": "A", "Cf": "C", "Gf": "G", "u": "U", "a": "A", "c": "C", "g": "G", "s": ""
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
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