import {lcmsToGcrs, MODIFICATIONS} from './map';
import * as DG from 'datagrok-api/dg';
import {DELIMITER} from './map';
import {sortByStringLengthInDescendingOrder} from '../utils/helpers';
//name: gcrsToLcms
//input: string nucleotides {semType: GCRS}
//output: string result {semType: LCMS}
export function gcrsToLcms(sequence: string): string {
  try {
    const df = DG.DataFrame.fromCsv(lcmsToGcrs);
    const arr1: string[] = df.getCol('GCRS').toList();
    const arr2: string[] = df.getCol('LCMS').toList();
    const obj: { [i: string]: string } = {};
    arr1.forEach((element, index) => obj[element] = arr2[index]);
    obj[DELIMITER] = DELIMITER;
    const codes = arr1
      .concat(DELIMITER)
      .concat(Object.keys(MODIFICATIONS));
    const sortedCodes = sortByStringLengthInDescendingOrder(codes);
    let i = 0;
    let r1 = '';
    while (i < sequence.length) {
      const matchedCode = sortedCodes.find((c) => c == sequence.slice(i, i + c.length))!;
      r1 += obj[sequence.slice(i, i + matchedCode.length)];
      i += matchedCode.length;
    }
    while (r1.indexOf('//') != -1)
      r1 = r1.replace('//', '/');
    return r1;
  } catch {
    return '<error>';
  }
}

//name: asoGapmersNucleotidesToBioSpring
//input: string nucleotides {semType: DNA nucleotides}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersNucleotidesToBioSpring(nucleotides: string): string {
  let count: number = -1;
  const objForEdges: { [index: string]: string } = {
    '(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'T': '5*', 'A': '6*', 'C': '7*', 'G': '8*'};
  const objForCenter: { [index: string]: string } = {
    '(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'T': 'T*', 'A': 'A*', 'C': '9*', 'G': 'G*'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|T|C|G)/g, function(x: string) {
    count++;
    return (count > 4 && count < 15) ? objForCenter[x] : objForEdges[x];
  }).slice(0, (nucleotides.endsWith('(invabasic)') || nucleotides.endsWith('(GalNAc-2-JNJ)')) ?
    nucleotides.length : 2 * count + 1);
}

//name: asoGapmersNucleotidesToGcrs
//input: string nucleotides {semType: DNA nucleotides}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersNucleotidesToGcrs(nucleotides: string): string {
  let count: number = -1;
  const objForEdges: { [index: string]: string } = {
    '(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'T': 'moeUnps',
    'A': 'moeAnps', 'C': 'moe5mCnps', 'G': 'moeGnps'};
  const objForCenter: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'C': '5mCps', 'A': 'Aps', 'T': 'Tps', 'G': 'Gps'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|T|C|G)/g, function(x: string) {
    count++;
    if (count < 5) return (count == 4) ? objForEdges[x].slice(0, -3) + 'ps' : objForEdges[x];
    if (count < 15) return (count == 14) ? objForCenter[x].slice(0, -2) + 'nps' : objForCenter[x];
    return objForEdges[x];
  }).slice(0, (nucleotides.endsWith('(invabasic)') || nucleotides.endsWith('(GalNAc-2-JNJ)')) ?
    nucleotides.length : -3);
}

//name: asoGapmersBioSpringToNucleotides
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: DNA nucleotides}
export function asoGapmersBioSpringToNucleotides(nucleotides: string): string {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    '*': '', '5': 'T', '6': 'A', '7': 'C', '8': 'G', '9': 'C'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|5|6|7|8|9)/g, function(x: string) {return obj[x];});
}

//name: asoGapmersBioSpringToGcrs
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersBioSpringToGcrs(nucleotides: string): string {
  let count: number = -1;
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    '5*': 'moeUnps', '6*': 'moeAnps', '7*': 'moe5mCnps', '8*': 'moeGnps', '9*': '5mCps', 'A*': 'Aps', 'T*': 'Tps',
    'G*': 'Gps', 'C*': 'Cps', '5': 'moeU', '6': 'moeA', '7': 'moe5mC', '8': 'moeG',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|5\*|6\*|7\*|8\*|9\*|A\*|T\*|G\*|C\*|5|6|7|8)/g,
    function(x: string) {
      count++;
      return (count == 4) ? obj[x].slice(0, -3) + 'ps' : (count == 14) ? obj[x].slice(0, -2) + 'nps' : obj[x];
    });
}


//name: asoGapmersGcrsToBioSpring
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersGcrsToBioSpring(nucleotides: string): string {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'fU': '1', 'fA': '2', 'fC': '3', 'fG': '4', 'mU': '5', 'mA': '6', 'mC': '7', 'mG': '8',
    'moeT': '5', 'moeA': '6', 'moe5mC': '7', 'moeG': '8', 'moeU': '5', '5mC': '9', 'nps': '*', 'ps': '*', 'U': 'T',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|fU|fA|fC|fG|mU|mA|mC|mG|moeT|moeA|moe5mC|moeG|moeU|5mC|nps|ps|U)/g,
    function(x: string) {return obj[x];});
}

//name: asoGapmersGcrsToNucleotides
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: DNA nucleotides}
export function asoGapmersGcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'moe': '', '5m': '', 'n': '', 'ps': '', 'U': 'T'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|moe|5m|n|ps|U)/g, function(x: string) {return obj[x];});
}

//name: siRnaBioSpringToNucleotides
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: RNA nucleotides}
export function siRnaBioSpringToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    '1': 'U', '2': 'A', '3': 'C', '4': 'G', '5': 'U', '6': 'A', '7': 'C', '8': 'G', '*': ''};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|1|2|3|4|5|6|7|8|\*)/g,
    function(x: string) {return obj[x];});
}

//name: siRnaBioSpringToAxolabs
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: Axolabs / siRNA}
export function siRnaBioSpringToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    '1': 'Uf', '2': 'Af', '3': 'Cf', '4': 'Gf', '5': 'u', '6': 'a', '7': 'c', '8': 'g', '*': 's'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|1|2|3|4|5|6|7|8|\*)/g,
    function(x: string) {return obj[x];});
}


//name: siRnaBioSpringToGcrs
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: GCRS}
export function siRnaBioSpringToGcrs(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    '1': 'fU', '2': 'fA', '3': 'fC', '4': 'fG', '5': 'mU', '6': 'mA', '7': 'mC', '8': 'mG', '*': 'ps'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|1|2|3|4|5|6|7|8|\*)/g,
    function(x: string) {return obj[x];});
}

//name: siRnaAxolabsToGcrs
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: GCRS}
export function siRnaAxolabsToGcrs(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'Uf': 'fU', 'Af': 'fA', 'Cf': 'fC', 'Gf': 'fG', 'u': 'mU', 'a': 'mA', 'c': 'mC', 'g': 'mG', 's': 'ps',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|Uf|Af|Cf|Gf|u|a|c|g|s)/g,
    function(x: string) {return obj[x];});
}

//name: siRnaAxolabsToBioSpring
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: BioSpring / siRNA}
export function siRnaAxolabsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'Uf': '1', 'Af': '2', 'Cf': '3', 'Gf': '4', 'u': '5', 'a': '6', 'c': '7', 'g': '8', 's': '*',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|Uf|Af|Cf|Gf|u|a|c|g|s)/g,
    function(x: string) {return obj[x];});
}

//name: siRnaAxolabsToNucleotides
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: RNA nucleotides}
export function siRnaAxolabsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'Uf': 'U', 'Af': 'A', 'Cf': 'C', 'Gf': 'G', 'u': 'U', 'a': 'A', 'c': 'C', 'g': 'G', 's': '',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|Uf|Af|Cf|Gf|u|a|c|g|s)/g,
    function(x: string) {return obj[x];});
}


//name: siRnaGcrsToNucleotides
//input: string nucleotides {semType: GCRS}
//output: string result {semType: RNA nucleotides}
export function siRnaGcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'fU': 'U', 'fA': 'A', 'fC': 'C', 'fG': 'G', 'mU': 'U', 'mA': 'A', 'mC': 'C', 'mG': 'G', 'ps': '',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|fU|fA|fC|fG|mU|mA|mC|mG|ps)/g,
    function(x: string) {return obj[x];});
}

//name: siRnaGcrsToBioSpring
//input: string nucleotides {semType: GCRS}
//output: string result {semType: BioSpring / siRNA}
export function siRnaGcrsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'fU': '1', 'fA': '2', 'fC': '3', 'fG': '4', 'mU': '5', 'mA': '6', 'mC': '7', 'mG': '8', 'ps': '*',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|fU|fA|fC|fG|mU|mA|mC|mG|ps)/g,
    function(x: string) {return obj[x];});
}

//name: siRnaGcrsToAxolabs
//input: string nucleotides {semType: GCRS}
//output: string result {semType: Axolabs / siRNA}
export function siRnaGcrsToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'fU': 'Uf', 'fA': 'Af', 'fC': 'Cf', 'fG': 'Gf', 'mU': 'u', 'mA': 'a', 'mC': 'c', 'mG': 'g', 'ps': 's',
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|fU|fA|fC|fG|mU|mA|mC|mG|ps)/g,
    function(x: string) {return obj[x];});
}

//name: siRnaNucleotideToBioSpringSenseStrand
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: BioSpring / siRNA}
export function siRnaNucleotideToBioSpringSenseStrand(nucleotides: string) {
  let count: number = -1;
  const objForLeftEdge: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'A': '6*', 'U': '5*', 'G': '8*', 'C': '7*'};
  const objForRightEdge: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'A': '*6', 'U': '*5', 'G': '*8', 'C': '*7'};
  const objForOddIndices: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'A': '6', 'U': '5', 'G': '8', 'C': '7'};
  const objForEvenIndices: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': '2', 'U': '1', 'G': '4', 'C': '3'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function(x: string) {
    count++;
    if (count < 2) return objForLeftEdge[x];
    if (count > nucleotides.length - 3) return objForRightEdge[x];
    return (count % 2 == 0) ? objForEvenIndices[x] : objForOddIndices[x];
  });
}

//name: siRnaNucleotidesToGcrs
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: GCRS}
export function siRnaNucleotidesToGcrs(nucleotides: string) {
  let count: number = -1;
  const objForLeftEdge: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'A': 'mAps', 'U': 'mUps', 'G': 'mGps', 'C': 'mCps'};
  const objForRightEdge: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'A': 'psmA', 'U': 'psmU', 'G': 'psmG', 'C': 'psmC'};
  const objForEvenIndices: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': 'fA', 'U': 'fU', 'G': 'fG', 'C': 'fC'};
  const objForOddIndices: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': 'mA', 'U': 'mU', 'G': 'mG', 'C': 'mC'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function(x: string) {
    count++;
    if (count < 2) return objForLeftEdge[x];
    if (count > nucleotides.length - 3) return objForRightEdge[x];
    return (count % 2 == 0) ? objForEvenIndices[x] : objForOddIndices[x];
  });
}

//name: siRnaNucleotideToAxolabsSenseStrand
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function siRnaNucleotideToAxolabsSenseStrand(nucleotides: string) {
  let count: number = -1;
  const objForLeftEdge: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'A': 'as', 'U': 'us', 'G': 'gs', 'C': 'cs'};
  const objForSomeIndices: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': 'Af', 'U': 'Uf', 'G': 'Gf', 'C': 'Cf'};
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'A': 'a', 'U': 'u', 'G': 'g', 'C': 'c'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function(x: string) {
    count++;
    if (count < 2) return objForLeftEdge[x];
    if (count == 6 || (count > 7 && count < 11)) return objForSomeIndices[x];
    if (count == nucleotides.length - 1) return 'a';
    return obj[x];
  });
}

//name: siRnaNucleotideToAxolabsAntisenseStrand
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function siRnaNucleotideToAxolabsAntisenseStrand(nucleotides: string) {
  let count: number = -1;
  const objForSmallLinkages: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': 'as', 'U': 'us', 'G': 'gs', 'C': 'cs'};
  const objForBigLinkages: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': 'Afs', 'U': 'Ufs', 'G': 'Gfs', 'C': 'Cfs'};
  const objForSomeIndices: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': 'Af', 'U': 'Uf', 'G': 'Gf', 'C': 'Cf'};
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)',
    '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)', 'A': 'a', 'U': 'u', 'G': 'g', 'C': 'c'};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function(x: string) {
    count++;
    if (count > 19 && count < 22) return objForSmallLinkages[x];
    if (count == 0) return 'us';
    if (count == 1) return objForBigLinkages[x];
    return (count == 5 || count == 7 || count == 8 || count == 13 || count == 15) ? objForSomeIndices[x] : obj[x];
  });
}

//name: gcrsToNucleotides
//input: string nucleotides {semType: GCRS}
//output: string result {semType: RNA nucleotides}
export function gcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'mAps': 'A', 'mUps': 'U', 'mGps': 'G', 'mCps': 'C', 'fAps': 'A', 'fUps': 'U', 'fGps': 'G', 'fCps': 'C',
    'fU': 'U', 'fA': 'A', 'fC': 'C', 'fG': 'G', 'mU': 'U', 'mA': 'A', 'mC': 'C', 'mG': 'G',
  };
  return nucleotides.replace(
    /(\(invabasic\)|\(GalNAc-2-JNJ\)|mAps|mUps|mGps|mCps|fAps|fUps|fGps|fCps|fU|fA|fC|fG|mU|mA|mC|mG)/g,
    function(x: string) {return obj[x];});
}

//name: gcrsToMermade12
//input: string nucleotides {semType: GCRS}
//output: string result {semType: Mermade 12 / siRNA}
export function gcrsToMermade12(nucleotides: string) {
  const obj: {[index: string]: string} = {'(invabasic)': '(invabasic)', '(GalNAc-2-JNJ)': '(GalNAc-2-JNJ)',
    'mAps': 'e', 'mUps': 'h', 'mGps': 'g', 'mCps': 'f', 'fAps': 'i', 'fUps': 'l', 'fGps': 'k', 'fCps': 'j', 'fU': 'L',
    'fA': 'I', 'fC': 'J', 'fG': 'K', 'mU': 'H', 'mA': 'E', 'mC': 'F', 'mG': 'G',
  };
  return nucleotides.replace(
    /(\(invabasic\)|\(GalNAc-2-JNJ\)|mAps|mUps|mGps|mCps|fAps|fUps|fGps|fCps|fU|fA|fC|fG|mU|mA|mC|mG)/g,
    function(x: string) {return obj[x];});
}
