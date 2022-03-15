import {map, stadardPhosphateLinkSmiles, SYNTHESIZERS, TECHNOLOGIES, MODIFICATIONS} from './map';
import {isValidSequence} from './sequence-codes-tools';
import {getMol} from './mol-transformations';

export function sequenceToMolV3000(sequence: string, inverted: boolean = false): string {
  const obj = getObjectWithCodesAndSmiles(sequence);
  let codes = sortByStringLengthInDescendingOrder(Object.keys(obj));
  let i = 0;
  const smilesCodes:string[] = [];
  const codesList = [];
  const links = ['s', 'ps', '*'];
  const includesStandardLinkAlready = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns);
  while (i < sequence.length) {
    const code = codes.find((s: string) => s == sequence.slice(i, i + s.length))!;
    i += code.length;
    inverted ? codesList.unshift(code) : codesList.push(code);
  }
  for (let i = 0; i < codesList.length; i++) {
    if (dropdowns.includes(codesList[i])) {
      if (i == codesList.length -1 || (i < codesList.length - 1 && links.includes(codesList[i + 1]))) {
        smilesCodes.push((i >= codesList.length / 2) ?
          MODIFICATIONS[codesList[i]].right:
          MODIFICATIONS[codesList[i]].left);
      } else if (i < codesList.length - 1) {
        smilesCodes.push((i >= codesList.length / 2) ?
          MODIFICATIONS[codesList[i]].right : MODIFICATIONS[codesList[i]].left);
        smilesCodes.push(stadardPhosphateLinkSmiles);
      }
    } else {
      if (links.includes(codesList[i]) ||
        includesStandardLinkAlready.includes(codesList[i]) ||
        (i < codesList.length - 1 && links.includes(codesList[i + 1]))
      )
        smilesCodes.push(obj[codesList[i]]);
      else {
        smilesCodes.push(obj[codesList[i]]);
        smilesCodes.push(stadardPhosphateLinkSmiles);
      }
    }
  }

  return getMol(smilesCodes);
}

export function sequenceToSmiles(sequence: string, inverted: boolean = false): string {
  const obj = getObjectWithCodesAndSmiles(sequence);
  let codes = sortByStringLengthInDescendingOrder(Object.keys(obj));
  let i = 0;
  let smiles = '';
  const codesList = [];
  const links = ['s', 'ps', '*'];
  const includesStandardLinkAlready = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns);
  while (i < sequence.length) {
    const code = codes.find((s: string) => s == sequence.slice(i, i + s.length))!;
    i += code.length;
    inverted ? codesList.unshift(code) : codesList.push(code);
  }
  for (let i = 0; i < codesList.length; i++) {
    if (dropdowns.includes(codesList[i])) {
      if (i == codesList.length -1 || (i < codesList.length - 1 && links.includes(codesList[i + 1]))) {
        smiles += (i >= codesList.length / 2) ?
          MODIFICATIONS[codesList[i]].right:
          MODIFICATIONS[codesList[i]].left;
      } else if (i < codesList.length - 1) {
        smiles += (i >= codesList.length / 2) ?
          MODIFICATIONS[codesList[i]].right + stadardPhosphateLinkSmiles:
          MODIFICATIONS[codesList[i]].left + stadardPhosphateLinkSmiles;
      }
    } else {
      if (links.includes(codesList[i]) ||
        includesStandardLinkAlready.includes(codesList[i]) ||
        (i < codesList.length - 1 && links.includes(codesList[i + 1]))
      )
        smiles += obj[codesList[i]];
      else
        smiles += obj[codesList[i]] + stadardPhosphateLinkSmiles;
    }
  }
  smiles = smiles.replace(/OO/g, 'O');
  return (
    (
      links.includes(codesList[codesList.length - 1]) &&
      codesList.length > 1 &&
      !includesStandardLinkAlready.includes(codesList[codesList.length - 2])
    ) ||
    dropdowns.includes(codesList[codesList.length - 1]) ||
    includesStandardLinkAlready.includes(codesList[codesList.length - 1])
  ) ?
    smiles :
    smiles.slice(0, smiles.length - stadardPhosphateLinkSmiles.length + 1);
}

function getObjectWithCodesAndSmiles(sequence: string) {
  const obj: { [code: string]: string } = {};
  for (const synthesizer of Object.keys(map)) {
    for (const technology of Object.keys(map[synthesizer])) {
      for (const code of Object.keys(map[synthesizer][technology]))
        obj[code] = map[synthesizer][technology][code].SMILES;
    }
  }
  // TODO: create object based from synthesizer type to avoid key(codes) duplicates
  const output = isValidSequence(sequence);
  if (output.expectedSynthesizer == SYNTHESIZERS.MERMADE_12)
    obj['g'] = map[SYNTHESIZERS.MERMADE_12][TECHNOLOGIES.SI_RNA]['g'].SMILES;
  else if (output.expectedSynthesizer == SYNTHESIZERS.AXOLABS)
    obj['g'] = map[SYNTHESIZERS.AXOLABS][TECHNOLOGIES.SI_RNA]['g'].SMILES;
  return obj;
}

function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a: string, b: string) {return b.length - a.length;});
}
