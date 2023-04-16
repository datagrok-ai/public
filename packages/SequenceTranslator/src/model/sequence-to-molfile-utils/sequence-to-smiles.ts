// import {map, SYNTHESIZERS, TECHNOLOGIES, MODIFICATIONS, DELIMITER} from './map';
import {map} from '../../hardcode-to-be-eliminated/map';
import {SYNTHESIZERS, DELIMITER, TECHNOLOGIES} from '../const';
import {isValidSequence} from '../code-converter/conversion-validation-tools';
import {sortByStringLengthInDescendingOrder} from '../helpers';

import {standardPhosphateLinkSmiles, MODIFICATIONS} from '../../hardcode-to-be-eliminated/const';

const LINKS = ['s', 'ps', '*', 'Rpn', 'Spn', 'Rps', 'Sps'];

export function sequenceToSmiles(sequence: string, inverted: boolean = false, format: string): string {
  const obj = getObjectWithCodesAndSmiles(sequence, format);
  let codes = sortByStringLengthInDescendingOrder(Object.keys(obj));
  let i = 0;
  let smiles = '';
  const codesList = [];
  const links = LINKS;
  const includesStandardLinkAlready = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns).concat(DELIMITER);
  while (i < sequence.length) {
    const code = codes.find((s: string) => s == sequence.slice(i, i + s.length))!;
    i += code.length;
    inverted ? codesList.unshift(code) : codesList.push(code);
  }
  for (let i = 0; i < codesList.length; i++) {
    if (dropdowns.includes(codesList[i])) {
      smiles += (i >= codesList.length / 2) ?
        MODIFICATIONS[codesList[i]].right + standardPhosphateLinkSmiles :
        MODIFICATIONS[codesList[i]].left + standardPhosphateLinkSmiles;
    } else {
      if (links.includes(codesList[i]) ||
        includesStandardLinkAlready.includes(codesList[i]) ||
        (i < codesList.length - 1 && links.includes(codesList[i + 1]))
      )
        smiles += obj[codesList[i]];
      else
        smiles += obj[codesList[i]] + standardPhosphateLinkSmiles;
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
    smiles.slice(0, smiles.length - standardPhosphateLinkSmiles.length + 1);
}


function getObjectWithCodesAndSmiles(sequence: string, format: string) {
  const obj: { [code: string]: string } = {};
  if (format == null) {
    for (const synthesizer of Object.keys(map)) {
      for (const technology of Object.keys(map[synthesizer])) {
        for (const code of Object.keys(map[synthesizer][technology]))
          obj[code] = map[synthesizer][technology][code].SMILES;
      }
    }
  } else {
    for (const technology of Object.keys(map[format])) {
      for (const code of Object.keys(map[format][technology]))
        obj[code] = map[format][technology][code].SMILES;
    }
  }
  obj[DELIMITER] = '';
  // TODO: create object based from synthesizer type to avoid key(codes) duplicates
  const output = isValidSequence(sequence, format);
  if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12))
    obj['g'] = map[SYNTHESIZERS.MERMADE_12][TECHNOLOGIES.SI_RNA]['g'].SMILES;
  else if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS))
    obj['g'] = map[SYNTHESIZERS.AXOLABS][TECHNOLOGIES.SI_RNA]['g'].SMILES;
  return obj;
}
