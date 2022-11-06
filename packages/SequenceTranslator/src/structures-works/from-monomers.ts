// import {map, SYNTHESIZERS, TECHNOLOGIES, delimiter} from './map';
import {map, SYNTHESIZERS, TECHNOLOGIES, delimiter} from './map';
import {isValidSequence} from './sequence-codes-tools';
import {getNucleotidesMol} from './mol-transformations';

import {standardPhosphateLinkSmiles, MODIFICATIONS} from './const';

export function sequenceToMolV3000(
  sequence: string, inverted: boolean = false, oclRender: boolean = false,
  format: string, monomersLib: string,
): string {
  const obj = getObjectWithCodesAndSmilesFromFile(sequence, format, monomersLib);
  console.log('obj', obj);
  let codes = sortByStringLengthInDescendingOrder(Object.keys(obj));
  let i = 0;
  const smilesCodes:string[] = [];
  const codesList = [];
  const links = ['s', 'ps', '*'];
  const includesStandardLinkAlready = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns).concat(delimiter);
  while (i < sequence.length) {
    const code = codes.find((s: string) => s === sequence.slice(i, i + s.length))!;
    i += code.length;
    inverted ? codesList.unshift(code) : codesList.push(code);
  }
  for (let i = 0; i < codesList.length; i++) {
    if (dropdowns.includes(codesList[i])) {
      smilesCodes.push((i >= codesList.length / 2) ?
        MODIFICATIONS[codesList[i]].right : MODIFICATIONS[codesList[i]].left);
      if (!(i < codesList.length - 1 && links.includes(codesList[i + 1])))
        smilesCodes.push(standardPhosphateLinkSmiles);
    } else {
      if (links.includes(codesList[i]) ||
        includesStandardLinkAlready.includes(codesList[i]) ||
        (i < codesList.length - 1 && links.includes(codesList[i + 1]))
      )
        smilesCodes.push(obj[codesList[i]]);
      else {
        smilesCodes.push(obj[codesList[i]]);
        smilesCodes.push(standardPhosphateLinkSmiles);
      }
    }
  }

  return getNucleotidesMol(smilesCodes);
}

export function sequenceToSmiles(sequence: string, inverted: boolean = false, format: string): string {
  const obj = getObjectWithCodesAndSmiles(sequence, format);
  let codes = sortByStringLengthInDescendingOrder(Object.keys(obj));
  let i = 0;
  let smiles = '';
  const codesList = [];
  const links = ['s', 'ps', '*'];
  const includesStandardLinkAlready = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns).concat(delimiter);
  while (i < sequence.length) {
    const code = codes.find((s: string) => s == sequence.slice(i, i + s.length))!;
    i += code.length;
    inverted ? codesList.unshift(code) : codesList.push(code);
  }
  for (let i = 0; i < codesList.length; i++) {
    if (dropdowns.includes(codesList[i])) {
      smiles += (i >= codesList.length / 2) ?
        MODIFICATIONS[codesList[i]].right + standardPhosphateLinkSmiles:
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
  obj[delimiter] = '';
  // TODO: create object based from synthesizer type to avoid key(codes) duplicates
  const output = isValidSequence(sequence, format);
  if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12))
    obj['g'] = map[SYNTHESIZERS.MERMADE_12][TECHNOLOGIES.SI_RNA]['g'].SMILES;
  else if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS))
    obj['g'] = map[SYNTHESIZERS.AXOLABS][TECHNOLOGIES.SI_RNA]['g'].SMILES;
  return obj;
}

// todo: remove
// const NAME = 'name';
const CODES = 'codes';
const SMILES = 'smiles';

function getObjectWithCodesAndSmilesFromFile(sequence: string, format: string, libFileContent: string) {
  const obj: { [code: string]: string } = {};
  // todo: type
  const lib: any[] = JSON.parse(libFileContent);

  if (format == null) {
    for (const item of lib) {
      for (const synthesizer of Object.keys(item[CODES])) {
        for (const technology of Object.keys(item[CODES][synthesizer])) {
          const codes = item[CODES][synthesizer][technology];
          let smiles: string = item[SMILES];
          smiles = smiles.replace(/\[OH:1\]/g, 'O');
          smiles = smiles.replace(/\[OH:2\]/g, 'O');
          for (const code of codes)
            obj[code] = smiles;
        }
      }
    }
  } else {
    for (const item of lib) {
      for (const synthesizer of Object.keys(item[CODES])) {
        if (synthesizer === format) {
          for (const technology of Object.keys(item[CODES][synthesizer])) {
            const codes = item[CODES][synthesizer][technology];
            let smiles: string = item[SMILES];
            smiles = smiles.replace(/\[OH:1\]/g, 'O');
            smiles = smiles.replace(/\[OH:2\]/g, 'O');
            for (const code of codes)
              obj[code] = smiles;
          }
        }
      }
    }
  }
  // todo: don't forget to replace [OH:1] and ... by O in smiles strings

  // if (format == null) {
  //   for (const synthesizer of Object.keys(map)) {
  //     for (const technology of Object.keys(map[synthesizer])) {
  //       for (const code of Object.keys(map[synthesizer][technology]))
  //         obj[code] = map[synthesizer][technology][code].SMILES;
  //     }
  //   }
  // } else {
  //   for (const technology of Object.keys(map[format])) {
  //     for (const code of Object.keys(map[format][technology]))
  //       obj[code] = map[format][technology][code].SMILES;
  //   }
  // }
  obj[delimiter] = '';
  // TODO: create object based from synthesizer type to avoid key(codes) duplicates
  const output = isValidSequence(sequence, format);
  if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12)) {
    // todo: remove as quickfix, optimize access to 'g'
    for (const item of lib) {
      for (const synthesizer of Object.keys(item[CODES])) {
        for (const technology of Object.keys(item[CODES][synthesizer])) {
          const codes = item[CODES][synthesizer][technology];
          for (const code of codes) {
            const condition =
              (code === 'g') &&
              (synthesizer === SYNTHESIZERS.MERMADE_12) &&
              (technology === TECHNOLOGIES.SI_RNA);
            if (condition) {
              let smiles: string = item[SMILES];
              smiles = smiles.replace(/\[OH:1\]/g, 'O');
              smiles = smiles.replace(/\[OH:2\]/g, 'O');
              obj[code] = smiles;
            }
          }
        }
      }
    }
  } else if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS)) {
    for (const item of lib) {
      for (const synthesizer of Object.keys(item[CODES])) {
        for (const technology of Object.keys(item[CODES][synthesizer])) {
          const codes = item[CODES][synthesizer][technology];
          for (const code of codes) {
            const condition =
              (code === 'g') &&
              (synthesizer === SYNTHESIZERS.AXOLABS) &&
              (technology === TECHNOLOGIES.SI_RNA);
            if (condition) {
              let smiles: string = item[SMILES];
              smiles = smiles.replace(/\[OH:1\]/g, 'O');
              smiles = smiles.replace(/\[OH:2\]/g, 'O');
              obj[code] = smiles;
            }
          }
        }
      }
    }
  }
  return obj;
}

function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a: string, b: string) {return b.length - a.length;});
}
