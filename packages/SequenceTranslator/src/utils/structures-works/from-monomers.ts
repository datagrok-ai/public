// import {map, SYNTHESIZERS, TECHNOLOGIES, MODIFICATIONS, DELIMITER} from './map';
import {map, SYNTHESIZERS, TECHNOLOGIES, DELIMITER} from '../../hardcode-to-be-eliminated/map';
import {isValidSequence} from '../../sdf-tab/sequence-codes-tools';
import {sortByStringLengthInDescendingOrder} from '../helpers';
import {getMonomerWorks} from '../../package';
import {getNucleotidesMol} from './mol-transformations';

import {standardPhosphateLinkSmiles, MODIFICATIONS} from '../../hardcode-to-be-eliminated/const';
import {getMonomerLib} from '../../package';
// todo: remove
// const NAME = 'name';
const CODES = 'codes';
// const SMILES = 'smiles';
const MOL = 'molfile';

export function sequenceToMolV3000(
  sequence: string, inverted: boolean = false, oclRender: boolean = false,
  format: string
): string {
  const monomerNameFromCode = getCodeToNameMap(sequence, format);
  let codes = sortByStringLengthInDescendingOrder(Object.keys(monomerNameFromCode));
  let i = 0;
  const codesList = [];
  const links = ['s', 'ps', '*'];
  const includesStandardLinkAlready = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns).concat(DELIMITER);
  while (i < sequence.length) {
    const code = codes.find((s: string) => s === sequence.slice(i, i + s.length))!;
    i += code.length;
    inverted ? codesList.unshift(code) : codesList.push(code);
  }

  const monomers: string[] = [];

  for (let i = 0; i < codesList.length; i++) {
    if (links.includes(codesList[i]) ||
      includesStandardLinkAlready.includes(codesList[i]) ||
      (i < codesList.length - 1 && links.includes(codesList[i + 1]))
    ) {
      const aa = monomerNameFromCode[codesList[i]];
      if (aa !== undefined)
        monomers.push(aa);
      else
        monomers.push(codesList[i]);
    } else {
      const aa = monomerNameFromCode[codesList[i]];
      if (aa !== undefined)
        monomers.push(aa);
      else
        monomers.push(codesList[i]);
      monomers.push('p linkage');
    }
  }

  const lib = getMonomerLib();
  const mols: string [] = [];
  for (let i = 0; i < monomers.length; i++) {
    const mnmr = lib?.getMonomer('RNA', monomers[i]);
    mols.push(mnmr?.molfile!);
  }


  return getNucleotidesMol(mols);
  //return getMonomerWorks()?.getAtomicLevel(monomers, 'RNA')!;
}

export function sequenceToMolV3000_new(
  sequence: string, inverted: boolean = false, oclRender: boolean = false,
  format: string,
): string {
  const monomerNameFromCode = getCodeToNameMap(sequence, format);
  let codes = sortByStringLengthInDescendingOrder(Object.keys(monomerNameFromCode));
  let i = 0;
  const codesList = [];
  const links = ['s', 'ps', '*'];
  const includesStandardLinkAlready = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns).concat(DELIMITER);
  while (i < sequence.length) {
    const code = codes.find((s: string) => s === sequence.slice(i, i + s.length))!;
    i += code.length;
    inverted ? codesList.unshift(code) : codesList.push(code);
  }

  const monomers: string[] = [];

  for (let i = 0; i < codesList.length; i++) {
    if (links.includes(codesList[i]) ||
      includesStandardLinkAlready.includes(codesList[i]) ||
      (i < codesList.length - 1 && links.includes(codesList[i + 1]))
    )
      monomers.push(monomerNameFromCode[codesList[i]]);
    else {
      monomers.push(monomerNameFromCode[codesList[i]]);
      monomers.push('p linkage');
    }
  }

  return getMonomerWorks()?.getAtomicLevel(monomers, 'RNA')!;
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

function getCodeToNameMap(sequence: string, format: string) {
  const obj: { [code: string]: string } = {};
  const NAME = 'name';
  if (format == null) {
    for (const synthesizer of Object.keys(map)) {
      for (const technology of Object.keys(map[synthesizer])) {
        for (const code of Object.keys(map[synthesizer][technology]))
          obj[code] = map[synthesizer][technology][code][NAME]!;
      }
    }
  } else {
    for (const technology of Object.keys(map[format])) {
      for (const code of Object.keys(map[format][technology]))
        obj[code] = map[format][technology][code][NAME]!;
      // obj[code] = map[format][technology][code].SMILES;
    }
  }
  obj[DELIMITER] = '';
  // TODO: create object based from synthesizer type to avoid key(codes) duplicates
  const output = isValidSequence(sequence, format);
  if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12))
    obj['g'] = map[SYNTHESIZERS.MERMADE_12][TECHNOLOGIES.SI_RNA]['g'][NAME]!;
  else if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS))
    obj['g'] = map[SYNTHESIZERS.AXOLABS][TECHNOLOGIES.SI_RNA]['g'][NAME]!;
  return obj;
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

function getObjectWithCodesAndMolsFromFile(sequence: string, format: string, libFileContent: string) {
  const obj: { [code: string]: string } = {};
  // todo: type
  const lib: any[] = JSON.parse(libFileContent); //consider using library

  for (const item of lib) {
    for (const synthesizer of Object.keys(item[CODES])) {
      if (synthesizer === format) {
        for (const technology of Object.keys(item[CODES][synthesizer])) {
          const codes = item[CODES][synthesizer][technology];
          let mol: string = item[MOL];
          // todo: find another solution
          mol = mol.replace(/ R /g, ' O ');

          for (const code of codes)
            obj[code] = mol;
        }
      }
    }
  }

  obj[DELIMITER] = '';
  // TODO: create object based on synthesizer type to avoid key(codes) duplicates
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
              let mol: string = item[MOL];
              // todo: find another solution
              mol = mol.replace(/ R /g, ' O ');
              obj[code] = mol;
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
              let mol: string = item[MOL];
              // todo: find another solution
              mol = mol.replace(/ R /g, ' O ');
              obj[code] = mol;
            }
          }
        }
      }
    }
  }
  return obj;
}
