import * as OCL from 'openchemlib/full.js';
import * as grok from 'datagrok-api/grok';

import {RDModule} from '@datagrok-libraries/bio/src/rdkit-api';

export async function getMacroMol(monomers: any[][]): Promise<string[]> {
  let result: string[] = [];
  const moduleRdkit = await grok.functions.call('Chem:getRdKitModule');
  for(let i = 0; i < monomers.length; i++) {
    for (let j = 0; j < monomers[i].length; j++){
      const mol = moduleRdkit.get_mol(monomers[i][j]['molfile']);
      const a = mol.get_v3Kmolblock();
      const indices = getIndices(monomers[i][j], a);
      monomers[i][j]['indices'] = indices;

      monomers[i][j]['molfile'] = await rotateBackboneV3000(a, indices);
      mol?.delete();
    }
    result.push(linkV3000(monomers[i]));
  }

  return result;
}

function getIndices(monomer: any, molV3000: string): {first: number, last: number, 
                                    remFirst: number, remLast: number, 
                                    remBondFirst: number, remBondLast: Number} {
  const molfile = monomer["molfile"];
  let indexStart = molfile.indexOf('M  RGP', 0) + 8;
  let indexEnd = molfile.indexOf('\n', indexStart);
  const indicesData = molfile.substring(indexStart, indexEnd).replaceAll('  ', ' ').replaceAll('  ', ' ');
  let parsedData = indicesData.split(' ')
  const remFirst = parsedData[2] == '1' ? parseInt(parsedData[1]) : parseInt(parsedData[3]);
  const remLast = parsedData[2] == '2' ? parseInt(parsedData[1]) : parseInt(parsedData[3]);
  
  const numbers = extractAtomsBondsNumbersV3000(molV3000);
  let indexBonds = molV3000.indexOf('M  V30 BEGIN BOND'); // V3000 index for bonds
  indexBonds = molV3000.indexOf('\n', indexBonds);
  indexStart = indexBonds;
  indexEnd = indexBonds;

  let first = 0;
  let last = 0;
  let remBondFirst = 0;
  let remBondLast = 0;

  for (let j = 0; j < numbers.nbond; j++) {
    if(first == 0 || last == 0){
      indexStart = molV3000.indexOf('V30', indexStart) + 4;
      indexEnd = molV3000.indexOf('\n', indexStart);
      const bondData = molV3000.substring(indexStart, indexEnd).replaceAll('  ', ' ').replaceAll('  ', ' ');
      parsedData = bondData.split(' ')

      if(parseInt(parsedData[2]) == remFirst){
        first = parseInt(parsedData[3]);
        remBondFirst = parseInt(parsedData[0]);
      }
      else if(parseInt(parsedData[3]) == remFirst){
        first = parseInt(parsedData[2]);
        remBondFirst = parseInt(parsedData[0]);
      }
      else if(parseInt(parsedData[2]) == remLast){
        last = parseInt(parsedData[3]);
        remBondLast = parseInt(parsedData[0]);
      }
      else if(parseInt(parsedData[3]) == remLast){
        last = parseInt(parsedData[2]);
        remBondLast = parseInt(parsedData[0]);
      }
    }
  }

  return {first, last, remFirst, remLast, remBondFirst, remBondLast};
}

async function rotateBackboneV3000(molBlock: string, indices:any): Promise<string> {
  const coordinates = extractAtomDataV3000(molBlock);
  const natom = coordinates.atomIndex.length;
  const first = indices['first'];
  const last = indices['last'];

  const xCenter = (coordinates.x[last] + coordinates.x[first])/2;
  const yCenter = (coordinates.y[last] + coordinates.y[first])/2;

  //place to center
  for (let i = 0; i < natom; i++) {
    coordinates.x[i] -= xCenter;
    coordinates.y[i] -= yCenter;
  }

  let angle = 0;
  if (coordinates.x[first] == 0)
    angle = coordinates.y[first] > coordinates.y[last] ? Math.PI/2 : 3*Math.PI/2;
  else if (coordinates.y[first] == 0)
    angle = coordinates.x[first] > coordinates.x[last] ? Math.PI : 0;
  else {
    const derivative = coordinates.y[first]/coordinates.x[first];
    if(coordinates.x[first] < coordinates.x[last])
      angle = derivative > 0 ? Math.PI - Math.atan(derivative) : Math.atan(derivative);
    else
      angle = derivative > 0 ? Math.atan(derivative) : Math.PI - Math.atan(derivative);
  }

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  for (let i = 0; i < natom; i++) {
    const xAdd = coordinates.x[i];
    coordinates.x[i] = xAdd*cos - coordinates.y[i]*sin;
    coordinates.y[i] = xAdd*sin + coordinates.y[i]*cos;
  }

  //place to right
  const xShift = coordinates.x[first];
  for (let i = 0; i < natom; i++)
    coordinates.x[i] -= xShift;

  //rewrite molBlock
  let index = molBlock.indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
  index = molBlock.indexOf('\n', index);
  let indexEnd = index;
  for (let i = 0; i < natom; i++) {
    index = molBlock.indexOf('V30', index) + 4;
    index = molBlock.indexOf(' ', index) + 1;
    index = molBlock.indexOf(' ', index) + 1;
    indexEnd = molBlock.indexOf(' ', index) + 1;
    indexEnd = molBlock.indexOf(' ', indexEnd);

    molBlock = molBlock.slice(0, index) +
      coordinates.x[i] + ' ' + coordinates.y[i] +
      molBlock.slice(indexEnd);

    index = molBlock.indexOf('\n', index) + 1;
  }

  return molBlock;
}

function extractAtomDataV3000(molBlock: string) {
  const numbers = extractAtomsBondsNumbersV3000(molBlock);
  let index = molBlock.indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
  index = molBlock.indexOf('\n', index);
  let indexEnd = index;

  const indexes: number[] = Array(numbers.natom);
  const types: string[] = Array(numbers.natom);
  const x: number[] = Array(numbers.natom);
  const y: number[] = Array(numbers.natom);

  for (let i = 0; i < numbers.natom; i++) {
    index = molBlock.indexOf('V30', index) + 4;
    indexEnd = molBlock.indexOf(' ', index);
    indexes[i] = parseInt(molBlock.substring(index, indexEnd));

    index = indexEnd + 1;
    indexEnd = molBlock.indexOf(' ', index);
    types[i] = molBlock.substring(index, indexEnd);

    index = indexEnd + 1;
    indexEnd = molBlock.indexOf(' ', index);
    x[i] = parseFloat(molBlock.substring(index, indexEnd));

    index = indexEnd + 1;
    indexEnd = molBlock.indexOf(' ', index);
    y[i] = parseFloat(molBlock.substring(index, indexEnd));

    index = molBlock.indexOf('\n', index) + 1;
  }

  return {atomIndex: indexes, atomType: types, x: x, y: y};
}

function linkV3000(monomers: any[]): string {
  let macroMolBlock = '\nDatagrok macromolecule handler\n\n';
  macroMolBlock += '  0  0  0  0  0  0            999 V3000\n';
  macroMolBlock += 'M  V30 BEGIN CTAB\n';
  let atomBlock = '';
  let bondBlock = '';
  let collectionBlock = '';
  const collection: number [] = [];
  let natom = 0;
  let nbond = 0;
  let xShift = 0;

  for (let i = 0; i < monomers.length; i++) {
    let molfile = monomers[i]['molfile'];
    const first = monomers[i]['indices']['first'];
    const last = monomers[i]['indices']['last'];
    const remFirst = monomers[i]['indices']['remFirst'];
    const remLast = monomers[i]['indices']['remLast'];
    const remBondFirst = monomers[i]['indices']['remBondFirst'];
    const remBondLast = monomers[i]['indices']['remBondLast'];
    molfile = molfile.replaceAll('(-\nM  V30 ', '(')
      .replaceAll('-\nM  V30 ', '').replaceAll(' )', ')');
    const numbers = extractAtomsBondsNumbersV3000(molfile);
    const coordinates = extractAtomDataV3000(molfile);

    let indexAtoms = molfile.indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
    indexAtoms = molfile.indexOf('\n', indexAtoms);
    let index = indexAtoms;
    let indexEnd = indexAtoms;
    const totalShift = xShift - coordinates.x[first - 1];

    for (let j = 0; j < numbers.natom; j++) {
      if (coordinates.atomIndex[j] != remFirst && coordinates.atomIndex[j] != remLast) { //|| i == 0) {
        //rewrite atom number
        index = molfile.indexOf('V30', index) + 4;
        indexEnd = molfile.indexOf(' ', index);
  
        let atomNumber = parseInt(molfile.substring(index, indexEnd))
        atomNumber = (atomNumber > remFirst && atomNumber > remLast) ? atomNumber - 2 :
                     (atomNumber > remFirst || atomNumber > remLast) ? atomNumber - 1 : atomNumber;
        atomNumber += natom;
        molfile = molfile.slice(0, index) + atomNumber + molfile.slice(indexEnd);

        //rewrite coordinates
        index = molfile.indexOf(' ', index) + 1;
        index = molfile.indexOf(' ', index) + 1;
        indexEnd = molfile.indexOf(' ', index);

        let coordinate = Math.round(10000*(parseFloat(molfile.substring(index, indexEnd)) + totalShift))/10000;
        molfile = molfile.slice(0, index) + coordinate + molfile.slice(indexEnd);

        index = molfile.indexOf(' ', index) + 1;
        indexEnd = molfile.indexOf(' ', index);
        coordinate = Math.round(10000*(parseFloat(molfile.substring(index, indexEnd))))/10000;
        molfile = molfile.slice(0, index) + coordinate + molfile.slice(indexEnd);

        index = molfile.indexOf('\n', index) + 1;
      } else {
        index = molfile.indexOf('M  V30', index) - 1;
        indexEnd = molfile.indexOf('\n', index + 1);
        molfile = molfile.slice(0, index) + molfile.slice(indexEnd);
      }
    }

    const indexAtomsEnd = molfile.indexOf('M  V30 END ATOM');
    atomBlock += molfile.substring(indexAtoms + 1, indexAtomsEnd);

    let indexBonds = molfile.indexOf('M  V30 BEGIN BOND'); // V3000 index for bonds
    indexBonds = molfile.indexOf('\n', indexBonds);
    index = indexBonds;
    indexEnd = indexBonds;
    let bondNumber = 0;

    for (let j = 0; j < numbers.nbond; j++) {
      //rewrite bond number
      index = molfile.indexOf('V30', index) + 4;
      indexEnd = molfile.indexOf(' ', index);
      bondNumber = parseInt(molfile.substring(index, indexEnd));

      if(bondNumber == remBondFirst || bondNumber == remBondLast){
        indexEnd = molfile.indexOf('\n', index) + 1;
        index -=7;
        molfile = molfile.slice(0, index) + molfile.slice(indexEnd);
        continue
      }

      bondNumber = (bondNumber > remBondFirst && bondNumber > remBondLast) ? bondNumber - 2 :
                   (bondNumber > remBondFirst || bondNumber > remBondLast) ? bondNumber - 1 : bondNumber;
      bondNumber += nbond;

      molfile = molfile.slice(0, index) + bondNumber + molfile.slice(indexEnd);

      //rewrite atom pair in bond
      index = molfile.indexOf(' ', index) + 1;
      index = molfile.indexOf(' ', index) + 1;
      indexEnd = molfile.indexOf(' ', index);
      let atomNumber = parseInt(molfile.substring(index, indexEnd))
      atomNumber = (atomNumber > remFirst && atomNumber > remLast) ? atomNumber - 2 :
                   (atomNumber > remFirst || atomNumber > remLast) ? atomNumber - 1 : atomNumber;
      atomNumber += natom;
      molfile = molfile.slice(0, index) + atomNumber + molfile.slice(indexEnd);
      index = molfile.indexOf(' ', index) + 1;
      indexEnd = Math.min(molfile.indexOf('\n', index), molfile.indexOf(' ', index));
      atomNumber = parseInt(molfile.substring(index, indexEnd))
      atomNumber = (atomNumber > remFirst && atomNumber > remLast) ? atomNumber - 2 :
                   (atomNumber > remFirst || atomNumber > remLast) ? atomNumber - 1 : atomNumber;
      atomNumber += natom;
      molfile = molfile.slice(0, index) + atomNumber + molfile.slice(indexEnd);

      index = molfile.indexOf('\n', index) + 1;
    }

    const indexBondEnd = molfile.indexOf('M  V30 END BOND');
    bondBlock += molfile.substring(indexBonds + 1, indexBondEnd);
    //let indexCollection = molfile.indexOf('M  V30 MDLV30/STEABS ATOMS=('); // V3000 index for collections

    // while (indexCollection != -1) {
    //   indexCollection += 28;
    //   const collectionEnd = molfile.indexOf(')', indexCollection);
    //   const collectionEntries = molfile.substring(indexCollection, collectionEnd).split(' ').slice(1);
    //   collectionEntries.forEach((e: string) => {
    //     collection.push(parseInt(e) + natom);
    //   });
    //   indexCollection = collectionEnd;
    //   indexCollection = molfile.indexOf('M  V30 MDLV30/STEABS ATOMS=(', indexCollection);
    // }

    natom += numbers.natom - 2;
    nbond += numbers.nbond - 2;
    xShift += coordinates.x[last] - coordinates.x[first] + 1;

    if(i == monomers.length -1){
      natom++;
      const shift = xShift + 0.2;
      atomBlock += 'M  V30 ' + natom + ' O ' + shift + ' 0 0.000000 0\n';
    }
    nbond++;
    if(i == monomers.length -1){
      const rightTerminal = (last > remFirst && last > remLast) ? last + natom - (numbers.natom - 2) - 3:
                            (last > remFirst || last > remLast) ? last + natom - (numbers.natom - 2) - 2 : 
                            last + natom - (numbers.natom - 2) - 1;
      bondBlock += 'M  V30 ' + nbond + ' 1 ' + rightTerminal + ' ' + natom + '\n';
    } else{
      const rightTerminal = (last > remFirst && last > remLast) ? last + natom - (numbers.natom - 2) - 2:
                            (last > remFirst || last > remLast) ? last + natom - (numbers.natom - 2) - 1 : 
                            last + natom - (numbers.natom - 2);
      
      const next = monomers[i + 1]['indices'];
      const nextFirst = next['first'];
      const nextRemFirst = next['remFirst'];
      const nextRemLast = next['remLast'];

      const leftTerminal = (nextFirst > nextRemFirst && nextFirst > nextRemLast) ? nextFirst + natom - 2 :
                           (nextFirst > nextRemFirst || nextFirst > nextRemLast) ? nextFirst + natom - 1 :
                            nextFirst + natom;

      bondBlock += 'M  V30 ' + nbond + ' 1 ' + rightTerminal + ' ' + leftTerminal + '\n';
    }
  }

  // const entries = 4;
  // const collNumber = Math.ceil(collection.length / entries);
  // collectionBlock += 'M  V30 MDLV30/STEABS ATOMS=(' + collection.length + ' -\n';
  // for (let i = 0; i < collNumber; i++) {
  //   collectionBlock += 'M  V30 ';
  //   const entriesCurrent = i + 1 == collNumber ? collection.length - (collNumber - 1)*entries : entries;
  //   for (let j = 0; j < entriesCurrent; j++) {
  //     collectionBlock += (j + 1 == entriesCurrent) ?
  //       (i == collNumber - 1 ? collection[entries*i + j] + ')\n' : collection[entries*i + j] + ' -\n') :
  //       collection[entries*i + j] + ' ';
  //   }
  // }

  //generate file
  macroMolBlock += 'M  V30 COUNTS ' + natom + ' ' + nbond + ' 0 0 0\n';
  macroMolBlock += 'M  V30 BEGIN ATOM\n';
  macroMolBlock += atomBlock;
  macroMolBlock += 'M  V30 END ATOM\n';
  macroMolBlock += 'M  V30 BEGIN BOND\n';
  macroMolBlock += bondBlock;
  macroMolBlock += 'M  V30 END BOND\n';
  //macroMolBlock += 'M  V30 BEGIN COLLECTION\n';
  //macroMolBlock += collectionBlock;
  //macroMolBlock += 'M  V30 END COLLECTION\n';
  macroMolBlock += 'M  V30 END CTAB\n';
  macroMolBlock += 'M  END\n';

  return macroMolBlock;
}

function extractAtomsBondsNumbersV3000(molBlock: string): {natom: number, nbond: number} {
  molBlock = molBlock.replaceAll('\r', ''); //equalize old and new sdf standards
  let index = molBlock.indexOf('COUNTS') + 7; // V3000 index for atoms and bonds number
  let indexEnd = molBlock.indexOf(' ', index);

  const atomsNumber = parseInt(molBlock.substring(index, indexEnd));
  index = indexEnd + 1;
  indexEnd = molBlock.indexOf(' ', index);
  const bondsNumber = parseInt(molBlock.substring(index, indexEnd));

  return {natom: atomsNumber, nbond: bondsNumber};
}

