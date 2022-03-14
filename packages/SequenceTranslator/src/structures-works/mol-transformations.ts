import * as OCL from 'openchemlib/full.js';

export function getMol(smilesCodes: string[]) {
  const molBlocks: string[] = [];

  for (let i = 0; i < smilesCodes.length; i++)
    molBlocks.push(rotateNucleotidesV3000(smilesCodes[i]));

  return linkV3000(molBlocks);
}

function linkV3000(molBlocks: string[]) {
  let macroMolBlock = '\nDatagrok macromolecule handler\n\n';
  macroMolBlock += '  0  0  0  0  0  0              0 V3000\n';
  macroMolBlock += 'M  V30 BEGIN CTAB\n';
  let atomBlock = '';
  let bondBlock = '';
  //let collectionBlock = '';
  let natom = 0;
  let nbond = 0;
  let xShift = 0;

  for (let i = 0; i < molBlocks.length; i++) {
    const numbers = extractAtomsBondsNumbersV3000(molBlocks[i]);
    const coordinates = extractAtomDataV3000(molBlocks[i]);
    let indexAtoms = molBlocks[i].indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
    indexAtoms = molBlocks[i].indexOf('\n', indexAtoms);
    let index = indexAtoms;
    let indexEnd = indexAtoms;

    for (let j = 0; j < numbers.natom; j++) {
      if (coordinates.atomIndex[j] != 1 || i == 0) {
        //rewrite atom number
        index = molBlocks[i].indexOf('V30', index) + 4;
        indexEnd = molBlocks[i].indexOf(' ', index);
        const atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
        molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);

        //rewrite coordinates
        index = molBlocks[i].indexOf(' ', index) + 1;
        index = molBlocks[i].indexOf(' ', index) + 1;
        indexEnd = molBlocks[i].indexOf(' ', index);
        const coordinate = parseFloat(molBlocks[i].substring(index, indexEnd)) + xShift;
        molBlocks[i] = molBlocks[i].slice(0, index) + coordinate + molBlocks[i].slice(indexEnd);

        index = molBlocks[i].indexOf('\n', index) + 1;
      } else {
        index = molBlocks[i].indexOf('M  V30', index) - 1;
        indexEnd = molBlocks[i].indexOf('\n', index + 1);
        molBlocks[i] = molBlocks[i].slice(0, index) + molBlocks[i].slice(indexEnd);
      }
    }

    const indexAtomsEnd = molBlocks[i].indexOf('M  V30 END ATOM');
    atomBlock += molBlocks[i].substring(indexAtoms + 1, indexAtomsEnd);

    let indexBonds = molBlocks[i].indexOf('M  V30 BEGIN BOND'); // V3000 index for bonds
    indexBonds = molBlocks[i].indexOf('\n', indexBonds);
    index = indexBonds;
    indexEnd = indexBonds;

    for (let j = 0; j < numbers.nbond; j++) {
      //rewrite bond number
      index = molBlocks[i].indexOf('V30', index) + 4;
      indexEnd = molBlocks[i].indexOf(' ', index);
      const bondNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + nbond;
      molBlocks[i] = molBlocks[i].slice(0, index) + bondNumber + molBlocks[i].slice(indexEnd);

      //rewrite atom pair in bond
      index = molBlocks[i].indexOf(' ', index) + 1;
      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = molBlocks[i].indexOf(' ', index);
      let atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
      molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);
      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = molBlocks[i].indexOf('\n', index);
      atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
      molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);

      index = molBlocks[i].indexOf('\n', index) + 1;
    }

    const indexBondEnd = molBlocks[i].indexOf('M  V30 END BOND');
    bondBlock += molBlocks[i].substring(indexBonds + 1, indexBondEnd);

    natom += numbers.natom - 1;
    nbond += numbers.nbond;
    xShift += coordinates.x[numbers.natom - 1];
  }

  natom++;
  macroMolBlock += 'M  V30 COUNTS ' + natom + ' ' + nbond + ' 0 0 0\n';
  macroMolBlock += 'M  V30 BEGIN ATOM\n';
  macroMolBlock += atomBlock;
  macroMolBlock += 'M  V30 END ATOM\n';
  macroMolBlock += 'M  V30 BEGIN BOND\n';
  macroMolBlock += bondBlock;
  macroMolBlock += 'M  V30 END BOND\n';
  macroMolBlock += 'M  V30 END CTAB\n';
  macroMolBlock += 'M  END\n';

  return macroMolBlock;
}

function rotateNucleotidesV3000(smiles: string) {
  let molBlock = OCL.Molecule.fromSmiles(smiles).toMolfileV3();
  const coordinates = extractAtomDataV3000(molBlock);
  const natom = coordinates.atomIndex.length;

  const indexFivePrime = coordinates.atomIndex.indexOf(1);
  const indexThreePrime = coordinates.atomIndex.indexOf(natom);

  //fix 5 prime if inadequate
  if (natom > 6)
    fix5Prime(coordinates, indexFivePrime, indexThreePrime);

  const xCenter = (coordinates.x[indexThreePrime] + coordinates.x[indexFivePrime])/2;
  const yCenter = (coordinates.y[indexThreePrime] + coordinates.y[indexFivePrime])/2;

  //place to center
  for (let i = 0; i < natom; i++) {
    coordinates.x[i] -= xCenter;
    coordinates.y[i] -= yCenter;
  }

  let angle = 0;
  if (coordinates.x[indexFivePrime] == 0)
    angle = coordinates.y[indexFivePrime] > coordinates.y[indexThreePrime] ? Math.PI/2 : 3*Math.PI/2;
  else if (coordinates.y[indexFivePrime] == 0)
    angle = coordinates.x[indexFivePrime] > coordinates.x[indexThreePrime] ? Math.PI : 0;
  else {
    const derivative = coordinates.y[indexFivePrime]/coordinates.x[indexFivePrime];
    angle = derivative > 0 ? Math.PI - Math.atan(derivative) : Math.atan(derivative);
  }

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  for (let i = 0; i < natom; i++) {
    const xAdd = coordinates.x[i];
    coordinates.x[i] = xAdd*cos - coordinates.y[i]*sin;
    coordinates.y[i] = xAdd*sin + coordinates.y[i]*cos;
  }

  //place to right
  const xShift = coordinates.x[indexFivePrime];
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

function fix5Prime(coordinates: {atomIndex: number[], atomType: string[], x: number[], y: number[]},
  indexFivePrime: number, indexThreePrime: number) {
  const indexFivePrimeNeighbour = indexFivePrime + 1;
  const xShift = coordinates.x[indexFivePrimeNeighbour];
  const yShift = coordinates.y[indexFivePrimeNeighbour];
  const base3PrimeX = coordinates.x[indexThreePrime] - xShift;
  const base3PrimeY = coordinates.y[indexThreePrime] - yShift;
  const base5PrimeX = coordinates.x[indexFivePrime] - xShift;
  const base5PrimeY = coordinates.y[indexFivePrime] - yShift;

  const rotated5PrimeX = base5PrimeX*Math.cos(Math.PI*2/3) - base5PrimeY*Math.sin(Math.PI*2/3);
  const rotated5PrimeY = base5PrimeX*Math.sin(Math.PI*2/3) + base5PrimeY*Math.cos(Math.PI*2/3);

  const dx = base5PrimeX - base3PrimeX;
  const dy = base5PrimeY - base3PrimeY;
  const dxRotated = rotated5PrimeX - base3PrimeX;
  const dyRotated = rotated5PrimeY - base3PrimeY;

  if (Math.sqrt(dyRotated*dyRotated + dxRotated*dxRotated) >= Math.sqrt(dy*dy + dx*dx)) {
    coordinates.x[indexFivePrime] = rotated5PrimeX + xShift;
    coordinates.y[indexFivePrime] = rotated5PrimeY + yShift;
  }
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
