import * as grok from 'datagrok-api/grok';

type Indices = {
  first: number, // node to which remFirst is attached
  last: number, // node to which remLast is attached
  remFirst: number, // "leftmost" r-group node of the monomer
  remLast: number, // "rightmost" r-group node of the monomer
  remBondFirst: number, // idx of the bond between first and remFirst
  remBondLast: number, // idx of the bond between last and remLast
}

// todo: improve
type AtomData = {
  atomIndices: number[],
  atomTypes: string[],
  x: number[],
  y: number[],
}

function extractAtomDataV3K(v3KMolblock: string): AtomData {
  const numbers = extractAtomAndBondCountsV3K(v3KMolblock);
  let begin = v3KMolblock.indexOf('M  V30 BEGIN ATOM'); // V3000 block for atom coordinates
  begin = v3KMolblock.indexOf('\n', begin);
  let end = begin;

  const atomIndices: number[] = Array(numbers.atomCount);
  const atomTypes: string[] = Array(numbers.atomCount);
  const x: number[] = Array(numbers.atomCount);
  const y: number[] = Array(numbers.atomCount);

  for (let i = 0; i < numbers.atomCount; i++) {
    begin = v3KMolblock.indexOf('V30', begin) + 4;
    end = v3KMolblock.indexOf(' ', begin);
    atomIndices[i] = parseInt(v3KMolblock.substring(begin, end));

    begin = end + 1;
    end = v3KMolblock.indexOf(' ', begin);
    atomTypes[i] = v3KMolblock.substring(begin, end);

    begin = end + 1;
    end = v3KMolblock.indexOf(' ', begin);
    x[i] = parseFloat(v3KMolblock.substring(begin, end));

    begin = end + 1;
    end = v3KMolblock.indexOf(' ', begin);
    y[i] = parseFloat(v3KMolblock.substring(begin, end));

    begin = v3KMolblock.indexOf('\n', begin) + 1;
  }

  return {atomIndices: atomIndices, atomTypes: atomTypes, x: x, y: y};
}

function extractAtomAndBondCountsV3K(v3KMolblock: string): {atomCount: number, bondCount: number} {
  v3KMolblock = v3KMolblock.replaceAll('\r', ''); // equalize old and new sdf standards

  // parse atom count
  let idxBegin = v3KMolblock.indexOf('COUNTS') + 7;
  let idxEnd = v3KMolblock.indexOf(' ', idxBegin);
  const numOfAtoms = parseInt(v3KMolblock.substring(idxBegin, idxEnd));

  // parse bond count
  idxBegin = idxEnd + 1;
  idxEnd = v3KMolblock.indexOf(' ', idxBegin);
  const numOfBonds = parseInt(v3KMolblock.substring(idxBegin, idxEnd));

  return {atomCount: numOfAtoms, bondCount: numOfBonds};
}

function getIndices(v2KMolblock: string, v3KMolblock: string): Indices {
  // one should also take into account that there can be multiple M RGP lines
  // with R-Groups specified
  let begin = v2KMolblock.indexOf('M  RGP', 0) + 8;
  let end = v2KMolblock.indexOf('\n', begin);

  // todo: maybe this part deserves a separate function
  // there may be situation when the rgp information is distributed among
  // multiple lines, this must be taken into account
  const rgpStringParsed = v2KMolblock.substring(begin, end).replaceAll('  ', ' ').replaceAll('  ', ' ').split(' ');
  const rgpData = rgpStringParsed.map((el) => parseInt(el));
  // rgpData[0] is the number of R-groups
  // the following code sets remFirst to node# to which RGP#1 is substituted

  // todo: handle the exceptional case when there is not enough rgroups
  // todo: handle the exceptional case when the order is different
  // rgpData[1] is the node to which rgp #1 gets substituted
  const remFirst = rgpData[2] == 1 ? rgpData[1] : -1;
  // rgpData[3] is the node to which rgp #2 gets substituted
  const remLast = rgpData[4] == 2 ? rgpData[3] : -1;
  if (remFirst === -1 || remLast === -1)
    throw new Error('RGP parsing: first and last groups have wrong format');

  // const remFirst = rgpData[2] == '1' ? parseInt(rgpData[1]) : parseInt(rgpData[3]);
  // const remLast = rgpData[2] == '2' ? parseInt(rgpData[1]) : parseInt(rgpData[3]);

  // todo: rename 'numbers'
  const numbers = extractAtomAndBondCountsV3K(v3KMolblock);
  let indexBonds = v3KMolblock.indexOf('M  V30 BEGIN BOND'); // V3000 bond block
  indexBonds = v3KMolblock.indexOf('\n', indexBonds);
  begin = indexBonds;
  end = indexBonds;

  let first = 0; // todo: improve notation
  let last = 0;
  let remBondFirst = 0;
  let remBondLast = 0;

  // iterate over edges of the graph and find those
  for (let j = 0; j < numbers.bondCount; j++) {
    if (first === 0 || last === 0) {
      begin = v3KMolblock.indexOf('V30', begin) + 4;
      end = v3KMolblock.indexOf('\n', begin);
      const bondStringParsed = v3KMolblock.substring(begin, end).replaceAll('  ', ' ').replaceAll('  ', ' ').split(' ');
      const bondData = bondStringParsed.map((el) => parseInt(el));

      if (bondData[2] === remFirst) { // bondData[2] is the 1st node/atom of the bond
        first = bondData[3]; // bondData[3] is the 2nd node/atom of the bond
        remBondFirst = bondData[0]; // bondData[0] is the idx of the associated bond/edge
      } else if (bondData[3] === remFirst) {
        first = bondData[2];
        remBondFirst = bondData[0];
      } else if (bondData[2] === remLast) {
        last = bondData[3];
        remBondLast = bondData[0];
      } else if (bondData[3] === remLast) {
        last = bondData[2];
        remBondLast = bondData[0];
      }
    }
  }

  return {first, last, remFirst, remLast, remBondFirst, remBondLast};
}

/* provide description */
async function rotateBackboneV3K(v3KMolblock: string, indices:any): Promise<string> {
  // todo: rename 'coordinates'
  const coordinates = extractAtomDataV3K(v3KMolblock);
  const atomCount = coordinates.atomIndices.length;
  const first = indices['first'];
  const last = indices['last'];

  const xCenter = (coordinates.x[last] + coordinates.x[first])/2;
  const yCenter = (coordinates.y[last] + coordinates.y[first])/2;

  //place to center
  for (let i = 0; i < atomCount; i++) {
    coordinates.x[i] -= xCenter;
    coordinates.y[i] -= yCenter;
  }

  let angle = 0;
  if (coordinates.x[first] == 0) {
    angle = coordinates.y[first] > coordinates.y[last] ? Math.PI/2 : 3*Math.PI/2;
  } else if (coordinates.y[first] == 0) {
    angle = coordinates.x[first] > coordinates.x[last] ? Math.PI : 0;
  } else {
    const tangent = coordinates.y[first]/coordinates.x[first];
    if (coordinates.x[first] < coordinates.x[last])
      angle = tangent > 0 ? Math.PI - Math.atan(tangent) : Math.atan(tangent);
    else
      angle = tangent > 0 ? Math.atan(tangent) : Math.PI - Math.atan(tangent);
  }

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  for (let i = 0; i < atomCount; i++) {
    const xAdd = coordinates.x[i];
    coordinates.x[i] = xAdd*cos - coordinates.y[i]*sin;
    coordinates.y[i] = xAdd*sin + coordinates.y[i]*cos;
  }

  //place to right
  const xShift = coordinates.x[first];
  for (let i = 0; i < atomCount; i++)
    coordinates.x[i] -= xShift;

  //rewrite v3KMolblock
  let index = v3KMolblock.indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
  index = v3KMolblock.indexOf('\n', index);
  let indexEnd = index;
  for (let i = 0; i < atomCount; i++) {
    index = v3KMolblock.indexOf('V30', index) + 4;
    index = v3KMolblock.indexOf(' ', index) + 1;
    index = v3KMolblock.indexOf(' ', index) + 1;
    indexEnd = v3KMolblock.indexOf(' ', index) + 1;
    indexEnd = v3KMolblock.indexOf(' ', indexEnd);

    v3KMolblock = v3KMolblock.slice(0, index) +
      coordinates.x[i] + ' ' + coordinates.y[i] +
      v3KMolblock.slice(indexEnd);

    index = v3KMolblock.indexOf('\n', index) + 1;
  }

  return v3KMolblock;
}

/* provide description */
function linkV3K(monomers: any[]): string {
  let macroMolBlock = '\nDatagrok macromolecule handler\n\n';
  macroMolBlock += '  0  0  0  0  0  0            999 V3000\n';
  macroMolBlock += 'M  V30 BEGIN CTAB\n';
  let atomBlock = '';
  let bondBlock = '';
  // const collectionBlock = '';
  // const collection: number [] = [];
  let atomCount = 0;
  let bondCount = 0;
  let xShift = 0; // ?

  for (let i = 0; i < monomers.length; i++) {
    let v3KMolfile = monomers[i]['molfile'];
    const first = monomers[i]['indices']['first'];
    const last = monomers[i]['indices']['last'];
    const remFirst = monomers[i]['indices']['remFirst'];
    const remLast = monomers[i]['indices']['remLast'];
    const remBondFirst = monomers[i]['indices']['remBondFirst'];
    const remBondLast = monomers[i]['indices']['remBondLast'];
    v3KMolfile = v3KMolfile.replaceAll('(-\nM  V30 ', '(')
      .replaceAll('-\nM  V30 ', '').replaceAll(' )', ')');
    // todo: improve naming
    const numbers = extractAtomAndBondCountsV3K(v3KMolfile);
    const coordinates = extractAtomDataV3K(v3KMolfile);

    let indexAtoms = v3KMolfile.indexOf('M  V30 BEGIN ATOM'); // V3000 atom block
    indexAtoms = v3KMolfile.indexOf('\n', indexAtoms);
    let index = indexAtoms;
    let indexEnd = indexAtoms;
    const totalShift = xShift - coordinates.x[first - 1];

    for (let j = 0; j < numbers.atomCount; j++) {
      if (coordinates.atomIndices[j] != remFirst && coordinates.atomIndices[j] != remLast) { //|| i == 0) {
        //rewrite atom number
        index = v3KMolfile.indexOf('V30', index) + 4;
        indexEnd = v3KMolfile.indexOf(' ', index);

        let atomNumber = parseInt(v3KMolfile.substring(index, indexEnd));
        atomNumber = (atomNumber > remFirst && atomNumber > remLast) ? atomNumber - 2 :
          (atomNumber > remFirst || atomNumber > remLast) ? atomNumber - 1 : atomNumber;
        atomNumber += atomCount;
        v3KMolfile = v3KMolfile.slice(0, index) + atomNumber + v3KMolfile.slice(indexEnd);

        //rewrite coordinates
        index = v3KMolfile.indexOf(' ', index) + 1;
        index = v3KMolfile.indexOf(' ', index) + 1;
        indexEnd = v3KMolfile.indexOf(' ', index);

        let coordinate = Math.round(10000*(parseFloat(v3KMolfile.substring(index, indexEnd)) + totalShift))/10000;
        v3KMolfile = v3KMolfile.slice(0, index) + coordinate + v3KMolfile.slice(indexEnd);

        index = v3KMolfile.indexOf(' ', index) + 1;
        indexEnd = v3KMolfile.indexOf(' ', index);
        coordinate = Math.round(10000*(parseFloat(v3KMolfile.substring(index, indexEnd))))/10000;
        v3KMolfile = v3KMolfile.slice(0, index) + coordinate + v3KMolfile.slice(indexEnd);

        index = v3KMolfile.indexOf('\n', index) + 1;
      } else {
        index = v3KMolfile.indexOf('M  V30', index) - 1;
        indexEnd = v3KMolfile.indexOf('\n', index + 1);
        v3KMolfile = v3KMolfile.slice(0, index) + v3KMolfile.slice(indexEnd);
      }
    }

    const indexAtomsEnd = v3KMolfile.indexOf('M  V30 END ATOM');
    atomBlock += v3KMolfile.substring(indexAtoms + 1, indexAtomsEnd);

    let indexBonds = v3KMolfile.indexOf('M  V30 BEGIN BOND'); // V3000 index for bonds
    indexBonds = v3KMolfile.indexOf('\n', indexBonds);
    index = indexBonds;
    indexEnd = indexBonds;
    let bondNumber = 0;

    for (let j = 0; j < numbers.bondCount; j++) {
      //rewrite bond number
      index = v3KMolfile.indexOf('V30', index) + 4;
      indexEnd = v3KMolfile.indexOf(' ', index);
      bondNumber = parseInt(v3KMolfile.substring(index, indexEnd));

      if (bondNumber == remBondFirst || bondNumber == remBondLast) {
        indexEnd = v3KMolfile.indexOf('\n', index) + 1;
        index -=7;
        v3KMolfile = v3KMolfile.slice(0, index) + v3KMolfile.slice(indexEnd);
        continue;
      }

      bondNumber = (bondNumber > remBondFirst && bondNumber > remBondLast) ? bondNumber - 2 :
        (bondNumber > remBondFirst || bondNumber > remBondLast) ? bondNumber - 1 : bondNumber;
      bondNumber += bondCount;

      v3KMolfile = v3KMolfile.slice(0, index) + bondNumber + v3KMolfile.slice(indexEnd);

      //rewrite atom pair in bond
      index = v3KMolfile.indexOf(' ', index) + 1;
      index = v3KMolfile.indexOf(' ', index) + 1;
      indexEnd = v3KMolfile.indexOf(' ', index);
      let atomNumber = parseInt(v3KMolfile.substring(index, indexEnd));
      atomNumber = (atomNumber > remFirst && atomNumber > remLast) ? atomNumber - 2 :
        (atomNumber > remFirst || atomNumber > remLast) ? atomNumber - 1 : atomNumber;
      atomNumber += atomCount;
      v3KMolfile = v3KMolfile.slice(0, index) + atomNumber + v3KMolfile.slice(indexEnd);
      index = v3KMolfile.indexOf(' ', index) + 1;
      indexEnd = Math.min(v3KMolfile.indexOf('\n', index), v3KMolfile.indexOf(' ', index));
      atomNumber = parseInt(v3KMolfile.substring(index, indexEnd));
      atomNumber = (atomNumber > remFirst && atomNumber > remLast) ? atomNumber - 2 :
        (atomNumber > remFirst || atomNumber > remLast) ? atomNumber - 1 : atomNumber;
      atomNumber += atomCount;
      v3KMolfile = v3KMolfile.slice(0, index) + atomNumber + v3KMolfile.slice(indexEnd);

      index = v3KMolfile.indexOf('\n', index) + 1;
    }

    const indexBondEnd = v3KMolfile.indexOf('M  V30 END BOND');
    bondBlock += v3KMolfile.substring(indexBonds + 1, indexBondEnd);
    //let indexCollection = v3KMolfile.indexOf('M  V30 MDLV30/STEABS ATOMS=('); // V3000 index for collections

    // while (indexCollection != -1) {
    //   indexCollection += 28;
    //   const collectionEnd = v3KMolfile.indexOf(')', indexCollection);
    //   const collectionEntries = v3KMolfile.substring(indexCollection, collectionEnd).split(' ').slice(1);
    //   collectionEntries.forEach((e: string) => {
    //     collection.push(parseInt(e) + atomCount);
    //   });
    //   indexCollection = collectionEnd;
    //   indexCollection = v3KMolfile.indexOf('M  V30 MDLV30/STEABS ATOMS=(', indexCollection);
    // }

    atomCount += numbers.atomCount - 2;
    bondCount += numbers.bondCount - 2;
    xShift += coordinates.x[last] - coordinates.x[first] + 1;

    if (i == monomers.length -1) {
      atomCount++;
      const shift = xShift + 0.2;
      atomBlock += 'M  V30 ' + atomCount + ' O ' + shift + ' 0 0.000000 0\n';
    }
    bondCount++;
    if (i == monomers.length -1) {
      const rightTerminal = (last > remFirst && last > remLast) ? last + atomCount - (numbers.atomCount - 2) - 3:
        (last > remFirst || last > remLast) ? last + atomCount - (numbers.atomCount - 2) - 2 :
          last + atomCount - (numbers.atomCount - 2) - 1;
      bondBlock += 'M  V30 ' + bondCount + ' 1 ' + rightTerminal + ' ' + atomCount + '\n';
    } else {
      const rightTerminal = (last > remFirst && last > remLast) ? last + atomCount - (numbers.atomCount - 2) - 2:
        (last > remFirst || last > remLast) ? last + atomCount - (numbers.atomCount - 2) - 1 :
          last + atomCount - (numbers.atomCount - 2);

      const next = monomers[i + 1]['indices'];
      const nextFirst = next['first'];
      const nextRemFirst = next['remFirst'];
      const nextRemLast = next['remLast'];

      const leftTerminal = (nextFirst > nextRemFirst && nextFirst > nextRemLast) ? nextFirst + atomCount - 2 :
        (nextFirst > nextRemFirst || nextFirst > nextRemLast) ? nextFirst + atomCount - 1 :
          nextFirst + atomCount;

      bondBlock += 'M  V30 ' + bondCount + ' 1 ' + rightTerminal + ' ' + leftTerminal + '\n';
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
  macroMolBlock += 'M  V30 COUNTS ' + atomCount + ' ' + bondCount + ' 0 0 0\n';
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

export async function getMacroMol(monomers: any[][]): Promise<string[]> {
  const result: string[] = [];
  const moduleRdkit = await grok.functions.call('Chem:getRdKitModule');
  for (let i = 0; i < monomers.length; i++) {
    for (let j = 0; j < monomers[i].length; j++) {
      const molblock = moduleRdkit.get_mol(monomers[i][j]['molfile']); // V2000
      const v3KMolblock = molblock.get_v3Kmolblock();
      const indices = getIndices(molblock, v3KMolblock);
      monomers[i][j]['indices'] = indices;

      // a new molfile for 'rotated' is obtained in v3k format
      monomers[i][j]['molfile'] = await rotateBackboneV3K(v3KMolblock, indices);
      molblock?.delete();
    }
    // seemingly, at this stage the bond is reconstructed
    result.push(linkV3K(monomers[i]));
  }

  return result;
}
