import * as DG from 'datagrok-api/dg';

import {randomInt} from '@datagrok-libraries/utils/src/random';
import {assert} from '@datagrok-libraries/utils/src/vector-operations';
import {SequenceGenerator} from './sequence-generator';

// eslint-disable-next-line max-len
const heavyColumns = 'Id,domain_no,hmm_species,chain_type,e-value,score,seqstart_index,seqend_index,identity_species,v_gene,v_identity,j_gene,j_identity,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,32A,33A,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,111A,111B,111C,111D,111E,112E,112D,112C,112B,112A,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128'.split(',');
// eslint-disable-next-line max-len
const lightColumns = 'Id,domain_no,hmm_species,chain_type,e-value,score,seqstart_index,seqend_index,identity_species,v_gene,v_identity,j_gene,j_identity,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,32A,33A,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112A,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127'.split(',');
const chainColumns = {
  'H': heavyColumns,
  'L': lightColumns,
};

function genChoices<T>(count: number, choices: T[]): T[] {
  const portions = choices.map(() => randomInt(count));
  const sum = portions.reduce((prev, curr) => prev + curr);
  const counts = portions.map((v) => v * count / sum);
  const vals: T[] = [];

  for (let i = 0; i < counts.length; ++i) {
    for (let j = 0; j < counts[i] && vals.length < count; ++j)
      vals.push(choices[i]);
  }
  return vals;
}

function genChainDataFrame(ids: string[], chainType: 'H' | 'L'): DG.DataFrame {
  const count = ids.length;
  const header = chainColumns[chainType];
  const columns: DG.Column[] = [];
  let index = 0;

  function addColumn(col: DG.Column) {
    assert(col.length == count, `${col.name} has wrong length (${col.length} != ${count}).`);
    columns.push(col);
  }

  addColumn(DG.Column.fromStrings(header[index++], ids)); // Id
  addColumn(DG.Column.fromInt32Array(header[index++], new Int32Array(count).fill(0))); // domain_no
  addColumn(DG.Column.fromStrings(header[index++], genChoices(count, [
    'alpaca',
    'mouse',
    'pig',
    'human',
    'rhesus',
  ]))); // hmm_species
  addColumn(DG.Column.fromStrings(header[index++], new Array<string>(count).fill(chainType))); // chain_type
  addColumn(DG.Column.fromFloat32Array(header[index++], new Float32Array(count).fill(0))); // e-value
  addColumn(DG.Column.fromFloat32Array(header[index++], new Float32Array(count).fill(123))); // score
  addColumn(DG.Column.fromInt32Array(header[index++], new Int32Array(count).fill(0))); // seqstart_index
  addColumn(DG.Column.fromInt32Array(header[index++], new Int32Array(count).fill(117))); // seqend_index
  addColumn(DG.Column.fromStrings(header[index++], new Array(count).fill(''))); // identity_species
  addColumn(DG.Column.fromStrings(header[index++], new Array(count).fill(''))); // v_gene
  addColumn(DG.Column.fromInt32Array(header[index++], new Int32Array(count).fill(0))); // v_identity
  addColumn(DG.Column.fromStrings(header[index++], new Array(count).fill(''))); // j_gene
  addColumn(DG.Column.fromInt32Array(header[index++], new Int32Array(count).fill(0))); // j_identity

  for (let i = index; i < header.length; ++i) {
    const gen = new SequenceGenerator(count);
    addColumn(DG.Column.fromStrings(header[i], gen.genList())); // pos
  }
  return DG.DataFrame.fromColumns(columns);
}

export function genData(ids: string[]): [DG.DataFrame, DG.DataFrame] {
  return [genChainDataFrame(ids, 'H'), genChainDataFrame(ids, 'L')];
}
