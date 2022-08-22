import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { expect } from '@datagrok-libraries/utils/src/test';
import { chemSpace } from '../analysis/chem-space';
import { getSimilaritiesMarix, getSimilaritiesMarixFromDistances } from '../utils/similarity-utils';
import { chemSpaceTopMenu } from '../package';
var { jStat } = require('jstat')


export async function _testChemSpaceReturnsResult(df: DG.DataFrame, algorithm: string) {
  const v = grok.shell.addTableView(df);
  const sp = await chemSpaceTopMenu(df, df.col('smiles')!, algorithm, 'Tanimoto', true);
  expect(sp != null, true);
  v.close();
}

export async function _testDimensionalityReducer(col: DG.Column, algorithm: string) {
  const chemSpaceParams = {
    seqCol: col,
    methodName: algorithm,
    similarityMetric: 'Tanimoto',
    embedAxesNames: ['Embed_X', 'Embed_Y']
  }
  const { distance, coordinates } = await chemSpace(chemSpaceParams);

  const dfSmiles = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'smiles', col.toList())]);
  const dim = col.length;
  const simArr: DG.Column[] = Array(dim - 1);

  if (!distance)
    await getSimilaritiesMarix(dim, col, dfSmiles, simArr);
  else
    getSimilaritiesMarixFromDistances(dim, distance, simArr);

  const nearestAndFarestNeighbours: IDistanceToPoint[][] = findNNearestAndFarestNeighbours(coordinates, col.length, 'Embed_X', 'Embed_Y');
  // const similaririesWithDistances = DG.DataFrame.create();
  const similaritiesArray: number[] = [];
  const distancesArray: number[] = [];
  /*     similaririesWithDistances.columns.addNewInt('mol1');
      similaririesWithDistances.columns.addNewInt('mol2');
      similaririesWithDistances.columns.addNewFloat('distance');
      similaririesWithDistances.columns.addNewFloat('similarity'); */
  nearestAndFarestNeighbours.forEach((molecule, index) => {
    molecule.forEach((molPair) => {
      const arrayInd = index > molPair.idx ? molPair.idx : index;
      const row = index < molPair.idx ? molPair.idx : index;
      const sim = simArr[arrayInd].get(row - arrayInd - 1);
      //   similaririesWithDistances.rows.addNew([index, molPair.idx, molPair.distance, sim]);
      similaritiesArray.push(sim);
      distancesArray.push(molPair.distance);
    })
  })
  //  grok.shell.addTableView(similaririesWithDistances);
  const corrCoef = jStat.corrcoeff(similaritiesArray, distancesArray);
  expect(corrCoef < -0.75, true);
}

interface IDistanceToPoint {
  idx: number;
  distance: number;
}

function findNNearestAndFarestNeighbours(coordinates: DG.ColumnList, nItems: number, xColName: string, yColName: string, n?: number) {

  const matrix: IDistanceToPoint[][] = [];
  const df = DG.DataFrame.create(nItems);
  for (const col of coordinates)
    df.columns.add(col);

  for (let i = 0; i < nItems; ++i) {
    const distances: IDistanceToPoint[] = [];
    for (let j = i + 1; j < nItems; ++j) {
      const x1 = df.col(xColName)!.get(i);
      const y1 = df.col(yColName)!.get(i);
      const x2 = df.col(xColName)!.get(j);
      const y2 = df.col(yColName)!.get(j);
      const d: number = euclideanDistance(x1, y1, x2, y2);
      /*       if (d === NaN) {
                console.log(x1, '***',y1, '***',x2, '***',y2)
            } */
      distances.push({ idx: j, distance: d });
    }
    if (n) {
      distances.sort((a, b) => a.distance - b.distance);
      matrix.push(distances.slice(0, n).concat(distances.slice(-5)));
    } else {
      matrix.push(distances);
    }
  }
  return matrix;
}


function euclideanDistance(x1: number, y1: number, x2: number, y2: number) {
  return Math.sqrt(Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2));
}
