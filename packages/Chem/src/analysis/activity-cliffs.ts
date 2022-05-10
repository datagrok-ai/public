import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {chemSpace} from './chem-space';
import * as chemSearches from '../chem-searches';
import {getSimilarityFromDistance} from '@datagrok-libraries/utils/src/similarity-metrics';
import {ScatterPlotViewer} from "datagrok-api/dg";

const options = {
  'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005}
}

// Searches for activity cliffs in a chemical dataset by selected cutoff
export async function getActivityCliffs(
  df: DG.DataFrame,
  smiles: DG.Column,
  activities: DG.Column,
  similarity: number,
  methodName: string) {
  const automaticSimilarityLimit = false;
  const MIN_SIMILARITY = 80;

  let initialSimilarityLimit = automaticSimilarityLimit ? MIN_SIMILARITY : similarity / 100;
  const axes = ['Embed_X', 'Embed_Y'];
  //@ts-ignore
  const colNameInd = df.columns.names().filter((it) => it.includes(axes[0])).length + 1;

  const {distance, coordinates} = await chemSpace(smiles, methodName, 'Tanimoto',
    axes.map((it) => `${it}_${colNameInd}`), (options as any)[methodName], true);

  for (const col of coordinates)
    df.columns.add(col);

  const dfSmiles = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'smiles', smiles.toList())]);
  const dim = smiles.length;
  const simArr: DG.Column[] = Array(dim - 1);

  if (!distance)
    await getSimilaritiesMarix(dim, smiles, dfSmiles, simArr);
  else
    getSimilaritiesMarixFromDistances(dim, distance, simArr);

  const optSimilarityLimit = initialSimilarityLimit;

  const simVals: number[] = [];
  const diffVals: number[] = [];
  const saliVals: number[] = [];
  const n1: number[] = [];
  const n2: number[] = [];

  for (let i = 0; i != dim - 1; ++i) {
    for (let j = 0; j != dim - 1 - i; ++j) {
      const sim: number = simArr[i].get(j);

      if (sim >= optSimilarityLimit) {
        n1.push(i);
        n2.push(i + j + 1);
        simVals.push(sim);
        const diff = Math.abs(activities.get(i) - activities.get(i + j + 1));
        diffVals.push(diff);
        if (sim != 1)
          saliVals.push(diff / (1 - sim));
        else
          saliVals.push(Infinity);
      }
    }
  }

  const neighboursCount = new Array(dim).fill(0);
  const similarityCount = new Array(dim).fill(0);
  const saliCount = new Array(dim).fill(0);

  for (let i = 0; i != n1.length; ++i) {
    neighboursCount[n1[i]] += 1;
    neighboursCount[n2[i]] += 1;
    similarityCount[n1[i]] += simVals[i];
    similarityCount[n2[i]] += simVals[i];
    if (saliVals[i] != Infinity) {
      if (activities.get(n1[i]) > activities.get(n2[i]))
        saliCount[n1[i]] += saliVals[i];
      else
        saliCount[n2[i]] += saliVals[i];
    }
  }

  const sali = DG.Column.fromList('double', `sali_${colNameInd}`, saliCount);

  df.columns.add(sali);

  const view = grok.shell.getTableView(df.name);
  const sp = view.addViewer(DG.Viewer.scatterPlot(df, {
    xColumnName: `${axes[0]}_${colNameInd}`,
    yColumnName: `${axes[1]}_${colNameInd}`,
    size: sali.name,
    color: activities.name,
    showXSelector: false,
    showYSelector: false,
    showSizeSelector: false,
    showColorSelector: false,
    markerMinSize: 5,
    markerMaxSize: 25,
  })) as ScatterPlotViewer;

  sp.onEvent('d4-before-draw-scene')
    .subscribe((_: any) => renderLines(sp, n1, n2, `${axes[0]}_${colNameInd}`, `${axes[1]}_${colNameInd}`));

  sp.addProperty('similarityLimit', 'double', optSimilarityLimit);
}

function renderLines(sp: DG.ScatterPlotViewer, n1: number[], n2: number[], xAxis: string, yAxis: string) {
  const ctx = sp.getInfo()['canvas'].getContext('2d');
  const x = sp.dataFrame!.columns.byName(xAxis);
  const y = sp.dataFrame!.columns.byName(yAxis);

  for (let i = 0; i < n1.length; i++) {
    ctx.beginPath();
    ctx.strokeStyle = 'green';
    ctx.lineWidth = 1;

    const num1 = n1[i];
    const num2 = n2[i];

    const pointFrom = sp.worldToScreen(x.get(num1), y.get(num1));
    ctx.lineTo(pointFrom.x, pointFrom.y);
    const pointTo = sp.worldToScreen(x.get(num2), sp.dataFrame.get(yAxis, num2));
    ctx.lineTo(pointTo.x, pointTo.y);

    ctx.stroke();
  }
}

async function getSimilaritiesMarix(dim: number, smiles: DG.Column, dfSmiles: DG.DataFrame, simArr: DG.Column[]) {
  for (let i = 0; i != dim - 1; ++i) {
    const mol = smiles.get(i);
    dfSmiles.rows.removeAt(0, 1, false);
    simArr[i] = (await chemSearches.chemGetSimilarities(dfSmiles.col('smiles')!, mol))!;
  }
  return simArr;
}

function getSimilaritiesMarixFromDistances(dim: number, distances: [], simArr: DG.Column[]) {
  for (let i = 0; i < dim - 1; ++i) {
    const similarityArr = [];
    for (let j = i + 1; j < dim; ++j)
      similarityArr.push(getSimilarityFromDistance(distances[i][j]));
    simArr[i] = DG.Column.fromFloat32Array('similarity', Float32Array.from(similarityArr));
  }
  return simArr;
}
