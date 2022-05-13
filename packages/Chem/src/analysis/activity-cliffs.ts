import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {chemSpace} from './chem-space';
import * as chemSearches from '../chem-searches';
import {getSimilarityFromDistance} from '@datagrok-libraries/utils/src/similarity-metrics';
import {renderMolecule} from '../rendering/render-molecule';
import {updateDivInnerHTML} from '../utils/ui-utils';
import {findMCS} from '../scripts-api';

const options = {
  'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005},
};

// Searches for activity cliffs in a chemical dataset by selected cutoff
export async function getActivityCliffs(
  df: DG.DataFrame,
  smiles: DG.Column,
  activities: DG.Column,
  similarity: number,
  methodName: string) {
  const automaticSimilarityLimit = false;
  const MIN_SIMILARITY = 80;

  const initialSimilarityLimit = automaticSimilarityLimit ? MIN_SIMILARITY : similarity / 100;
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
  }));
  const div = ui.div('');
  view.root.append(div);

  //@ts-ignore
  const canvas = sp.getInfo()['canvas'];
  let lines: any[] = [];
  //@ts-ignore
  canvas.addEventListener('mousedown', function(event) {
    setTimeout(() => {
      const rect = canvas.getBoundingClientRect();
      const x = event.clientX - rect.left;
      const y = event.clientY - rect.top;
      let clickedOnLine = false;
      for (const line of lines) {
        const dist =
      Math.abs(Math.hypot(line.a[0] - x, line.a[1] - y) +
      Math.hypot(line.b[0] - x, line.b[1] - y) - Math.hypot(line.a[0] - line.b[0], line.a[1] - line.b[1]));
        if (dist < 2) {
          const tooltipDiv = ui.divV([]);
          //@ts-ignore
          line.mols.forEach((mol) => {
            tooltipDiv.append(renderMolecule(
            //@ts-ignore
              df.get('smiles', mol), {width: 120, height: 60}));
          });
          tooltipDiv.append(ui.divText('MCS:'));
          const mscDiv = ui.div('Loading mcs...');
          tooltipDiv.append(mscDiv);
          const mcsDf = DG.DataFrame.create(2);
          //@ts-ignore
          mcsDf.columns.addNewString('smiles').init((i) => df.get('smiles', line.mols[i]));
          findMCS('smiles', mcsDf).then((mcs) => {
            updateDivInnerHTML(mscDiv, mcs);
          });
          updateDivInnerHTML(div, tooltipDiv);
          div.style.top = (event.clientY) + 'px';
          div.style.left = (event.clientX) + 'px';
          div.style.position = 'fixed';
          div.style.display = 'block';
          div.style.background = 'snow';
          div.style.border = 'solid 1px';
          clickedOnLine = true;
          break;
        }
      }
      if (!clickedOnLine) {
        div.style.border = 'none';
        div.innerHTML = '';
      }

      event.stopPropagation();
    }, 500);
  });

  sp.onEvent('d4-before-draw-scene')
    .subscribe((_) => {
      div.innerHTML = '';
      div.style.border = 'none';
      lines = renderLines(sp, n1, n2, `${axes[0]}_${colNameInd}`, `${axes[1]}_${colNameInd}`);
    });

  sp.addProperty('similarityLimit', 'double', optSimilarityLimit);
}

function renderLines(sp: DG.Viewer, n1: number[], n2: number[], xAxis: string, yAxis: string) {
  //@ts-ignore
  const canvas = sp.getInfo()['canvas'];
  const ctx = canvas.getContext('2d') as CanvasRenderingContext2D;

  const x = sp.dataFrame!.columns.byName(xAxis);
  const y = sp.dataFrame!.columns.byName(yAxis);

  const lines = [];

  for (let i = 0; i < n1.length; i++) {
    const line = new Path2D();
    const num1 = n1[i];
    const num2 = n2[i];
    //@ts-ignore
    const pointFrom = sp.worldToScreen(x.get(num1), y.get(num1));
    //@ts-ignore
    const pointTo = sp.worldToScreen(x.get(num2), sp.dataFrame.get(yAxis, num2));
    line.moveTo(pointFrom.x, pointFrom.y);
    ctx.strokeStyle = 'green';
    ctx.lineWidth = 1;
    line.lineTo(pointTo.x, pointTo.y);
    ctx.stroke(line);
    lines.push({a: [pointFrom.x, pointFrom.y], b: [pointTo.x, pointTo.y], mols: [num1, num2]});
  }
  return lines;
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
