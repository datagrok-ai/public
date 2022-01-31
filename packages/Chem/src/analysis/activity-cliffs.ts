import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {chemSpace} from './chem-space'
import * as chemSearches from '../chem-searches';

/** Searches for activity cliffs in a chemical dataset by selected cutoff*/
export async function getActivityCliffs(df: DG.DataFrame, smiles: DG.Column, activities: DG.Column, similarity: number, methodName: string) {
  const automaticSimilarityLimit = false;
  const MIN_SIMILARITY = 80;

  let initialSimilarityLimit = 0;
  if (automaticSimilarityLimit)
    initialSimilarityLimit = MIN_SIMILARITY;
  else
    initialSimilarityLimit = similarity / 100;


  const colSmiles = DG.Column.fromList('string', 'smiles', smiles.toList());
  const dfSmiles = DG.DataFrame.fromColumns([colSmiles]);

  const arrSmiles = smiles.toList();
  const dim = arrSmiles.length;
  let simArr: DG.Column[] = Array(dim -1);//[(dim*dim - dim/2)];


  for (let i = 0; i != dim - 1; ++i) {
    const mol = arrSmiles[i];
    const dfSmilesNew = dfSmiles.clone();
    dfSmilesNew.rows.removeAt(0, i + 1, false);
   //const tmp = await grok.chem.getSimilarities(dfSmilesNew.col('smiles')!, mol);
 
    // if(i == 4){
    //   let a =5;
    // } 
    simArr[i] = (await chemSearches.chemGetSimilarities(dfSmilesNew.col('smiles')!, mol))!;
  }

  const optSimilarityLimit = initialSimilarityLimit;

  const simVals: number[] = [];
  const diffVals: number[] = [];
  const saliVals: number[] = [];
  const n1: number[] = [];
  const n2: number[] = [];

  for (let i = 0; i != dim - 1; ++i) {
    for (let j = 0; j != dim - 1 - i; ++j) {

      if(i == 4 && j == 0){
        let a =5;
      } 
  

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

  const sali = DG.Column.fromList('double', 'sali', saliCount);
  const coords = await chemSpace(smiles, methodName, "Tanimoto");

  for (const col of coords)
    df.columns.add(col);
  
  df.columns.add(sali);

  const view = grok.shell.getTableView(df.name);
  const sp = view.addViewer(DG.Viewer.scatterPlot(df, {
    xColumnName: 'Embed_X',
    yColumnName: 'Embed_Y',
    size: 'sali',
  }));

  sp.props.showXSelector = false;
  sp.props.showYSelector = false;
  sp.props.showSizeSelector = false;
  sp.props.showColorSelector = false;
  sp.props.markerMinSize = 5;
  sp.props.markerMaxSize = 25;
  sp.props.colorColumnName = 'activity';

  sp.onEvent('d4-before-draw-scene').subscribe((_) => renderLines(sp, n1, n2));
}

function renderLines(sp: DG.Viewer, n1: number[], n2: number[]) {
  //@ts-ignore
  const ctx = sp.getInfo()['canvas'].getContext('2d');

  const x = sp.dataFrame!.columns.byName('Embed_X');
  const y = sp.dataFrame!.columns.byName('Embed_Y');

  for (let i = 0; i < n1.length; i++) {
    ctx.beginPath();
    ctx.strokeStyle = 'green';
    ctx.lineWidth = 1;

    const num1 = n1[i];
    const num2 = n2[i];
    
    //@ts-ignore
    const pointFrom = sp.worldToScreen(x.get(num1), y.get(num1));
    ctx.lineTo(pointFrom.x, pointFrom.y);
    //@ts-ignore
    const pointTo = sp.worldToScreen(x.get(num2), sp.dataFrame.get('Embed_Y', num2));
    ctx.lineTo(pointTo.x, pointTo.y);

    ctx.stroke();
  }
}