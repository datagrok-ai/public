import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function getActivityCliffs(df: DG.DataFrame, smiles: DG.Column, activities: DG.Column, similarity: number) {

  function renderLines(sp: DG.Viewer, n1: number[], n2: number[]) {

    //@ts-ignore
    let ctx = sp.getInfo()['canvas'].getContext('2d');

    for (let i = 0; i < n1.length; i++) {
      ctx.beginPath();
      ctx.strokeStyle = 'green';
      ctx.lineWidth = 1;

      let num1 = n1[i];
      let num2 = n2[i];

      //@ts-ignore
      let pointFrom = sp.worldToScreen(sp.dataFrame.get('~x_coord', num1), sp.dataFrame.get('~y_coord', num1));
      ctx.lineTo(pointFrom.x, pointFrom.y);
      //@ts-ignore
      let pointTo = sp.worldToScreen(sp.dataFrame.get('~x_coord', num2), sp.dataFrame.get('~y_coord', num2));
      ctx.lineTo(pointTo.x, pointTo.y);

      ctx.stroke();
    }
  }

  let automaticSimilarityLimit = false;

  let AVERAGE_NEIGHBOR_COUNT = 6;
  let MIN_SIMILARITY = 80;
  let DEFAULT_SIMILARITY = 95;
  let VIEW_CYCLE_COUNT = 50000;

  let initialSimilarityLimit = 0;
  if (automaticSimilarityLimit) { initialSimilarityLimit = MIN_SIMILARITY; }
  else { initialSimilarityLimit = similarity / 100; }

  let colSmiles = DG.Column.fromList("string", "smiles", smiles.toList());
  let dfSmiles = DG.DataFrame.fromColumns([colSmiles]);

  let arrSmiles = smiles.toList();
  let simArr: number[] = [];
  let dim = arrSmiles.length;

  for (let i = 0; i != dim - 1; ++i) {
    let mol = arrSmiles[i];

    let dfSmilesNew = dfSmiles.clone();
    dfSmilesNew.rows.removeAt(0, i + 1, false);
    let tmp = await grok.chem.getSimilarities(dfSmilesNew.col('smiles')!, mol);
    simArr = simArr.concat(tmp!.toList());
  }

  let optSimilarityLimit = initialSimilarityLimit;

  let values = activities.toList();
  let simVals = [];
  let diffVals = [];
  let saliVals = [];
  let n1: number[] = [];
  let n2: number[] = [];

  for (let i = 0; i != dim - 1; ++i) {
    for (let j = i + 1; j != dim; ++j) {
      let sim = simArr[i * dim + j];

      if (sim >= optSimilarityLimit) {
        n1.push(i);
        n2.push(j);
        simVals.push(sim);
        let diff = Math.abs(values[i] - values[j]);
        diffVals.push(diff);
        if (sim != 1) { saliVals.push(diff / (1 - sim)); }
        else { saliVals.push(Infinity); }
      }
    }
  }

  let neighboursCount = new Array(dim).fill(0);
  let similarityCount = new Array(dim).fill(0);
  let saliCount = new Array(dim).fill(0);

  for (let i = 0; i != n1.length; ++i) {
    neighboursCount[n1[i]] += 1;
    neighboursCount[n2[i]] += 1;
    similarityCount[n1[i]] += simVals[i];
    similarityCount[n2[i]] += simVals[i];
    if (saliVals[i] != Infinity) {
      if (values[n1[i]] > values[n2[i]]) { saliCount[n1[i]] += saliVals[i]; }
      else { saliCount[n2[i]] += saliVals[i]; }
    }
  }

  let mx = new Array(dim).fill(0);
  let my = new Array(dim).fill(0);
  let mDx = new Array(dim).fill(0);
  let mDy = new Array(dim).fill(0);

  for (let i = 0; i != dim; ++i) {
    mx[i] = Math.random();
    my[i] = Math.random();
  }

  let sorted_index = Array.from(Array(dim).keys()).sort(function (a, b) {
    if (mx[a] < mx[b])
      return -1;
    if (mx[a] > mx[b])
      return 1;
    return 0;
  });

  let minDistance = 1 / Math.sqrt(dim);

  for (let cycle = 0; cycle != VIEW_CYCLE_COUNT; ++cycle) {
    mDx = new Array(dim).fill(0);
    mDy = new Array(dim).fill(0);

    let cycleState = cycle / VIEW_CYCLE_COUNT;
    let cycleMinDistance = cycleState * (2 - cycleState) * minDistance;
    let attractionCycleFactor = 0;
    let repulsionCycleFactor = 0;

    if (cycleState < 0.5) {
      attractionCycleFactor = 0.8 * (1 - cycleState) * 0.5
      repulsionCycleFactor = 0.5
    } else {
      attractionCycleFactor = 0.8 * (1 - cycleState) * (1 - cycleState)
      repulsionCycleFactor = 1 - cycleState
    }

    for (let i = 0; i != n1.length; ++i) {
      let num1 = n1[i];
      let num2 = n2[i];
      let dx = mx[num2] - mx[num1];
      let dy = my[num2] - my[num1];
      let distance = Math.sqrt(dx * dx + dy * dy);
      let shift = distance - minDistance;
      let neigbourFactor1 = 0;
      let neigbourFactor2 = 0;

      if (shift > 0) {
        dx *= shift / distance;
        dy *= shift / distance;
        if (neighboursCount[num1] > 4) { neigbourFactor1 = 4.0 / neighboursCount[num1]; }
        else { neigbourFactor1 = 1.0; }
        if (neighboursCount[num2] > 4) { neigbourFactor2 = 4.0 / neighboursCount[num2]; }
        else { neigbourFactor2 = 1 }

        mDx[num1] = mDx[num1] + dx * attractionCycleFactor * neigbourFactor1;
        mDy[num1] = mDy[num1] + dy * attractionCycleFactor * neigbourFactor1;
        mDx[num2] = mDx[num2] - dx * attractionCycleFactor * neigbourFactor2;
        mDy[num2] = mDy[num2] - dy * attractionCycleFactor * neigbourFactor2;
      }
    }

    for (let i = 0; i != dim; ++i) {
      for (let j = i + 1; j != dim; ++j) {
        let num1 = sorted_index[i];
        let num2 = sorted_index[j];
        let dx = mx[num2] - mx[num1];

        if (Math.abs(dx) >= cycleMinDistance) { break; }
        let dy = my[num2] - my[num1];
        if (Math.abs(dy) < cycleMinDistance) {

          let distance = Math.sqrt(dx * dx + dy * dy);
          if (distance < cycleMinDistance) {

            if (distance == 0) {
              let angle = 2 * Math.PI * Math.random();
              dx = Math.sin(angle) * cycleMinDistance;
              dy = Math.cos(angle) * cycleMinDistance;
            } else {
              let shift = cycleMinDistance - distance;
              dx *= shift / distance;
              dy *= shift / distance;
            }

            mDx[num1] = mDx[num1] - dx * repulsionCycleFactor;
            mDy[num1] = mDy[num1] - dy * repulsionCycleFactor;
            mDx[num2] = mDx[num2] + dx * repulsionCycleFactor;
            mDy[num2] = mDy[num2] + dy * repulsionCycleFactor;
          }
        }
      }
    }

    for (let i = 0; i != dim; ++i) {
      mx[i] = Math.min(Math.max(mx[i] + mDx[i], 0), 1);
      my[i] = Math.min(Math.max(my[i] + mDy[i], 0), 1);
    }
  }

  let vec_lengths = new Array(dim).fill(0);
  for (let i = 0; i != dim; ++i) {
    vec_lengths[i] = Math.sqrt((mx[i] - 0.5) * (mx[i] - 0.5) + (my[i] - 0.5) * (my[i] - 0.5));
  }


  sorted_index = Array.from(Array(dim).keys()).sort(function (a, b) {
    if (vec_lengths[a] < vec_lengths[b])
      return -1;
    if (vec_lengths[a] > vec_lengths[b])
      return 1;
    return 0;
  });

  for (let i = 0; i != dim; ++i) {
    let n = sorted_index[i];
    let x = mx[n] - 0.5;
    let y = my[n] - 0.5;
    let a = 0;

    if (y != 0) {
      a = Math.atan(x / y)
      if (y < 0) {
        if (x < 0)
          a -= Math.PI;
        else
          a += Math.PI;
      }
    } else {
      if (x > 0)
        a = Math.PI / 2;
      else
        a = -Math.PI / 2;
    }

    mx[n] = mx[n] * 2 - 1;
    my[n] = my[n] * 2 - 1;
  }


  let x_coord = DG.Column.fromList("double", "x_coord", mx);
  let y_coord = DG.Column.fromList("double", "y_coord", my);
  let sali = DG.Column.fromList("double", "sali", saliCount);

  df.columns.insert(x_coord);
  df.columns.insert(y_coord);
  df.columns.insert(sali);

  df.columns.byName('x_coord').name = '~x_coord';
  df.columns.byName('y_coord').name = '~y_coord';

  let view = grok.shell.getTableView(df.name);
  let sp = view.addViewer(DG.Viewer.scatterPlot(df, {
    xColumnName: '~x_coord',
    yColumnName: '~y_coord',
    size: 'sali'
  }));

  sp.props.showXSelector = false;
  sp.props.showYSelector = false;
  sp.props.showSizeSelector = false;
  sp.props.showColorSelector = false;
  sp.props.markerMinSize = 5;
  sp.props.markerMaxSize = 25;
  sp.props.colorColumnName = 'activity';

  sp.onEvent('d4-before-draw-scene').subscribe(_ => renderLines(sp, n1, n2));
}  