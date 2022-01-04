import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Searches for activity cliffs in a chemical dataset by selected cutoff*/
export async function getActivityCliffs(df: DG.DataFrame, smiles: DG.Column, activities: DG.Column, similarity: number) {
  function renderLines(sp: DG.Viewer, n1: number[], n2: number[]) {
    //@ts-ignore
    const ctx = sp.getInfo()['canvas'].getContext('2d');

    for (let i = 0; i < n1.length; i++) {
      ctx.beginPath();
      ctx.strokeStyle = 'green';
      ctx.lineWidth = 1;

      const num1 = n1[i];
      const num2 = n2[i];

      //@ts-ignore
      const pointFrom = sp.worldToScreen(sp.dataFrame.get('~x_coord', num1), sp.dataFrame.get('~y_coord', num1));
      ctx.lineTo(pointFrom.x, pointFrom.y);
      //@ts-ignore
      const pointTo = sp.worldToScreen(sp.dataFrame.get('~x_coord', num2), sp.dataFrame.get('~y_coord', num2));
      ctx.lineTo(pointTo.x, pointTo.y);

      ctx.stroke();
    }
  }

  const automaticSimilarityLimit = false;

  const AVERAGE_NEIGHBOR_COUNT = 6;
  const MIN_SIMILARITY = 80;
  const DEFAULT_SIMILARITY = 95;
  const VIEW_CYCLE_COUNT = 50000;

  let initialSimilarityLimit = 0;
  if (automaticSimilarityLimit) {
    initialSimilarityLimit = MIN_SIMILARITY;
  } else {
    initialSimilarityLimit = similarity / 100;
  }

  const colSmiles = DG.Column.fromList('string', 'smiles', smiles.toList());
  const dfSmiles = DG.DataFrame.fromColumns([colSmiles]);

  const arrSmiles = smiles.toList();
  let simArr: number[] = [];
  const dim = arrSmiles.length;

  for (let i = 0; i != dim - 1; ++i) {
    const mol = arrSmiles[i];

    const dfSmilesNew = dfSmiles.clone();
    dfSmilesNew.rows.removeAt(0, i + 1, false);
    const tmp = await grok.chem.getSimilarities(dfSmilesNew.col('smiles')!, mol);
    simArr = simArr.concat(tmp!.toList());
  }

  const optSimilarityLimit = initialSimilarityLimit;

  const values = activities.toList();
  const simVals = [];
  const diffVals = [];
  const saliVals = [];
  const n1: number[] = [];
  const n2: number[] = [];

  for (let i = 0; i != dim - 1; ++i) {
    for (let j = i + 1; j != dim; ++j) {
      const sim = simArr[i * dim + j];

      if (sim >= optSimilarityLimit) {
        n1.push(i);
        n2.push(j);
        simVals.push(sim);
        const diff = Math.abs(values[i] - values[j]);
        diffVals.push(diff);
        if (sim != 1) {
          saliVals.push(diff / (1 - sim));
        } else {
          saliVals.push(Infinity);
        }
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
      if (values[n1[i]] > values[n2[i]]) {
        saliCount[n1[i]] += saliVals[i];
      } else {
        saliCount[n2[i]] += saliVals[i];
      }
    }
  }

  const mx = new Array(dim).fill(0);
  const my = new Array(dim).fill(0);
  let mDx = new Array(dim).fill(0);
  let mDy = new Array(dim).fill(0);

  for (let i = 0; i != dim; ++i) {
    mx[i] = Math.random();
    my[i] = Math.random();
  }

  let sorted_index = Array.from(Array(dim).keys()).sort(function(a, b) {
    if (mx[a] < mx[b]) {
      return -1;
    }
    if (mx[a] > mx[b]) {
      return 1;
    }
    return 0;
  });

  const minDistance = 1 / Math.sqrt(dim);

  for (let cycle = 0; cycle != VIEW_CYCLE_COUNT; ++cycle) {
    mDx = new Array(dim).fill(0);
    mDy = new Array(dim).fill(0);

    const cycleState = cycle / VIEW_CYCLE_COUNT;
    const cycleMinDistance = cycleState * (2 - cycleState) * minDistance;
    let attractionCycleFactor = 0;
    let repulsionCycleFactor = 0;

    if (cycleState < 0.5) {
      attractionCycleFactor = 0.8 * (1 - cycleState) * 0.5;
      repulsionCycleFactor = 0.5;
    } else {
      attractionCycleFactor = 0.8 * (1 - cycleState) * (1 - cycleState);
      repulsionCycleFactor = 1 - cycleState;
    }

    for (let i = 0; i != n1.length; ++i) {
      const num1 = n1[i];
      const num2 = n2[i];
      let dx = mx[num2] - mx[num1];
      let dy = my[num2] - my[num1];
      const distance = Math.sqrt(dx * dx + dy * dy);
      const shift = distance - minDistance;
      let neigbourFactor1 = 0;
      let neigbourFactor2 = 0;

      if (shift > 0) {
        dx *= shift / distance;
        dy *= shift / distance;
        if (neighboursCount[num1] > 4) {
          neigbourFactor1 = 4.0 / neighboursCount[num1];
        } else {
          neigbourFactor1 = 1.0;
        }
        if (neighboursCount[num2] > 4) {
          neigbourFactor2 = 4.0 / neighboursCount[num2];
        } else {
          neigbourFactor2 = 1;
        }

        mDx[num1] = mDx[num1] + dx * attractionCycleFactor * neigbourFactor1;
        mDy[num1] = mDy[num1] + dy * attractionCycleFactor * neigbourFactor1;
        mDx[num2] = mDx[num2] - dx * attractionCycleFactor * neigbourFactor2;
        mDy[num2] = mDy[num2] - dy * attractionCycleFactor * neigbourFactor2;
      }
    }

    for (let i = 0; i != dim; ++i) {
      for (let j = i + 1; j != dim; ++j) {
        const num1 = sorted_index[i];
        const num2 = sorted_index[j];
        let dx = mx[num2] - mx[num1];

        if (Math.abs(dx) >= cycleMinDistance) {
          break;
        }
        let dy = my[num2] - my[num1];
        if (Math.abs(dy) < cycleMinDistance) {
          const distance = Math.sqrt(dx * dx + dy * dy);
          if (distance < cycleMinDistance) {
            if (distance == 0) {
              const angle = 2 * Math.PI * Math.random();
              dx = Math.sin(angle) * cycleMinDistance;
              dy = Math.cos(angle) * cycleMinDistance;
            } else {
              const shift = cycleMinDistance - distance;
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

  const vec_lengths = new Array(dim).fill(0);
  for (let i = 0; i != dim; ++i) {
    vec_lengths[i] = Math.sqrt((mx[i] - 0.5) * (mx[i] - 0.5) + (my[i] - 0.5) * (my[i] - 0.5));
  }


  sorted_index = Array.from(Array(dim).keys()).sort(function(a, b) {
    if (vec_lengths[a] < vec_lengths[b]) {
      return -1;
    }
    if (vec_lengths[a] > vec_lengths[b]) {
      return 1;
    }
    return 0;
  });

  for (let i = 0; i != dim; ++i) {
    const n = sorted_index[i];
    const x = mx[n] - 0.5;
    const y = my[n] - 0.5;
    let a = 0;

    if (y != 0) {
      a = Math.atan(x / y);
      if (y < 0) {
        if (x < 0) {
          a -= Math.PI;
        } else {
          a += Math.PI;
        }
      }
    } else {
      if (x > 0) {
        a = Math.PI / 2;
      } else {
        a = -Math.PI / 2;
      }
    }

    mx[n] = mx[n] * 2 - 1;
    my[n] = my[n] * 2 - 1;
  }


  const x_coord = DG.Column.fromList('double', 'x_coord', mx);
  const y_coord = DG.Column.fromList('double', 'y_coord', my);
  const sali = DG.Column.fromList('double', 'sali', saliCount);

  df.columns.insert(x_coord);
  df.columns.insert(y_coord);
  df.columns.insert(sali);

  df.columns.byName('x_coord').name = '~x_coord';
  df.columns.byName('y_coord').name = '~y_coord';

  const view = grok.shell.getTableView(df.name);
  const sp = view.addViewer(DG.Viewer.scatterPlot(df, {
    xColumnName: '~x_coord',
    yColumnName: '~y_coord',
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
