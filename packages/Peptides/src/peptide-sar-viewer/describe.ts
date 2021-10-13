import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore
import * as jStat from 'jstat';

import {splitAlignedPeptides} from '../split-aligned';
import {decimalAdjust, tTest} from '../utils/misc';
import {ChemPalette} from '../utils/chem-palette';

const cp = new ChemPalette('grok');

export async function describe(
  df: DG.DataFrame,
  activityColumn: string,
  activityScaling: string,
  filterMode: boolean,
  sourceGrid: DG.Grid,
): Promise<[DG.Grid, DG.DataFrame] | [null, null]> {
  //Split the aligned sequence into separate AARs
  let splitSeqDf: DG.DataFrame | undefined;
  const col: DG.Column = df.columns.bySemType('alignedSequence');
  if (col) {
    splitSeqDf = splitAlignedPeptides(col);
    splitSeqDf.name = 'Split sequence';
  }

  if (typeof splitSeqDf === 'undefined') {
    return [null, null];
  }

  const positionColumns = splitSeqDf.columns.names();
  const activityColumnScaled = `${activityColumn}Scaled`;

  splitSeqDf.columns.add(df.getCol(activityColumn));

  if (df.col(activityColumnScaled)) {
    df.columns.remove(activityColumnScaled);
  }

  //FIXME: this column usually duplicates, so remove it then
  if (df.col(`${activityColumnScaled} (2)`)) {
    df.columns.remove(`${activityColumnScaled} (2)`);
  }

  // append splitSeqDf columns to source table and make sure columns are not added more than once
  const dfColsSet = new Set(df.columns.names());
  if (!positionColumns.every((col: string) => dfColsSet.has(col))) {
    df.join(splitSeqDf, [activityColumn], [activityColumn], df.columns.names(), positionColumns, 'inner', true);
  }

  for (const col of df.columns) {
    if (splitSeqDf.col(col.name) && col.name != activityColumn) {
      col.semType = 'aminoAcids';
      col.setTag('cell.renderer', 'aminoAcids');
    }
  }

  // scale activity
  switch (activityScaling) {
  case 'lg':
    await df.columns.addNewCalculated(activityColumnScaled, 'Log10(${' + activityColumn + '})');
    splitSeqDf.columns.add(df.getCol(activityColumnScaled));
    sourceGrid.col(activityColumnScaled)!.name = `Log10(${activityColumn})`;
    sourceGrid.columns.setOrder([`Log10(${activityColumn})`]);
    break;
  case '-lg':
    await df.columns.addNewCalculated(activityColumnScaled, '-1*Log10(${' + activityColumn + '})');
    splitSeqDf.columns.add(df.getCol(activityColumnScaled));
    sourceGrid.col(activityColumnScaled)!.name = `-Log10(${activityColumn})`;
    sourceGrid.columns.setOrder([`-Log10(${activityColumn})`]);
    break;
  default:
    await df.columns.addNewCalculated(activityColumnScaled, '${' + activityColumn + '}');
    splitSeqDf.columns.add(df.getCol(activityColumnScaled));
    sourceGrid.col(activityColumnScaled)!.name = `${activityColumn}`;
    sourceGrid.columns.setOrder([`${activityColumn}`]);
    break;
  }

  const positionColName = 'position';
  const aminoAcidResidue = 'aminoAcidResidue';
  // const medianColName = 'MAD';

  //unpivot a table and handle duplicates
  splitSeqDf = splitSeqDf.groupBy(positionColumns)
    .add('med', activityColumnScaled, activityColumnScaled)
    .aggregate();

  const peptidesCount = splitSeqDf.getCol(activityColumnScaled).length;

  let matrixDf = splitSeqDf.unpivot([activityColumnScaled], positionColumns, positionColName, aminoAcidResidue);

  //this table contains overall statistics on activity
  const totalStats = matrixDf.groupBy()
    .add('med', activityColumnScaled, 'med')
    .aggregate();

  //preparing for mad
  const formula = 'Abs(${' + activityColumnScaled + '}-' + totalStats.get('med', 0) + ')';
  await matrixDf.columns.addNewCalculated('innerMAD', formula);

  //statistics for specific AAR at a specific position
  matrixDf = matrixDf.groupBy([positionColName, aminoAcidResidue])
    .add('count', activityColumnScaled, 'Count')
    //.add('med', 'innerMAD', medianColName) //final step of MAD calculation
    //.add('q1', activityColumnScaled, 'Q1')
    //.add('q2', activityColumnScaled, 'Median')
    //.add('q3', activityColumnScaled, 'Q3')
    .aggregate();
  
  const countThreshold = 4;
  //@ts-ignore: never gets old
  matrixDf.rows.filter((row) => row.Count >= countThreshold && row.Count <= peptidesCount - countThreshold);
  matrixDf = matrixDf.clone(matrixDf.filter);

  // calculate additional stats
  //await matrixDf.columns.addNewCalculated('IQR', '${q3}-${q1}');
  //await matrixDf.columns.addNewCalculated('CQV', '(${q3}-${q1})/(${q3}+${q1})');
  await matrixDf.columns.addNewCalculated('Ratio', '${count}/'.concat(`${peptidesCount}`));

  //calculate p-values based on t-test
  const pValues: number[] = [];
  const mDiff: number[] = [];
  let position: string;
  let AAR: string;
  let currentActivity: number[];
  let otherActivity: number[];
  let testResult;
  let currentMeanDiff: number;

  for (let i = 0; i < matrixDf.rowCount; i++) {
    position = matrixDf.get(positionColName, i);
    AAR = matrixDf.get(aminoAcidResidue, i);

    //@ts-ignore
    splitSeqDf.rows.select((row) => row[position] === AAR);
    currentActivity = splitSeqDf
      .clone(splitSeqDf.selection, [activityColumnScaled])
      .getCol(activityColumnScaled)
      .toList();

    //@ts-ignore
    splitSeqDf.rows.select((row) => row[position] !== AAR);
    otherActivity = splitSeqDf
      .clone(splitSeqDf.selection, [activityColumnScaled])
      .getCol(activityColumnScaled)
      .toList();

    testResult = tTest(currentActivity, otherActivity);
    // testResult = uTest(currentActivity, otherActivity);
    currentMeanDiff = testResult['Mean difference']!;

    mDiff.push(currentMeanDiff);
    pValues.push(testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less']);
  }
  matrixDf.columns.add(DG.Column.fromList(DG.TYPE.FLOAT, 'Mean difference', mDiff));
  matrixDf.columns.add(DG.Column.fromList(DG.TYPE.FLOAT, 'p-value', pValues));


  const statsDf = matrixDf.clone();

  //pivot a table to make it matrix-like
  matrixDf = matrixDf.groupBy([aminoAcidResidue])
    .pivot(positionColName)
    .add('first', 'Mean difference', '')
    .aggregate();
  matrixDf.name = 'SAR';

  // Setting category order
  // let aarCategoryOrder: string[] = ['-'];
  // for (const group of ChemPalette.grokGroups) {
  //   aarCategoryOrder = aarCategoryOrder.concat(group[0]);
  // }
  // matrixDf.getCol(aminoAcidResidue).setCategoryOrder(aarCategoryOrder);
  let aarList: string[] = statsDf.getCol(aminoAcidResidue).toList();
  const aarWeights: {[index: string]: number} = {};
  aarList.forEach((value) => {
    aarWeights[value] = (aarWeights[value] || 0) + 1;
  });
  aarList = aarList.filter((value, index, self) => self.indexOf(value) === index);
  aarList.sort((first, second) => aarWeights[second] - aarWeights[first]);
  matrixDf.getCol(aminoAcidResidue).setCategoryOrder(aarList);

  // !!! DRAWING PHASE !!!
  //find min and max MAD across all of the dataframe
  // const dfMin = jStat.min(mDiff);
  const dfMax = jStat.max(mDiff);
  const grid = matrixDf.plot.grid();

  grid.sort([aminoAcidResidue]);
  grid.columns.setOrder([aminoAcidResidue].concat(positionColumns));

  for (const col of matrixDf.columns) {
    if (col.name === aminoAcidResidue) {
      col.semType = 'aminoAcids';
      col.setTag('cell.renderer', 'aminoAcids');
      // let maxLen = 0;
      // col.categories.forEach( (ent:string)=>{
      //   if ( ent.length > maxLen) {
      //     maxLen = ent.length;
      //   }
      // });
      // grid.columns.byName(aminoAcidResidue)!.width = maxLen * 15;
    }
  }

  //render column headers and AAR symbols centered
  grid.onCellRender.subscribe(function(args: DG.GridCellRenderArgs) {
    args.g.save();
    args.g.beginPath();
    args.g.rect(args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height);
    args.g.clip();

    if (args.cell.isRowHeader && args.cell.gridColumn.visible) {
      args.cell.gridColumn.visible = false;
      args.preventDefault();
      return;
    }

    if (args.cell.isColHeader) {
      if (args.cell.gridColumn.name != aminoAcidResidue) {
        const textSize = args.g.measureText(args.cell.gridColumn.name);
        args.g.fillStyle = '#4b4b4a';
        args.g.fillText(
          args.cell.gridColumn.name,
          args.bounds.x + (args.bounds.width - textSize.width) / 2,
          args.bounds.y + (textSize.actualBoundingBoxAscent + textSize.actualBoundingBoxDescent),
        );
      }
      args.preventDefault();
    }

    if (args.cell.isTableCell && args.cell.tableRowIndex !== null && args.cell.tableColumn !== null) {
      if (args.cell.cell.value !== null && args.cell.tableColumn.name !== aminoAcidResidue) {
        const query =
          `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, args.cell.tableRowIndex)} ` +
          `and ${positionColName} = ${args.cell.tableColumn.name}`;

        //don't draw AAR that too little appearnces at this position
        const count: number = statsDf.groupBy(['Count']).where(query).aggregate().get('Count', 0);
        // if (count < countThreshold || count > peptidesCount - countThreshold) {
        //   args.preventDefault();
        //   args.g.restore();
        //   return;
        // }

        const pVal: number = statsDf.groupBy(['p-value']).where(query).aggregate().get('p-value', 0);

        let coef;
        const variant = args.cell.cell.value >= 0;
        if (pVal < 0.01) {
          coef = variant ? '#299617' : '#722F37';
        } else if (pVal < 0.05) {
          coef = variant ? '#32CD32' : '#D10000';
        } else if (pVal < 0.1) {
          coef = variant ? '#98FF98' : '#FF8A8A';
        } else {
          coef = DG.Color.toHtml(DG.Color.lightLightGray);
        }

        const rCoef = Math.abs(args.cell.cell.value) / dfMax;

        const maxRadius = 0.9 * (args.bounds.width > args.bounds.height ? args.bounds.height : args.bounds.width) / 2;
        const radius = Math.ceil(maxRadius * rCoef);

        args.g.beginPath();
        // args.g.fillStyle = DG.Color.toHtml(DG.Color.scaleColor(
        //   coef,
        //   0,
        //   1,
        //   undefined,
        //   [DG.Color.darkGray, args.cell.cell.value >= 0 ? DG.Color.green : DG.Color.red],
        // ));
        args.g.fillStyle = coef;
        args.g.arc(
          args.bounds.x + args.bounds.width / 2,
          args.bounds.y + args.bounds.height / 2,
          radius < 3 ? 3 : radius,
          0,
          Math.PI * 2,
          true,
        );
        args.g.closePath();

        args.g.fill();
        args.preventDefault();
      }
    }
    args.g.restore();
  });

  // show all the statistics in a tooltip over cell
  grid.onCellTooltip(function(cell, x, y) {
    if (
      !cell.isRowHeader &&
        !cell.isColHeader &&
        cell.tableColumn !== null &&
        cell.tableColumn.name !== aminoAcidResidue &&
        cell.cell.value !== null &&
        cell.tableRowIndex !== null
    ) {
      const tooltipMap: { [index: string]: string } = {};

      for (const col of statsDf.columns.names()) {
        if (col !== aminoAcidResidue && col !== positionColName) {
          const query =
              `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, cell.tableRowIndex)} ` +
              `and ${positionColName} = ${cell.tableColumn.name}`;
          let text = `${decimalAdjust('floor', statsDf.groupBy([col]).where(query).aggregate().get(col, 0), -5)}`;

          //@ts-ignore: I'm sure it's gonna be fine, text contains a number
          if (col === 'Count') {
            // if (parseInt(text) < countThreshold || parseInt(text) > peptidesCount - countThreshold) {
            //   return true;
            // }
            text += ` / ${peptidesCount}`;
          } else if (col === 'p-value') {
            text = parseFloat(text) !== 0 ? text : '<0.01'; 
          }
          
          tooltipMap[col] = text;
        }
      }

      ui.tooltip.show(ui.tableFromMap(tooltipMap), x, y);
      // } else if (cell.isColHeader && !cell.isRowHeader) {
      //   ui.tooltip.show((await df.plot.fromType('peptide-logo-viewer')).root, x, y);
    }
    if (
      !cell.isColHeader &&
        cell.tableColumn !== null &&
        cell.tableColumn.name == aminoAcidResidue &&
        cell.cell.value !== null &&
        cell.tableRowIndex !== null
    ) {
      cp.showTooltip(cell, x, y);
    }
    return true;
  });

  // Select columns in source table that correspond to the currently clicked cell
  grid.table.onCurrentCellChanged.subscribe((_: any) => {
    if (grid.table.currentCell.value && grid.table.currentCol.name !== aminoAcidResidue) {
      const currentAAR: string = grid.table.get(aminoAcidResidue, grid.table.currentRowIdx);
      const currentPosition = grid.table.currentCol.name;

      const query = `aminoAcidResidue = ${currentAAR} and position = ${currentPosition}`;
      // const text = statsDf.groupBy(['Count']).where(query).aggregate().get('Count', 0);
      // if (text < countThreshold || text > peptidesCount - countThreshold) {
      //   return;
      // }

      const splitColName = '~splitCol';
      const otherLabel = 'Other';
      const aarLabel = `${currentAAR === '-' ? 'Empty' : currentAAR} - ${currentPosition}`;

      const bitset = filterMode ? df.filter : df.selection;

      if (!df.col(splitColName)) {
        df.columns.addNew(splitColName, 'string');
      }

      let isChosen: boolean;
      for (let i = 0; i < df.rowCount; i++) {
        isChosen = df.get(currentPosition, i) === currentAAR;
        bitset.set(i, isChosen);
        df.getCol(splitColName).set(i, isChosen ? aarLabel : otherLabel);
      }

      // bitset.init((i) => df.get(currentPosition, i) === currentAAR);

      // const splitArray: string[] = [];
      // for (let i = 0; i < bitset.length; i++) {
      //   splitArray.push(bitset.get(i) ? aarLabel : otherLabel);
      // }

      // //FIXME: it resets filter; don't replace, init instead
      // const splitCol = DG.Column.fromStrings(splitColName, splitArray);
      // if (!df.col(splitColName)) {
      //   df.columns.add(splitCol);
      // } else {
      //   df.columns.replace(splitColName, splitCol);
      // }

      // df.getCol(splitColName).setCategoryOrder([otherLabel, aarLabel]);
      const colorMap: {[index: string]: string | number} = {};
      colorMap[otherLabel] = DG.Color.blue;
      colorMap[aarLabel] = DG.Color.orange;
      // colorMap[currentAAR] = cp.getColor(currentAAR);
      df.getCol(splitColName).colors.setCategorical(colorMap);
    }
  });

  for (const col of matrixDf.columns.names()) {
    console.log(grid.props['rowHeight']);
    grid.col(col)!.width = grid.props['rowHeight']; 
  }

  return [grid, statsDf];
}
