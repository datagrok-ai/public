import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore
import * as jStat from 'jstat';

import {splitAlignedPeptides} from '../split-aligned';
import {decimalAdjust, tTest} from '../utils/misc';
import {ChemPalette} from '../utils/chem-palette';
import * as grok from 'datagrok-api/src/chem';

const cp = new ChemPalette('grok');

export async function describe(
  df: DG.DataFrame,
  activityColumn: string,
  activityScaling: string,
  filterMode: boolean,
): Promise<[DG.Grid, DG.DataFrame] | [null, null]> {
  //Split the aligned sequence into separate AARs
  let splitSeqDf: DG.DataFrame | undefined;
  for (const col of df.columns) {
    if (col.semType === 'alignedSequence') {
      splitSeqDf = splitAlignedPeptides(col);
      splitSeqDf.name = 'Split sequence';
      break;
    }
  }

  if (typeof splitSeqDf === 'undefined') {
    return [null, null];
  }

  const positionColumns = splitSeqDf.columns.names();
  const activityColumnScaled = `~${activityColumn}Scaled`;

  splitSeqDf.columns.add(df.getCol(activityColumn));

  if (df.col(activityColumnScaled)) {
    df.columns.remove(activityColumnScaled);
  }

  //FIXME: this column usually duplicates, so remove it then
  if (df.col('~IC50Scaled (2)')) {
    df.columns.remove('~IC50Scaled (2)');
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
    break;
  case '-lg':
    await df.columns.addNewCalculated(activityColumnScaled, '-1*Log10(${' + activityColumn + '})');
    splitSeqDf.columns.add(df.getCol(activityColumnScaled));
    break;
  default:
    await df.columns.addNewCalculated(activityColumnScaled, '${' + activityColumn + '}');
    splitSeqDf.columns.add(df.getCol(activityColumnScaled));
    break;
  }

  const positionColName = 'position';
  const aminoAcidResidue = 'aminoAcidResidue';
  const medianColName = 'MAD';

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
    pValues.push(testResult['p-value more']);
    mDiff.push(testResult['Mean difference']!);
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

  // const aarCategoryOrder = [
  //   '-', //black I guess
  //   'C', 'U', //yellow
  //   'G', 'P', //red
  //   'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', //all_green
  //   'R', 'H', 'K', //light_blue
  //   'D', 'E', //dark_blue
  //   'S', 'T', 'N', 'Q', //orange
  // ];
  let aarCategoryOrder: string[] = ['-'];
  for (const group of ChemPalette.grokGroups) {
    aarCategoryOrder = aarCategoryOrder.concat(group[0]);
  }
  matrixDf.getCol(aminoAcidResidue).setCategoryOrder(aarCategoryOrder);

  // !!! DRAWING PHASE !!!
  //find min and max MAD across all of the dataframe
  const dfMin = jStat.min(mDiff);
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
        const count = statsDf.groupBy(['Count']).where(query).aggregate().get('Count', 0);
        if (count < 5) {
          args.preventDefault();
          args.g.restore();
          return;
        }

        const pVal = statsDf.groupBy(['p-value']).where(query).aggregate().get('p-value', 0);

        let coef;
        if (pVal < 0.01) {
          coef = 1;
        } else if (pVal < 0.05) {
          coef = 2/3;
        } else if (pVal < 0.1) {
          coef = 1/3;
        } else {
          coef = 0.01;
        }

        const rCoef = (args.cell.cell.value - dfMin) / (dfMax - dfMin);

        const maxRadius = 0.9 * (args.bounds.width > args.bounds.height ? args.bounds.height : args.bounds.width) / 2;
        const radius = Math.ceil(maxRadius * rCoef);

        args.g.beginPath();
        args.g.fillStyle = DG.Color.toHtml(DG.Color.scaleColor(
          coef,
          0,
          1,
          undefined,
          [DG.Color.lightLightGray, DG.Color.green],
        ));
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
          if (col === 'Count' && text < 5) {
            return true;
          }

          text = col === 'Count' ? text + ` / ${peptidesCount}` : text;
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
      const toDisplay = [ui.divText(cell.cell.value as string)];
      // eslint-disable-next-line no-unused-vars

      const [_c, aar, _p] = cp.getColorAAPivot(cell.cell.value as string);
      if (aar in ChemPalette.AASmiles) {
        const sketch = grok.chem.svgMol(ChemPalette.AASmiles[aar]);
        toDisplay.push(sketch);
      }
      ui.tooltip.show(ui.divV(toDisplay), x, y);
    }
    return true;
  });

  // Select columns in source table that correspond to the currently clicked cell
  grid.table.onCurrentCellChanged.subscribe((_: any) => {
    if (grid.table.currentCell.value && grid.table.currentCol.name !== aminoAcidResidue) {
      const currentAAR: string = grid.table.get(aminoAcidResidue, grid.table.currentRowIdx);
      const currentPosition = grid.table.currentCol.name;

      const query = `aminoAcidResidue = ${currentAAR} and position = ${currentPosition}`;
      const text = statsDf.groupBy(['Count']).where(query).aggregate().get('Count', 0);
      if (text < 5) {
        return;
      }

      const splitColName = '~splitCol';
      const otherColName = 'Other';

      const bitset = filterMode ? df.filter : df.selection;
      bitset.init((i) => df.get(currentPosition, i) === currentAAR);

      const splitArray: string[] = [];
      for (let i = 0; i < bitset.length; i++) {
        //TODO: generate better label
        splitArray.push(bitset.get(i) ?
          `${currentAAR === '-' ? 'Empty' : currentAAR} - ${currentPosition}` : otherColName);
      }

      const splitCol = DG.Column.fromStrings(splitColName, splitArray);
      //TODO: use replace as soon as it is ready
      if (!df.col(splitColName)) {
        df.columns.add(splitCol);
      } else {
        df.columns.remove(splitColName);
        df.columns.add(splitCol);
      }

      //FIXME: coloring doesn't work now
      // const colorMap: {[index: string]: string | number} = {otherColName: DG.Color.lightGray};
      // colorMap[currentAAR] = cp.getColor(currentAAR);
      // df.getCol(splitColName).colors.setCategorical(colorMap);
      // df.getCol(splitColName).setCategoryOrder([otherColName]);
    }
  });

  return [grid, statsDf];
}
