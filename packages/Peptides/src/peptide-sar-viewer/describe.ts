import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {splitAlignedPeptides} from '../split-aligned';


function decimalAdjust(type: 'floor' | 'ceil' | 'round', value: number, exp: number): number {
  // If the exp is undefined or zero...
  if (typeof exp === 'undefined' || +exp === 0) {
    return Math[type](value);
  }
  value = +value;
  exp = +exp;
  // If the value is not a number or the exp is not an integer...
  if (isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0)) {
    return NaN;
  }
  // Shift
  let valueArr = value.toString().split('e');
  value = Math[type](+(valueArr[0] + 'e' + (valueArr[1] ? (+valueArr[1] - exp) : -exp)));
  // Shift back
  valueArr = value.toString().split('e');
  return +(valueArr[0] + 'e' + (valueArr[1] ? (+valueArr[1] + exp) : exp));
}

export async function describe(
  df: DG.DataFrame,
  activityColumn: string,
  activityScaling: string,
  filterMode: boolean,
): Promise<DG.Grid | null> {
  //Split the aligned sequence into separate AARs
  let splitSeqDf: DG.DataFrame | undefined;
  for (const col of df.columns) {
    if (col.semType === 'alignedSequence') {
      splitSeqDf = splitAlignedPeptides(col);
      // splitSeqDf.name = 'splitSeq';
      break;
    }
  }

  if (typeof splitSeqDf === 'undefined') {
    return null;
  }

  const positionColumns = splitSeqDf.columns.names();

  // splitSeqDf.columns.add(df.getCol(activityColumn));

  // append splitSeqDf columns to source table and make sure columns are not added more than once
  const dfColsSet = new Set(df.columns.names());
  if (!positionColumns.every((col: string) => dfColsSet.has(col))) {
    df.join(splitSeqDf, [activityColumn], [activityColumn], df.columns.names(), positionColumns, 'inner', true);
  }
  positionColumns.forEach((name:string)=> {
    const col = df.getCol(name);
    col.semType = 'aminoAcids';
    col.setTag('cell.renderer', 'aminoAcids');
  });


  const activityColumnScaled = activityColumn + 'Scaled';

  // scale activity
  //TODO: how to NOT render these?
  switch (activityScaling) {
  case 'lg':
    await df.columns.addNewCalculated(activityColumnScaled, 'Log10(${' + activityColumn + '})');
    splitSeqDf.columns.add(df.getCol(activityColumnScaled));
    // splitSeqDf.columns.remove(activityColumn);
    // splitSeqDf.getCol('lg').name = activityColumn;
    break;
  case '-lg':
    await df.columns.addNewCalculated(activityColumnScaled, '-1*Log10(${' + activityColumn + '})');
    splitSeqDf.columns.add(df.getCol(activityColumnScaled));
    // splitSeqDf.columns.remove(activityColumn);
    // splitSeqDf.getCol('-lg').name = activityColumn;
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
  let matrixDf = splitSeqDf.groupBy(positionColumns)
    .add('med', activityColumnScaled, activityColumnScaled)
    .aggregate()
    .unpivot([activityColumnScaled], positionColumns, positionColName, aminoAcidResidue);

  const peptidesCount = splitSeqDf.getCol(activityColumnScaled).length;

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
    .add('med', 'innerMAD', medianColName) //final step of MAD calculation
    .add('q1', activityColumnScaled, 'Q1')
    .add('q2', activityColumnScaled, 'Median')
    .add('q3', activityColumnScaled, 'Q3')
    .aggregate();

  // calculate additional stats
  await matrixDf.columns.addNewCalculated('IQR', '${q3}-${q1}');
  await matrixDf.columns.addNewCalculated('CQV', '(${q3}-${q1})/(${q3}+${q1})');
  await matrixDf.columns.addNewCalculated('Ratio', '${count}/'.concat(`${peptidesCount}`));

  const statsDf = matrixDf.clone();

  //pivot a table to make it matrix-like
  matrixDf = matrixDf.groupBy([aminoAcidResidue])
    .pivot(positionColName)
    .add('first', medianColName, '')
    .aggregate();
  matrixDf.name = 'SAR';

  // !!! DRAWING PHASE !!!
  //find min and max MAD across all of the dataframe
  const dfMinMedian = statsDf.getCol(medianColName).min;
  const dfMaxMedian = statsDf.getCol(medianColName).max;
  const grid = matrixDf.plot.grid();

  for (const col of matrixDf.columns) {
    if (col.name === aminoAcidResidue) {
      col.semType = 'aminoAcids';
      col.setTag('cell.renderer', 'aminoAcids');
    }
  }
  grid.columns.setOrder([aminoAcidResidue].concat(positionColumns));

  //render column headers and AAR symbols centered
  grid.onCellRender.subscribe(function(args: DG.GridCellRenderArgs) {
    if (args.cell.isColHeader) {
      const textSize = args.g.measureText(args.cell.gridColumn.name);
      if ( args.cell.gridColumn.name != aminoAcidResidue) {
        args.g.fillText(
          args.cell.gridColumn.name,
          args.bounds.x + (args.bounds.width - textSize.width) / 2,
          args.bounds.y + (textSize.actualBoundingBoxAscent + textSize.actualBoundingBoxDescent),
        );
        args.g.fillStyle = '#4b4b4a';
      }
      args.preventDefault();
    }

    if (args.cell.isTableCell && args.cell.tableRowIndex !== null && args.cell.tableColumn !== null) {
      if (args.cell.tableColumn.name === aminoAcidResidue) {

      } else if (args.cell.cell.value !== null) {
        const query =
          `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, args.cell.tableRowIndex)} ` +
          `and ${positionColName} = ${args.cell.tableColumn.name}`;
        const ratio = statsDf.groupBy(['ratio']).where(query).aggregate().get('ratio', 0);
        const maxRadius = 0.95 * (args.bounds.width > args.bounds.height ? args.bounds.height : args.bounds.width) / 2;
        const radius = Math.ceil(maxRadius * ratio);

        args.g.beginPath();
        args.g.fillStyle = DG.Color.toHtml(DG.Color.scaleColor(
          args.cell.cell.value,
          dfMinMedian,
          dfMaxMedian,
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
      const tooltipMap = {};

      for (const col of statsDf.columns.names()) {
        if (col !== aminoAcidResidue && col !== positionColName) {
          const query =
            `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, cell.tableRowIndex)} ` +
            `and ${positionColName} = ${cell.tableColumn.name}`;
          const text = `${decimalAdjust('floor', statsDf.groupBy([col]).where(query).aggregate().get(col, 0), -5)}`;
          // @ts-ignore: idk what's wrong with indexing object with a string :/
          tooltipMap[col] = text;
        }
      }

      ui.tooltip.show(ui.tableFromMap(tooltipMap), x, y);
    }
    return true;
  });

  // Select columns in source table that correspond to the currently clicked cell
  grid.table.onCurrentCellChanged.subscribe((_: any) => {
    if (grid.table.currentCell.value !== null && grid.table.currentCol.name !== aminoAcidResidue) {
      const currentAAR = grid.table.get(aminoAcidResidue, grid.table.currentCell.rowIndex);
      const currentPosition = grid.table.currentCol.name;

      // @ts-ignore: I'd love to use row.get(), but unfortunately there's no column 'get' :(
      splitSeqDf!.rows.select((row) => row[currentPosition] === currentAAR);
      const bitset = filterMode ? df.filter : df.selection;
      // bitset.init((i) => splitSeqDf!.selection.get(i));
      bitset.copyFrom(splitSeqDf!.selection);

      if (!df.col('splitCol')) {
        df.columns.add(DG.Column.fromBitSet('splitCol', bitset));
      }
    }
  });

  return grid;
}
