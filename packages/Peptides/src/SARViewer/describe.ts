import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {splitAlignedPeptides} from '../splitAligned';


function decimalAdjust(type: 'floor' | 'ceil' | 'round', value: number, exp: number) {
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

export async function describe(df: DG.DataFrame, activityColumn: string, activityScaling: string) {
  //Split the aligned sequence into separate AARs
  let splitSeqDf: DG.DataFrame | undefined;
  for (const col of df.columns) {
    //FIXME: semType is still undefined at this point                          ?
    if (col.semType === 'alignedSequence' || col.name === 'AlignedSequence') {
      splitSeqDf = splitAlignedPeptides(col);
      break;
    }
  }

  if (typeof splitSeqDf === 'undefined') {
    return null;
  }

  const positionColumns = splitSeqDf.columns.names();

  splitSeqDf.columns.add(df.getCol(activityColumn));

  switch (activityScaling) {
  case 'ln':
    await splitSeqDf.columns.addNewCalculated('ln', 'Ln(${' + activityColumn + '})');
    splitSeqDf.columns.remove(activityColumn);
    splitSeqDf.getCol('ln').name = activityColumn;
    break;
  case '-ln':
    await splitSeqDf.columns.addNewCalculated('-ln', '-1*Ln(${' + activityColumn + '})');
    splitSeqDf.columns.remove(activityColumn);
    splitSeqDf.getCol('-ln').name = activityColumn;
    break;
  default:
    break;
  }

  const positionColName = 'position';
  const aminoAcidResidue = 'aminoAcidResidue';
  const medianColName = 'med';

  //unpivot a table and handle duplicates
  let matrixDf = splitSeqDf.groupBy(positionColumns)
    .add('med', activityColumn, activityColumn)
    .aggregate()
    .unpivot([activityColumn], positionColumns, positionColName, aminoAcidResidue);

  const peptidesCount = splitSeqDf.getCol(activityColumn).length;

  //this table contains overall statistics on activity
  const totalStats = matrixDf.groupBy()
    .add('med', activityColumn, medianColName)
    .aggregate();

  //statistics for specific AAR at a specific position
  matrixDf = matrixDf.groupBy([positionColName, aminoAcidResidue])
    .add('med', activityColumn, medianColName)
    .add('avg', activityColumn, 'avg')
    .add('min', activityColumn, 'min')
    .add('max', activityColumn, 'max')
    .add('q1', activityColumn, 'q1')
    .add('q2', activityColumn, 'q2')
    .add('q3', activityColumn, 'q3')
    .add('count', activityColumn, 'count')
    .add('variance', activityColumn, 'variance')
    .add('skew', activityColumn, 'skew')
    .add('kurt', activityColumn, 'kurt')
    .add('stdev', activityColumn, 'stdev')
    .add('sum', activityColumn, 'sum')
    .aggregate();

  // calculate additional stats
  await matrixDf.columns.addNewCalculated('ratio', '${count}/'.concat(`${peptidesCount}`));
  await matrixDf.columns.addNewCalculated('cv', '${stdev}/${avg}');
  await matrixDf.columns.addNewCalculated('iqr', '${q3}-${q1}');
  await matrixDf.getCol(medianColName).applyFormula('${med}-'.concat(`${totalStats.get(medianColName, 0)}`));

  const statsDf = matrixDf.clone();

  //pivot a table to make it matrix-like
  matrixDf = matrixDf.groupBy([aminoAcidResidue])
    .pivot(positionColName)
    .add('first', 'ratio', '')
    .aggregate();
  matrixDf.name = 'SAR';

  // !!! DRAWING PHASE !!!
  const grid = matrixDf.plot.grid();

  // render column headers and AAR symbols centered
  grid.onCellRender.subscribe(function(args: DG.GridCellRenderArgs) {
    if (args.cell.isColHeader) {
      const textSize = args.g.measureText(args.cell.gridColumn.name);
      args.g.fillText(
        args.cell.gridColumn.name,
        args.bounds.x + (args.bounds.width - textSize.width) / 2,
        args.bounds.y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent),
      );
      args.g.fillStyle = '#4b4b4a';
      args.preventDefault();
    }

    if (args.cell.isTableCell && args.cell.tableRowIndex !== null && args.cell.tableColumn !== null) {
      if (args.cell.tableColumn.name === aminoAcidResidue) {
        const textSize = args.g.measureText(args.cell.cell.value);
        args.g.fillText(
          args.cell.cell.value,
          args.bounds.x + (args.bounds.width - textSize.width) / 2,
          args.bounds.y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent),
        );
        args.g.fillStyle = '#4b4b4a';
        args.preventDefault();
      } else {
        args.g.beginPath();
        args.g.arc(
          args.bounds.x + args.bounds.width / 2,
          args.bounds.y + args.bounds.height / 2,
          Math.ceil(10 * args.cell.cell.value),
          0,
          Math.PI * 2,
          true,
        );
        args.g.closePath();
        //TODO: set color based on activity medians
        args.g.fillStyle = 'green';
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
      const textDivs = [];

      for (const col of statsDf.columns.names()) {
        if (col !== aminoAcidResidue && col !== positionColName) {
          const query = `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, cell.tableRowIndex)} and ${positionColName} = ${cell.tableColumn.name}`;
          const text = `${decimalAdjust('floor', statsDf.groupBy([col]).where(query).aggregate().get(col, 0), -5)}`;
          textDivs.push(ui.divText(`${col}: ${text}`));
        }
      }

      ui.tooltip.show(ui.divV(textDivs), x, y);
      return true;
    }
  });

  // Select columns in source table that correspond to the currently clicked cell
  grid.table.onCurrentCellChanged.subscribe((_: any) => {
    if (grid.table.currentCell.value !== null && grid.table.currentCol.name !== aminoAcidResidue) {
      const currentAAR = grid.table.get(aminoAcidResidue, grid.table.currentCell.rowIndex);
      const currentPosition = grid.table.currentCol.name;

      // @ts-ignore: I'd love to use row.get(), but unfortunately there's no column 'get' :(
      splitSeqDf!.rows.select((row) => row[currentPosition] === currentAAR);
      df.selection.init((i) => splitSeqDf!.selection.get(i));
    }
  });

  return grid;
}
