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
  activityScaling: string
  ): Promise<DG.Grid | null> {
  //Split the aligned sequence into separate AARs
  let splitSeqDf: DG.DataFrame | undefined;
  for (const col of df.columns) {
    if (col.semType === 'alignedSequence') {
      splitSeqDf = splitAlignedPeptides(col);
      break;
    }
  }

  if (typeof splitSeqDf === 'undefined') {
    return null;
  }

  const positionColumns = splitSeqDf.columns.names();

  splitSeqDf.columns.add(df.getCol(activityColumn));

  // append splitSeqDf columns to source table and make sure columns are not added more than once
  const dfColsSet = new Set(df.columns.names());
  if (!positionColumns.every((col: string) => dfColsSet.has(col))) {
    df.join(splitSeqDf, [activityColumn], [activityColumn], df.columns.names(), positionColumns, 'inner', true);
  }

  // scale activity
  switch (activityScaling) {
  case 'lg':
    await splitSeqDf.columns.addNewCalculated('lg', 'Log10(${' + activityColumn + '})');
    splitSeqDf.columns.remove(activityColumn);
    splitSeqDf.getCol('lg').name = activityColumn;
    break;
  case '-lg':
    await splitSeqDf.columns.addNewCalculated('-lg', '-1*Log10(${' + activityColumn + '})');
    splitSeqDf.columns.remove(activityColumn);
    splitSeqDf.getCol('-lg').name = activityColumn;
    break;
  default:
    break;
  }

  const positionColName = 'position';
  const aminoAcidResidue = 'aminoAcidResidue';
  const medianColName = 'MAD';

  //unpivot a table and handle duplicates
  let matrixDf = splitSeqDf.groupBy(positionColumns)
    .add('med', activityColumn, activityColumn)
    .aggregate()
    .unpivot([activityColumn], positionColumns, positionColName, aminoAcidResidue);

  const peptidesCount = splitSeqDf.getCol(activityColumn).length;

  //this table contains overall statistics on activity
  const totalStats = matrixDf.groupBy()
    .add('med', activityColumn, 'med')
    .aggregate();

  //preparing for mad
  await matrixDf.columns.addNewCalculated('innerMAD', 'Abs(${' + activityColumn + '}-' + totalStats.get('med', 0) + ')');

  //statistics for specific AAR at a specific position
  matrixDf = matrixDf.groupBy([positionColName, aminoAcidResidue])
    .add('count', activityColumn, 'Count')
    .add('med', 'innerMAD', medianColName)  //final step of MAD calculation
    .add('q1', activityColumn, 'Q1')
    .add('q2', activityColumn, 'Median')
    .add('q3', activityColumn, 'Q3')
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

  // color-coding cells based on the entire dataframe
  // colors are hard-coded and can be either gray, light-green or green
  const parts = 3;
  const colPartLength = (dfMaxMedian - dfMinMedian) / parts;
  // LightGray, LightGreen, Green
  const colors = ['#D3D3D3', '#90EE90', '#008000'];

  let range;
  let condition = '{';
  for (let part = 0; part < parts; part++) {
    if (part !== 0) { condition = condition.concat(',') }
    // add to upper boundary a bit so it colors the highest values too
    range = `"${dfMinMedian+colPartLength*part}-${dfMinMedian+colPartLength*(part+1)+(part === (parts-1) ? 1 : 0)}"`;
    condition = condition.concat(`${range}:"${colors[part]}"`);
  }
  condition = condition.concat('}');

  //setting color coding tags
  for (const col of matrixDf.columns) {
    if (col.name === aminoAcidResidue) { continue; }

    col.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
    col.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = condition;
  }

  const grid = matrixDf.plot.grid();

  grid.columns.setOrder([aminoAcidResidue].concat(positionColumns));

  // // render column headers and AAR symbols centered
  grid.onCellRender.subscribe(function(args: DG.GridCellRenderArgs) {
    if (args.cell.isColHeader) {
      const textSize = args.g.measureText(args.cell.gridColumn.name);
      args.g.fillText(
        args.cell.gridColumn.name,
        args.bounds.x + (args.bounds.width - textSize.width) / 2,
        args.bounds.y + (textSize.actualBoundingBoxAscent + textSize.actualBoundingBoxDescent),
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
          //FIXME: the text is too high in the cell
          args.bounds.y + (textSize.actualBoundingBoxAscent + textSize.actualBoundingBoxDescent + args.bounds.height) / 2,
        );
        args.g.fillStyle = '#4b4b4a';
        args.preventDefault();
      } else if (args.cell.cell.value !== null) {
        const query = `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, args.cell.tableRowIndex)} and ${positionColName} = ${args.cell.tableColumn.name}`;
        const ratio = statsDf.groupBy(['ratio']).where(query).aggregate().get('ratio', 0);

        args.g.beginPath();
        const maxRadius = 0.95 * (args.bounds.width > args.bounds.height ? args.bounds.height : args.bounds.width) / 2;
        const radius = Math.ceil(maxRadius * ratio);
        args.g.arc(
          args.bounds.x + args.bounds.width / 2,
          args.bounds.y + args.bounds.height / 2,
          radius < 3 ? 3 : radius,
          0,
          Math.PI * 2,
          true,
        );
        args.g.closePath();

        args.g.fillStyle = DG.Color.getCellColorHtml(args.cell.cell);
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
          const query = `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, cell.tableRowIndex)} and ${positionColName} = ${cell.tableColumn.name}`;
          const text = `${decimalAdjust('floor', statsDf.groupBy([col]).where(query).aggregate().get(col, 0), -5)}`;
          // textDivs.push(ui.divText(`${col}: ${text}`));
          tooltipMap[<string>col] = text;
        }
      }

      ui.tooltip.show(ui.tableFromMap(tooltipMap), x, y);
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
