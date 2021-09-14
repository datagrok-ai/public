import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { splitAlignedPeptides } from '../splitAligned';


function decimalAdjust(type: 'floor' | 'ceil', value: number, exp: number) {
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

export async function describe(df: DG.DataFrame, activityColumn: string) {
  //Split the aligned sequence into separate AARs
  let splitSeqDf;
  for (let col of df.columns) {
    //FIXME: semType is still undefined at this point                          ?
    if (col.semType === 'alignedSequence' || col.name === 'AlignedSequence') {
      splitSeqDf = splitAlignedPeptides(col);
      break;
    }
  }

  if (typeof splitSeqDf === 'undefined') { return null; }

  let positionColumns = splitSeqDf.columns.names();

  splitSeqDf.columns.add(df.getCol(activityColumn));

  let positionColName = 'position';
  let aminoAcidResidue = 'aminoAcidResidue';

  //unpivot a table and handle duplicates
  let matrixDf = splitSeqDf.groupBy(positionColumns)
    .add('med', activityColumn, activityColumn)
    .aggregate()
    .unpivot([activityColumn], positionColumns, positionColName, aminoAcidResidue);
    
  let peptidesCount = splitSeqDf.getCol(activityColumn).length;

  //this table contains overall statistics on activity
  let totalStats = matrixDf.groupBy()
    .add('med', activityColumn, 'med')
    .aggregate()

  //statistics for specific AAR at a specific position
  matrixDf = matrixDf.groupBy([positionColName, aminoAcidResidue])
    .add('med', activityColumn, 'med')
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
    .aggregate()

  // calculate additional stats
  await matrixDf.columns.addNewCalculated('ratio', '100*${count}/'.concat(`${peptidesCount}`));
  await matrixDf.columns.addNewCalculated('cv', '${stdev}/${avg}');
  await matrixDf.columns.addNewCalculated('iqr', '${q3}-${q1}');
  await matrixDf.getCol('med').applyFormula('${med}-'.concat(`${totalStats.get('med', 0)}`));

  let statsDf = matrixDf.clone();

  //pivot a table to make it matrix-like
  matrixDf = matrixDf.groupBy([aminoAcidResidue])
    .pivot(positionColName)
    .add('first', 'ratio', '')
    .aggregate();  
  matrixDf.name = 'SAR';

  // !!! DRAWING PHASE !!!
  let medianColName = 'med';
  
  //find min and max across all of the dataframe
  let dfMin = Infinity;
  let dfMax = -Infinity;
  for (let col of matrixDf.columns) {
    if (col.name === aminoAcidResidue) { continue; }
    dfMin = dfMin > col.min ? col.min : dfMin;
    dfMax = dfMax < col.max ? col.max : dfMax;
  }

  // color-coding cells based on the entire dataframe
  // colors are hard-coded and can be either white, light-green or green
  for (let col of matrixDf.columns) {
    if (col.name === aminoAcidResidue) { continue; }
    let parts = 3;
    let colPartLength = (dfMax - dfMin) / parts;
    let range;
    let color;
    
    let condition = '{';
    for (let part = 0; part < parts; part++) {
      if (part !== 0) { condition = condition.concat(',') }
      // add to upper boundary a bit so it colors the highest values too
      range = `"${dfMin+colPartLength*part}-${dfMin+colPartLength*(part+1) + (part === (parts-1) ? 1 : 0)}"`;
      color = `"rgb(${Math.floor(255-part*(255/(parts-1)))},${Math.floor(255-part*((255-122)/(parts-1)))},${Math.floor(255-part*(255/(parts-1)))})"`
      condition = condition.concat(`${range}:${color}`);
    }
    condition = condition.concat('}');

    // Unable to choose custom colors for linear                               ?
    // col.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
    col.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';

    col.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = condition;
  }

  let grid = matrixDf.plot.grid();

  // set custom text and color in cells
  grid.onCellPrepare(function (gc) {
    if (
        gc.isTableCell
        && gc.tableColumn !== null
        && gc.tableColumn.name !== aminoAcidResidue
        && gc.tableRowIndex !== null
        && gc.cell.value !== null
      ) {
      let cellColor = DG.Color.getCellColor(gc.cell!);
      gc.style.backColor = cellColor;
      
      gc.style.textColor = DG.Color.getContrastColor(cellColor);
      let query = `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, gc.tableRowIndex)} and ${positionColName} = ${gc.tableColumn.name}`;
      gc.customText = `${decimalAdjust('floor', statsDf.groupBy([medianColName]).where(query).aggregate().get(medianColName, 0), -2)}`;
    }
  });

  // render column headers and AAR symbols centered
  grid.onCellRender.subscribe(function (args) {
    if(args.cell.isColHeader){
      let textSize = args.g.measureText(args.cell.gridColumn.name);
      args.g.fillText(
        args.cell.gridColumn.name,
        args.bounds.x + (args.bounds.width - textSize.width)/2,
        args.bounds.y + (textSize.fontBoundingBoxAscent+textSize.fontBoundingBoxDescent)
      );
      args.g.fillStyle = '#4b4b4a';
      args.preventDefault();
    }

    if (args.cell.isTableCell && args.cell.tableColumn !== null && args.cell.tableColumn.name === aminoAcidResidue) {
      let textSize = args.g.measureText(args.cell.cell.value);
      args.g.fillText(
        args.cell.cell.value,
        args.bounds.x + (args.bounds.width - textSize.width)/2,
        args.bounds.y + (textSize.fontBoundingBoxAscent+textSize.fontBoundingBoxDescent)
      );
      args.g.fillStyle = '#4b4b4a';
      args.preventDefault();
    }
  });

  // show all the statistics in a tooltip over cell
  grid.onCellTooltip(function (cell, x, y) {
    if (
        !cell.isRowHeader && !cell.isColHeader 
        && cell.tableColumn !== null 
        && cell.tableColumn.name !== aminoAcidResidue 
        && cell.cell.value !== null 
        && cell.tableRowIndex !== null
      ) {
      let textDivs = [];

      for (let col of statsDf.columns.names()) {
        if (col !== aminoAcidResidue && col !== positionColName) {
          let query = `${aminoAcidResidue} = ${matrixDf.get(aminoAcidResidue, cell.tableRowIndex)} and ${positionColName} = ${cell.tableColumn.name}`;
          let text = `${decimalAdjust('floor', statsDf.groupBy([col]).where(query).aggregate().get(col, 0), -5)}`;
          textDivs.push(ui.divText(`${col}: ${text}`));
        }
      }

      ui.tooltip.show(ui.divV(textDivs), x, y);
      return true;
    }
  });

  return grid;
}
