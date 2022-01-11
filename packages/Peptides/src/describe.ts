import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {splitAlignedPeptides} from './utils/split-aligned';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {fdrcorrection} from '@datagrok-libraries/statistics/src/multiple-tests';
import {ChemPalette} from './utils/chem-palette';
import {setAARRenderer} from './utils/cell-renderer';

const cp = new ChemPalette('grok');

export const aarGroups = {
  'R': 'PC',
  'H': 'PC',
  'K': 'PC',
  'D': 'NC',
  'E': 'NC',
  'S': 'U',
  'T': 'U',
  'N': 'U',
  'Q': 'U',
  'C': 'SC',
  'U': 'SC',
  'G': 'SC',
  'P': 'SC',
  'A': 'H',
  'V': 'H',
  'I': 'H',
  'L': 'H',
  'M': 'H',
  'F': 'H',
  'Y': 'H',
  'W': 'H',
  '-': '-',
};

const groupDescription: {[key: string]: {'description': string, 'aminoAcids': string[]}} = {
  'PC': {'description': 'Positive Amino Acids, with Electrically Charged Side Chains', 'aminoAcids': ['R', 'H', 'K']},
  'NC': {'description': 'Negative Amino Acids, with Electrically Charged Side Chains', 'aminoAcids': ['D', 'E']},
  'U': {'description': 'Amino Acids with Polar Uncharged Side Chains', 'aminoAcids': ['S', 'T', 'N', 'Q']},
  'SC': {'description': 'Special Cases', 'aminoAcids': ['C', 'U', 'G', 'P']},
  'H': {
    'description': 'Amino Acids with Hydrophobic Side Chain',
    'aminoAcids': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
  },
  '-': {'description': 'Unknown Amino Acid', 'aminoAcids': ['-']},
};

/*function customGridColumnHeader(cell: DG.GridCell) {
  if (cell.isColHeader && cell.tableColumn != null) {
    if (highlightedColumns.includes(parseInt(cell.tableColumn.name))) {
      cell.style.backColor = 0xff1f77b4;
    }
  }
}*/

function joinDataFrames(
  activityColumnScaled: string,
  df: DG.DataFrame,
  positionColumns: string[],
  splitSeqDf: DG.DataFrame,
  activityColumn: string,
) {
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
}

function sortSourceGrid(sourceGrid: DG.Grid) {
  if (sourceGrid) {
    const colNames:string[] = [];
    for (let i = 0; i < sourceGrid.columns.length; i++) {
      colNames.push(sourceGrid.columns.byIndex(i)!.name);
    }
    colNames.sort((a, b)=>{
      if (sourceGrid.columns.byName(a)?.column?.semType == 'aminoAcids') {
        if (sourceGrid.columns.byName(b)?.column?.semType == 'aminoAcids') {
          return 0;
        }
        return -1;
      }
      if (sourceGrid.columns.byName(b)?.column?.semType == 'aminoAcids') {
        return 1;
      }
      return 0;
    });
    sourceGrid?.columns.setOrder(colNames);
  }
}

async function scaleActivity(
  activityScaling: string,
  activityColumn: string,
  activityColumnScaled: string,
  sourceGrid: DG.Grid,
  splitSeqDf: DG.DataFrame,
) {
  const df = sourceGrid.dataFrame!;
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
}

async function calculateStatistics(
  matrixDf: DG.DataFrame,
  positionColName: string,
  aminoAcidResidue: string,
  activityColumnScaled: string,
  peptidesCount: number,
  splitSeqDf: DG.DataFrame,
  groupMapping: {[key: string]: string},
) {
  matrixDf = matrixDf.groupBy([positionColName, aminoAcidResidue])
    .add('count', activityColumnScaled, 'Count')
    .aggregate();

  const countThreshold = 4;
  //@ts-ignore: never gets old
  matrixDf.rows.filter((row) => row.Count >= countThreshold && row.Count <= peptidesCount - countThreshold);
  matrixDf = matrixDf.clone(matrixDf.filter);

  // calculate additional stats
  await matrixDf.columns.addNewCalculated('Ratio', '${count}/'.concat(`${peptidesCount}`));

  //calculate p-values based on t-test
  let pvalues: Float32Array = new Float32Array(matrixDf.rowCount).fill(1);
  const mdCol: DG.Column = matrixDf.columns.addNewFloat('Mean difference');
  const pValCol: DG.Column = matrixDf.columns.addNewFloat('pValue');
  for (let i = 0; i < matrixDf.rowCount; i++) {
    const position = matrixDf.get(positionColName, i);
    const aar = matrixDf.get(aminoAcidResidue, i);

    //@ts-ignore
    splitSeqDf.rows.select((row) => groupMapping[row[position]] === aar);
    const currentActivity: number[] = splitSeqDf
      .clone(splitSeqDf.selection, [activityColumnScaled])
      .getCol(activityColumnScaled)
      .toList();

    //@ts-ignore
    splitSeqDf.rows.select((row) => groupMapping[row[position]] !== aar);
    const otherActivity: number[] = splitSeqDf
      .clone(splitSeqDf.selection, [activityColumnScaled])
      .getCol(activityColumnScaled)
      .toList();

    const testResult = tTest(currentActivity, otherActivity);
    // testResult = uTest(currentActivity, otherActivity);
    const currentMeanDiff = testResult['Mean difference']!;
    const pvalue = testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'];

    mdCol.set(i, currentMeanDiff);
    pvalues[i] = pvalue;
  }

  if (true) {
    pvalues = fdrcorrection(pvalues)[1];
  }

  for (let i = 0; i < pvalues.length; ++i) {
    pValCol.set(i, pvalues[i]);
  }

  return matrixDf.clone();
}

async function setCategoryOrder(
  twoColorMode: boolean, statsDf: DG.DataFrame, aminoAcidResidue: string, matrixDf: DG.DataFrame,
) {
  const sortArgument = twoColorMode ? 'Absolute Mean difference' : 'Mean difference';
  if (twoColorMode) {
    await statsDf.columns.addNewCalculated('Absolute Mean difference', 'Abs(${Mean difference})');
  }
  const aarWeightsDf = statsDf.groupBy([aminoAcidResidue]).sum(sortArgument, 'weight').aggregate();
  const aarList = aarWeightsDf.getCol(aminoAcidResidue).toList();
  const getWeight = (aar: string) => aarWeightsDf
    .groupBy(['weight'])
    .where(`${aminoAcidResidue} = ${aar}`)
    .aggregate()
    .get('weight', 0);
  aarList.sort((first, second) => getWeight(second) - getWeight(first));

  matrixDf.getCol(aminoAcidResidue).setCategoryOrder(aarList);
}

function createVerticalTable(
  statsDf: DG.DataFrame,
  aminoAcidResidue: string,
  positionColName: string,
  twoColorMode: boolean,
) {
  // TODO: aquire ALL of the positions
  let sequenceDf = statsDf.groupBy(['Mean difference', aminoAcidResidue, positionColName, 'Count', 'Ratio', 'pValue'])
    .where('pValue <= 0.1')
    .aggregate();

  let tempStats: DG.Stats;
  const maxAtPos: {[index: string]: number} = {};
  for (const pos of sequenceDf.getCol(positionColName).categories) {
    tempStats = DG.Stats.fromColumn(
      sequenceDf.getCol('Mean difference'),
      DG.BitSet.create(sequenceDf.rowCount, (i) => sequenceDf.get(positionColName, i) === pos),
    );
    maxAtPos[pos] = twoColorMode ?
      (tempStats.max > Math.abs(tempStats.min) ? tempStats.max : tempStats.min) : tempStats.max;
  }
  sequenceDf = sequenceDf.clone(DG.BitSet.create(sequenceDf.rowCount, (i) => {
    return sequenceDf.get('Mean difference', i) === maxAtPos[sequenceDf.get(positionColName, i)];
  }));

  return sequenceDf;
}

function createGrids(
  matrixDf: DG.DataFrame,
  aminoAcidResidue: string,
  positionColumns: string[],
  sequenceDf: DG.DataFrame,
  positionColName: string,
  grouping: boolean,
) {
  const sarGrid = matrixDf.plot.grid();
  sarGrid.sort([aminoAcidResidue]);
  sarGrid.columns.setOrder([aminoAcidResidue].concat(positionColumns));

  const sarVGrid = sequenceDf.plot.grid();
  sarVGrid.sort([positionColName]);
  sarVGrid.col('pValue')!.format = 'four digits after comma';
  sarVGrid.col('pValue')!.name = 'P-Value';

  if (!grouping) {
    let tempCol = matrixDf.columns.byName(aminoAcidResidue);
    if (tempCol) {
      setAARRenderer(tempCol, sarGrid);
    }
    tempCol = sequenceDf.columns.byName(aminoAcidResidue);
    if (tempCol) {
      setAARRenderer(tempCol, sarGrid);
    }
  }

  return [sarGrid, sarVGrid];
}

function setCellRendererFunc(
  renderColNames: string[],
  positionColName: string,
  aminoAcidResidue: string,
  statsDf: DG.DataFrame,
  twoColorMode: boolean,
  sarGrid: DG.Grid,
  sarVGrid: DG.Grid,
) {
  const mdCol = statsDf.getCol('Mean difference');
  const cellRendererFunc = function(args: DG.GridCellRenderArgs) {
    args.g.save();
    args.g.beginPath();
    args.g.rect(args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height);
    args.g.clip();

    if (args.cell.isRowHeader && args.cell.gridColumn.visible) {
      args.cell.gridColumn.visible = false;
      args.preventDefault();
      return;
    }

    if (
      args.cell.isTableCell &&
      args.cell.tableRowIndex !== null &&
      args.cell.tableColumn !== null &&
      args.cell.cell.value !== null
    ) {
      if (renderColNames.indexOf(args.cell.tableColumn.name) !== -1) {
        const currentPosition = args.cell.tableColumn.name !== 'Mean difference' ?
          args.cell.tableColumn.name : args.cell.grid.table.get(positionColName, args.cell.tableRowIndex);
        const query =
          `${aminoAcidResidue} = ${args.cell.grid.table.get(aminoAcidResidue, args.cell.tableRowIndex)} ` +
          `and ${positionColName} = ${currentPosition}`;

        const pVal: number = statsDf.groupBy(['pValue']).where(query).aggregate().get('pValue', 0);

        let coef;
        const variant = args.cell.cell.value < 0;
        if (pVal < 0.01) {
          coef = variant && twoColorMode ? '#FF7900' : '#299617';
        } else if (pVal < 0.05) {
          coef = variant && twoColorMode ? '#FFA500' : '#32CD32';
        } else if (pVal < 0.1) {
          coef = variant && twoColorMode ? '#FBCEB1' : '#98FF98';
        } else {
          coef = DG.Color.toHtml(DG.Color.lightLightGray);
        }

        const chooseMin = () => twoColorMode ? 0 : mdCol.min;
        const chooseMax = () => twoColorMode ? Math.max(Math.abs(mdCol.min), mdCol.max) : mdCol.max;
        const chooseCurrent = () => twoColorMode ? Math.abs(args.cell.cell.value) : args.cell.cell.value;

        const rCoef = (chooseCurrent() - chooseMin()) / (chooseMax() - chooseMin());

        const maxRadius = 0.9 * (args.bounds.width > args.bounds.height ? args.bounds.height : args.bounds.width) / 2;
        const radius = Math.floor(maxRadius * rCoef);

        args.g.beginPath();
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
  };
  sarGrid.onCellRender.subscribe(cellRendererFunc);
  sarVGrid.onCellRender.subscribe(cellRendererFunc);
}

function setTooltipFunc(
  renderColNames: string[],
  statsDf: DG.DataFrame,
  aminoAcidResidue: string,
  positionColName: string,
  peptidesCount: number,
  grouping: boolean,
  sarGrid: DG.Grid,
  sarVGrid: DG.Grid,
) {
  const onCellTooltipFunc = function(cell: DG.GridCell, x: number, y: number) {
    if (
      !cell.isRowHeader &&
      !cell.isColHeader &&
      cell.tableColumn !== null &&
      cell.cell.value !== null &&
      cell.tableRowIndex !== null &&
      renderColNames.indexOf(cell.tableColumn.name) !== -1
    ) {
      const tooltipMap: { [index: string]: string } = {};

      for (const col of statsDf.columns.names()) {
        if (col !== aminoAcidResidue && col !== positionColName) {
          const currentPosition = cell.tableColumn.name !== 'Mean difference' ?
            cell.tableColumn.name : cell.grid.table.get(positionColName, cell.tableRowIndex);
          const query =
            `${aminoAcidResidue} = ${cell.grid.table.get(aminoAcidResidue, cell.tableRowIndex)} ` +
            `and ${positionColName} = ${currentPosition}`;
          const textNum = statsDf.groupBy([col]).where(query).aggregate().get(col, 0);
          let text = `${col === 'Count' ? textNum : textNum.toFixed(5)}`;

          if (col === 'Count') {
            text += ` / ${peptidesCount}`;
          } else if (col === 'pValue') {
            text = parseFloat(text) !== 0 ? text : '<0.01';
          }

          tooltipMap[col === 'pValue' ? 'p-value' : col] = text;
        }
      }

      ui.tooltip.show(ui.tableFromMap(tooltipMap), x, y);
    }
    if (
      !cell.isColHeader &&
      cell.tableColumn !== null &&
      cell.tableColumn.name == aminoAcidResidue &&
      cell.cell.value !== null &&
      cell.tableRowIndex !== null
    ) {
      if (grouping) {
        const currentGroup = groupDescription[cell.cell.value];
        const divText = ui.divText('Amino Acids in this group: ' + currentGroup['aminoAcids'].join(', '));
        ui.tooltip.show(ui.divV([ui.h3(currentGroup['description']), divText]), x, y);
      } else {
        cp.showTooltip(cell, x, y);
      }
    }
    return true;
  };
  sarGrid.onCellTooltip(onCellTooltipFunc);
  sarVGrid.onCellTooltip(onCellTooltipFunc);
}

function postProcessGrids(
  sourceGrid: DG.Grid,
  invalidIndexes: number[],
  matrixDf: DG.DataFrame,
  grouping: boolean,
  aminoAcidResidue: string,
  sarGrid: DG.Grid,
  sarVGrid: DG.Grid,
) {
  sourceGrid.onCellPrepare((cell: DG.GridCell) => {
    const currentRowIndex = cell.tableRowIndex;
    if (currentRowIndex && invalidIndexes.includes(currentRowIndex) && !cell.isRowHeader) {
      cell.style.backColor = DG.Color.lightLightGray;
    }
  });

  for (const col of matrixDf.columns.names()) {
    sarGrid.col(col)!.width = sarGrid.props.rowHeight;
  }

  if (grouping) {
    sarGrid.col(aminoAcidResidue)!.name = 'Groups';
    sarVGrid.col(aminoAcidResidue)!.name = 'Groups';
  }

  sarGrid.props.allowEdit = false;
  sarVGrid.props.allowEdit = false;
}

export async function describe(
  df: DG.DataFrame,
  activityColumn: string,
  activityScaling: string,
  sourceGrid: DG.Grid,
  twoColorMode: boolean,
  initialBitset: DG.BitSet | null,
  grouping: boolean,
): Promise<[DG.Grid, DG.Grid, DG.DataFrame, {[key: string]: string}]> {
  //Split the aligned sequence into separate AARs
  let splitSeqDf: DG.DataFrame | undefined;
  let invalidIndexes: number[];
  const col: DG.Column = df.columns.bySemType('alignedSequence');
  [splitSeqDf, invalidIndexes] = splitAlignedPeptides(col);
  splitSeqDf.name = 'Split sequence';

  const positionColumns = splitSeqDf.columns.names();
  const activityColumnScaled = `${activityColumn}Scaled`;
  const renderColNames: string[] = splitSeqDf.columns.names();
  const positionColName = 'Position';
  const aminoAcidResidue = 'AAR';

  splitSeqDf.columns.add(df.getCol(activityColumn));

  joinDataFrames(activityColumnScaled, df, positionColumns, splitSeqDf, activityColumn);

  for (const col of df.columns) {
    if (splitSeqDf.col(col.name) && col.name != activityColumn) {
      setAARRenderer(col, sourceGrid);
    }
  }

  sortSourceGrid(sourceGrid);

  await scaleActivity(activityScaling, activityColumn, activityColumnScaled, sourceGrid, splitSeqDf);
  splitSeqDf = splitSeqDf.clone(initialBitset);

  //unpivot a table and handle duplicates
  splitSeqDf = splitSeqDf.groupBy(positionColumns)
    .add('med', activityColumnScaled, activityColumnScaled)
    .aggregate();

  const peptidesCount = splitSeqDf.getCol(activityColumnScaled).length;

  let matrixDf = splitSeqDf.unpivot([activityColumnScaled], positionColumns, positionColName, aminoAcidResidue);

  //TODO: move to chem palette
  let groupMapping: {[key: string]: string} = {};
  if (grouping) {
    groupMapping = aarGroups;
    const aarCol = matrixDf.getCol(aminoAcidResidue);
    aarCol.init((index) => groupMapping[aarCol.get(index)[0]] ?? '-');
    aarCol.compact();
  } else {
    Object.keys(aarGroups).forEach((value) => groupMapping[value] = value);
  }

  //statistics for specific AAR at a specific position
  const statsDf = await calculateStatistics(
    matrixDf, positionColName, aminoAcidResidue, activityColumnScaled, peptidesCount, splitSeqDf, groupMapping,
  );

  // SAR matrix table
  //pivot a table to make it matrix-like
  matrixDf = statsDf.groupBy([aminoAcidResidue])
    .pivot(positionColName)
    .add('first', 'Mean difference', '')
    .aggregate();
  matrixDf.name = 'SAR';

  // Setting category order
  await setCategoryOrder(twoColorMode, statsDf, aminoAcidResidue, matrixDf);

  // SAR vertical table (naive, choose best Mean difference from pVals <= 0.01)
  const sequenceDf = createVerticalTable(statsDf, aminoAcidResidue, positionColName, twoColorMode);
  renderColNames.push('Mean difference');

  const [sarGrid, sarVGrid] = createGrids(
    matrixDf, aminoAcidResidue, positionColumns, sequenceDf, positionColName, grouping,
  );

  setCellRendererFunc(
    renderColNames, positionColName, aminoAcidResidue, statsDf, twoColorMode, sarGrid, sarVGrid,
  );

  // show all the statistics in a tooltip over cell
  setTooltipFunc(
    renderColNames, statsDf, aminoAcidResidue, positionColName, peptidesCount, grouping, sarGrid, sarVGrid,
  );

  postProcessGrids(sourceGrid, invalidIndexes, matrixDf, grouping, aminoAcidResidue, sarGrid, sarVGrid);

  //TODO: return class instead
  return [sarGrid, sarVGrid, statsDf, groupMapping];
}
