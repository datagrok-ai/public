import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from '../utils/constants';
import * as type from '../utils/types';
import {addExpandIcon, getSeparator, setGridProps} from '../utils/misc';
import {renderCellSelection} from '../utils/cell-renderer';

export type MutationCliffsOptions = {
  mutationCliffs: type.MutationCliffs, mutationCliffsSelection: type.Selection, sequenceColumnName: string,
  positionColumns: DG.Column<string>[], gridColumns: DG.GridColumnList, activityCol: DG.Column<number>,
};

/**
 * Creates mutation cliffs widget that shows grid of sequences that form mutation cliffs pairs as well as
 * grid of mutation cliffs pairs themselves.
 * @param table - table with sequences.
 * @param options - options for mutation cliffs widget.
 * @return - mutation cliffs widget.
 */
export function mutationCliffsWidget(table: DG.DataFrame, options: MutationCliffsOptions): DG.Widget {
  const filteredIndexes = table.filter.getSelectedIndexes();
  const positions = Object.keys(options.mutationCliffsSelection);

  if (!positions.length || options.mutationCliffs === null)
    return new DG.Widget(ui.label('No mutations table generated'));


  const substitutionsArray: string[] = [];
  const deltaArray: number[] = [];
  const substitutedToArray: string[] = [];
  const fromIdxArray: number[] = [];
  const toIdxArray: number[] = [];
  const alignedSeqCol = table.getCol(options.sequenceColumnName!);
  const alignedSeqColCategories = alignedSeqCol.categories;
  const alignedSeqColData = alignedSeqCol.getRawData();
  // const activityScaledCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY);
  const activityScaledColData = options.activityCol.getRawData();
  const seenIndexes = new Map<number, number[]>();
  const uniqueSequencesBitSet = DG.BitSet.create(table.rowCount);

  const positionColumns: { [colName: string]: DG.Column<string> } =
    Object.fromEntries(options.positionColumns.map((col) => [col.name, col]));
  for (const pos of positions) {
    const posCol = positionColumns[pos];
    const posColCategories = posCol.categories;
    const posColData = posCol.getRawData();

    for (const monomer of options.mutationCliffsSelection[pos]) {
      const substitutionsMap = options.mutationCliffs.get(monomer)
        ?.get(pos) as Map<number, type.UTypedArray> | undefined;
      if (typeof substitutionsMap === 'undefined')
        continue;


      for (const [referenceIdx, indexArray] of substitutionsMap.entries()) {
        if (!filteredIndexes.includes(referenceIdx))
          continue;


        const forbiddentIndexes = seenIndexes.get(referenceIdx) ?? [];
        const baseSequence = alignedSeqColCategories[alignedSeqColData[referenceIdx]];
        const baseActivity = activityScaledColData[referenceIdx];

        for (const subIdx of indexArray) {
          if (forbiddentIndexes.includes(subIdx) || !filteredIndexes.includes(subIdx))
            continue;


          if (!seenIndexes.has(subIdx))
            seenIndexes.set(subIdx, []);


          const subSeq = alignedSeqColCategories[alignedSeqColData[subIdx]];

          seenIndexes.get(subIdx)!.push(referenceIdx);
          substitutionsArray.push(`${baseSequence}#${subSeq}`);
          deltaArray.push(baseActivity - activityScaledColData[subIdx]);
          substitutedToArray.push(posColCategories[posColData[subIdx]]);
          fromIdxArray.push(referenceIdx);
          toIdxArray.push(subIdx);
          uniqueSequencesBitSet.set(referenceIdx, true);
          uniqueSequencesBitSet.set(subIdx, true);
        }
      }
    }
  }

  if (substitutionsArray.length === 0)
    return new DG.Widget(ui.label('No mutations table generated'));


  const substCol = DG.Column.fromStrings('Mutation', substitutionsArray);
  const activityDeltaCol = DG.Column.fromList('double', 'Delta', deltaArray);
  const hiddenSubstToAarCol = DG.Column.fromStrings('~to', substitutedToArray);
  const toIdxCol = DG.Column.fromList(DG.COLUMN_TYPE.INT, '~toIdx', toIdxArray);
  const fromIdxCol = DG.Column.fromList(DG.COLUMN_TYPE.INT, '~fromIdx', fromIdxArray);
  const pairsTable = DG.DataFrame.fromColumns([substCol, activityDeltaCol, hiddenSubstToAarCol, toIdxCol, fromIdxCol]);
  pairsTable.name = 'Mutation Cliff pairs';

  const aminoToInput = ui.stringInput('Mutated to:', '', () => {
    const substitutedToAar = aminoToInput.stringValue;
    if (substitutedToAar !== '')
      pairsTable.filter.init((idx) => hiddenSubstToAarCol.get(idx) === substitutedToAar);
    else
      pairsTable.filter.setAll(true);
  });
  aminoToInput.setTooltip('Filter the rows by the monomer that the mutation was substituted to');

  const pairsGrid = pairsTable.plot.grid();
  setGridProps(pairsGrid);
  substCol.semType = C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;
  substCol.tags[C.TAGS.SEPARATOR] = getSeparator(alignedSeqCol);
  substCol.tags[DG.TAGS.UNITS] = alignedSeqCol.tags[DG.TAGS.UNITS];
  substCol.tags[DG.TAGS.CELL_RENDERER] = 'MacromoleculeDifference';

  let keyPress = false;
  let lastSelectedIndex: number | null = null;
  const pairsSelectedIndexes: number[] = [];
  pairsGrid.onCurrentCellChanged.subscribe((gridCell: DG.GridCell) => {
    try {
      const rowIdx = gridCell.tableRowIndex;
      if (!keyPress)
        return;


      if (rowIdx === null)
        return;


      if (lastSelectedIndex !== null)
        pairsSelectedIndexes.splice(pairsSelectedIndexes.indexOf(lastSelectedIndex), 1);


      if (!pairsSelectedIndexes.includes(rowIdx)) {
        pairsSelectedIndexes.push(rowIdx);
        pairsGrid.invalidate();
      }
      uniqueSequencesTable.filter.fireChanged();
    } finally {
      keyPress = false;
      lastSelectedIndex = gridCell.tableRowIndex;
    }
  });
  pairsGrid.root.addEventListener('keydown', (event) => keyPress = event.key.startsWith('Arrow'));
  pairsGrid.root.addEventListener('click', (event) => {
    const gridCell = pairsGrid.hitTest(event.offsetX, event.offsetY);
    if (!gridCell || gridCell.tableRowIndex === null)
      return;


    const rowIdx = gridCell.tableRowIndex;
    if (!event.shiftKey) {
      pairsSelectedIndexes.length = 0;
      pairsSelectedIndexes.push(rowIdx);
    } else {
      const rowIdxIdx = pairsSelectedIndexes.indexOf(rowIdx);
      if (rowIdxIdx === -1)
        pairsSelectedIndexes.push(rowIdx);
      else
        pairsSelectedIndexes.splice(rowIdxIdx, 1);
    }
    uniqueSequencesTable.filter.fireChanged();
    pairsGrid.invalidate();
  });
  pairsGrid.onCellRender.subscribe((gcArgs) => {
    if (gcArgs.cell.tableColumn?.name !== substCol.name || !pairsSelectedIndexes.includes(gcArgs.cell.tableRowIndex!))
      return;


    renderCellSelection(gcArgs.g, gcArgs.bounds);
  });

  const originalGridColCount = options.gridColumns.length;
  const columnNames: string[] = [];
  for (let colIdx = 1; colIdx < originalGridColCount; colIdx++) {
    const gridCol = options.gridColumns.byIndex(colIdx);
    if (gridCol?.name === options.sequenceColumnName ||
      (gridCol?.visible === true && !Object.hasOwn(positionColumns, gridCol.name)))
      columnNames.push(gridCol!.name);
  }

  const uniqueSequencesTable = table.clone(uniqueSequencesBitSet, columnNames);
  uniqueSequencesTable.name = 'Unique sequences that form Mutation Cliffs pairs';
  const seqIdxCol = uniqueSequencesTable.columns.addNewInt('~seqIdx');
  const seqIdxColData = seqIdxCol.getRawData();
  const selectedIndexes = uniqueSequencesBitSet.getSelectedIndexes();
  seqIdxCol.init((idx) => selectedIndexes[idx]);
  const uniqueSequencesGrid = uniqueSequencesTable.plot.grid();
  setGridProps(uniqueSequencesGrid);
  uniqueSequencesGrid.props.rowHeight = 20;
  uniqueSequencesTable.filter.onChanged.subscribe(() => {
    const uniqueSelectedIndexes: number[] = [];
    for (const idx of pairsSelectedIndexes) {
      uniqueSelectedIndexes.push(fromIdxCol.get(idx)!);
      uniqueSelectedIndexes.push(toIdxCol.get(idx)!);
    }
    uniqueSequencesTable.filter.init(
      (idx) => pairsSelectedIndexes.length === 0 || uniqueSelectedIndexes.includes(seqIdxColData[idx]), false);
  });

  addExpandIcon(pairsGrid);
  addExpandIcon(uniqueSequencesGrid);

  return new DG.Widget(ui.divV([aminoToInput.root, pairsGrid.root, uniqueSequencesGrid.root],
    {style: {width: '100%'}}));
}
