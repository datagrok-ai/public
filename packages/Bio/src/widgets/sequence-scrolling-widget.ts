/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import { MSAScrollingHeader } from '@datagrok-libraries/bio/src/utils/sequence-position-scroller';
import { MonomerPlacer } from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import { ALPHABET, TAGS as bioTAGS, SplitterFunc } from '@datagrok-libraries/bio/src/utils/macromolecule';
import { _package } from '../package';
import { ISeqSplitted } from '@datagrok-libraries/bio/src/utils/macromolecule/types';


/**
 * Calculate conservation scores for each position in a multiple sequence alignment
 * @param {string[]} sequences - Array of sequence strings
 * @param {SplitterFunc} splitter - Function to split sequences into monomers
 * @return {number[]} Array of conservation scores (0-1) for each position
 */
function calculateConservation(sequences: string[], splitter: SplitterFunc): number[] {
  if (!sequences || sequences.length <= 1) return [];

  // Split sequences using the specialized splitter function
  const splitSequences: ISeqSplitted[] = sequences
    .map((seq) => splitter(seq))
    .filter(Boolean);

  if (splitSequences.length <= 1) return [];

  // Find the maximum sequence length
  const maxLen = Math.max(...splitSequences.map((seq) => seq.length));
  const result: number[] = new Array(maxLen).fill(0);

  // For each position in the alignment
  for (let pos = 0; pos < maxLen; pos++) {
    const residueCounts: Record<string, number> = {};
    let totalResidues = 0;

    // Count the residues at this position
    for (const seq of splitSequences) {
      if (pos < seq.length) {
        // Use getCanonical method to get standardized residue representation
        const residue = seq.getCanonical(pos);

        // Skip gaps - using the interface's method to check
        if (residue && !seq.isGap(pos)) {
          residueCounts[residue] = (residueCounts[residue] || 0) + 1;
          totalResidues++;
        }
      }
    }

    if (totalResidues === 0) {
      result[pos] = 0;
      continue;
    }

    // Calculate conservation score - frequency of most common residue
    const mostCommonCount = Math.max(...Object.values(residueCounts), 0);
    const conservation = mostCommonCount / totalResidues;

    result[pos] = conservation;
  }

  return result;
}


export function handleSequenceHeaderRendering() {
  const handleGrid = (grid: DG.Grid) => {
    setTimeout(() => {
      if (grid.isDetached)
        return;
      const df = grid.dataFrame;
      if (!df)
        return;
      const seqCols = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);

      for (const seqCol of seqCols) {
        // TODO: Extend this to non-canonical monomers and non msa
        const sh = _package.seqHelper.getSeqHandler(seqCol);
        if (!sh)
          continue;
        if (sh.isHelm() || sh.alphabet === ALPHABET.UN)
          continue;

        const gCol = grid.col(seqCol.name);
        if (!gCol)
          continue;

        let positionStatsViewerAddedOnce = !!grid.tableView && Array.from(grid.tableView.viewers).some((v) => v.type === 'Sequence Position Statistics');
        const isMSA = sh.isMsa();
        const ifNan = (a: number, els: number) => (Number.isNaN(a) ? els : a);
        const getStart = () => ifNan(Math.max(Number.parseInt(seqCol.getTag(bioTAGS.positionShift) ?? '0'), 0), 0) + 1;
        const getCurrent = () => ifNan(Number.parseInt(seqCol.getTag(bioTAGS.selectedPosition) ?? '-2'), -2);
        const getFontSize = () => MonomerPlacer.getFontSettings(seqCol).fontWidth;
        // get the maximum length of seqs by getting the primitive length first and then splitting it with correct splitter;
        let pseudoMaxLenIndex = 0;
        let pseudoMaxLength = 0;
        const cats = seqCol.categories;
        for (let i = 0; i < cats.length; i++) {
          const seq = cats[i];
          if (seq && seq.length > pseudoMaxLength) {
            pseudoMaxLength = seq.length;
            pseudoMaxLenIndex = i;
          }
        }
        const seq = cats[pseudoMaxLenIndex];
        const split = sh.splitter(seq);
        const maxSeqLen = split ? split.length : 30;

        // makes no sense to have scroller if we have shorter than 50 positions
        if (maxSeqLen < 50)
          continue;
        let conservationData: number[] = [];
        if (cats.length > 1) conservationData = calculateConservation(cats.filter(Boolean), sh.splitter);

        const defaultHeaderHeight = 40;
        const scroller = new MSAScrollingHeader({
          canvas: grid.overlay,
          headerHeight: defaultHeaderHeight,
          totalPositions: maxSeqLen + 1,
          conservationData: conservationData, // Pass conservation data directly
          conservationHeight: 40, // Set height for conservation barplot
          onPositionChange: (scrollerCur, scrollerRange) => {
            setTimeout(() => {
              const start = getStart();
              const cur = getCurrent();
              if (start !== scrollerRange.start)
                seqCol.setTag(bioTAGS.positionShift, (scrollerRange.start - 1).toString());
              if (cur !== scrollerCur) {
                seqCol.setTag(bioTAGS.selectedPosition, (scrollerCur).toString());
                if (scrollerCur >= 0 && !positionStatsViewerAddedOnce && grid.tableView) {
                  positionStatsViewerAddedOnce = true;
                  const v = grid.tableView.addViewer('Sequence Position Statistics', { sequenceColumnName: seqCol.name });
                  grid.tableView.dockManager.dock(v, DG.DOCK_TYPE.DOWN, null, 'Sequence Position Statistics', 0.4);
                }
              }
            });
          },
        });

        grid.props.colHeaderHeight = 65;
        setTimeout(() => { if (grid.isDetached) return; gCol.width = 400; }, 300); // needed because renderer sets its width
        grid.sub(grid.onCellRender.subscribe((e) => {
          const cell = e.cell;
          if (!cell || !cell.isColHeader || cell?.gridColumn?.name !== gCol?.name)
            return;
          const cellBounds = e.bounds;
          if (!cellBounds || cellBounds.height <= 50) {
            scroller.headerHeight = 0;
            return;
          }

          const headerTitleSpace = 20;
          const availableTitleSpace = cellBounds.height - defaultHeaderHeight;
          // Only draw title if we have enough space for it
          if (availableTitleSpace > headerTitleSpace) {
            // update header height to fit whole header
            scroller.headerHeight = cellBounds.height - headerTitleSpace;

            const ctx = grid.overlay.getContext('2d');
            if (!ctx)
              return;
            // Save context state
            ctx.save();
            ctx.rect(cellBounds.x, cellBounds.y, cellBounds.width, headerTitleSpace);
            ctx.clip();
            // Draw title text
            ctx.font = grid.props.colHeaderFont ?? 'bold 13px Roboto, Roboto Local';
            ctx.fillStyle = '#4a4a49';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';

            // Position the title in the middle of the available space above the header
            const titleX = cellBounds.x + cellBounds.width / 2;
            const titleY = cellBounds.y + headerTitleSpace / 2;

            // Draw the text
            ctx.fillText(seqCol.name ?? '', titleX, titleY);

            ctx.restore();
          } else
            scroller.headerHeight = Math.max(defaultHeaderHeight, cellBounds.height);

          const titleShift = Math.max(0, availableTitleSpace ?? 0) > headerTitleSpace ? headerTitleSpace : 0;
          const font = getFontSize();
          scroller.positionWidth = font + (isMSA ? 8 : 0);
          const start = getStart();
          const startPadding = isMSA ? 0 : 4;
          scroller.draw(cellBounds.x + startPadding, cellBounds.y + titleShift, cellBounds.width - startPadding, cellBounds.height - titleShift, getCurrent(), start, e);
        }));
      }
    }, 1000);
  };
  // handle all new grids
  const _ = grok.events.onViewerAdded.subscribe((e) => {
    if (!e.args || !(e.args.viewer instanceof DG.Grid))
      return;
    const grid = e.args.viewer as DG.Grid;
    handleGrid(grid);
  });

  const openTables = grok.shell.tableViews;
  for (const tv of openTables) {
    const grid = tv?.grid;
    if (grid)
      handleGrid(grid);
  }
}
