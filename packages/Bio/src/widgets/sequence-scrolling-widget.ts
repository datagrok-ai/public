/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import { MSAScrollingHeader } from '@datagrok-libraries/bio/src/utils/sequence-position-scroller';
import { MonomerPlacer } from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import { ALPHABET, TAGS as bioTAGS } from '@datagrok-libraries/bio/src/utils/macromolecule';
import { _package } from '../package';
import { ISeqSplitted, SplitterFunc } from '@datagrok-libraries/bio/src/utils/macromolecule/types';


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
        // (handles modifications, non-standard residues, etc.)
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

        // Create a proper ILogger implementation

        let positionStatsViewerAddedOnce = !!grid.tableView && Array.from(grid.tableView.viewers).some((v) => v.type === 'Sequence Position Statistics');
        const isMSA = sh.isMsa();
        const ifNan = (a: number, els: number) => (Number.isNaN(a) ? els : a);
        const getStart = () => ifNan(Math.max(Number.parseInt(seqCol.getTag(bioTAGS.positionShift) ?? '0'), 0), 0) + 1;
        const getCurrent = () => ifNan(Number.parseInt(seqCol.getTag(bioTAGS.selectedPosition) ?? '-2'), -2);
        const getFontSize = () => MonomerPlacer.getFontSettings(seqCol).fontWidth;

        // Provide props for MonomerPlacer
        const fontSettings = MonomerPlacer.getFontSettings(seqCol);
        const propsProvider = () => {
          const { font, fontWidth } = fontSettings;
          return {
            separatorWidth: sh.isMsa() ? 8 : 5,
            monomerToShort: (mon:any, limit:any) => mon.length > limit ? mon.substring(0, limit) : mon,
            font,
            fontCharWidth: fontWidth
          };
        };

        // Create and initialize MonomerPlacer
        const monomerPlacer = new MonomerPlacer(gCol, seqCol, _package.logger, 50, propsProvider);
        monomerPlacer.init();  // Initialize it

        // Override the onClick method with our fixed version
        monomerPlacer.onClick = function (gridCell: DG.GridCell, e: MouseEvent): void {
          // First check if seqHelper is initialized
          if (!_package.seqHelper || gridCell.tableRowIndex == null) return;

          // Calculate position with the standard method
          const positionShift = this.positionShift;
          const gridCellBounds = gridCell.bounds;
          const argsX = e.offsetX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCellBounds.x);
          const position = this.getPosition(gridCell.tableRowIndex!, argsX, gridCellBounds.width);

          if (position !== null && position >= 0) {
            // Add 1 to fix the off-by-one error
            const correctedPosition = position + 1;

            // Update the selected position with correction
            this.tableCol.setTag(bioTAGS.selectedPosition, (correctedPosition + positionShift).toString());

            // Handle stats viewer
            if (correctedPosition + positionShift >= 0 && !positionStatsViewerAddedOnce && grid.tableView) {
              positionStatsViewerAddedOnce = true;
              const v = grid.tableView.addViewer('Sequence Position Statistics', { sequenceColumnName: seqCol.name });
              grid.tableView.dockManager.dock(v, DG.DOCK_TYPE.DOWN, null, 'Sequence Position Statistics', 0.4);
            }

            // Check if we should center the view
            const font = getFontSize();
            const charWidth = font + (isMSA ? 8 : 0);
            const visiblePositions = Math.floor(gridCellBounds.width / charWidth);
            const start = getStart();

            if ((correctedPosition + positionShift) >= 1 &&
              ((correctedPosition + positionShift) < start ||
                (correctedPosition + positionShift) > start + visiblePositions - 1)) {
              // Center view on clicked position if possible
              const newStart = Math.max(1, (correctedPosition + positionShift) - Math.floor(visiblePositions / 2));
              this.tableCol.setTag(bioTAGS.positionShift, (newStart - 1).toString());
            }
          }
        };

        // Store it for later access
        seqCol.temp.monomerPlacer = monomerPlacer;

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

        const defaultHeaderHeight = 40;
        let conservationData: number[] = [];
        if (isMSA && cats.length > 1) {
          // Get all sequences and calculate conservation
          const sequences = cats.filter(Boolean);
          conservationData = calculateConservation(sequences, sh.splitter);
        }

        const positionMarkersHeight = 30;
        const sliderHeight = 8;
        const padding = 5;
        const conservationHeight = 40;

        // Minimum height needed for basic header (positions + slider)
        const minHeaderHeight = positionMarkersHeight + sliderHeight + padding;

        // Calculate header height based on current grid header height
        // This ensures we start with just the basic header
        const headerTitleSpace = 20;
        const initialHeaderHeight = minHeaderHeight;


        const scroller = new MSAScrollingHeader({
          canvas: grid.overlay,
          headerHeight: initialHeaderHeight,
          totalPositions: maxSeqLen + 1,
          conservationData: conservationData as number[],
          conservationHeight: conservationHeight,
          sliderHeight: sliderHeight,
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
        grid.props.colHeaderHeight = minHeaderHeight + headerTitleSpace;

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
          const availableTitleSpace = cellBounds.height - minHeaderHeight;

          // Only draw title if we have enough space for it
          if (availableTitleSpace > headerTitleSpace) {
            // Available height for header components (excluding title)
            const availableComponentHeight = cellBounds.height - headerTitleSpace;

            // Now determine how much of this height to use based on thresholds
            let headerHeight;

            const minHeightForConservation = minHeaderHeight + conservationHeight;

            if (availableComponentHeight >= minHeightForConservation) {
              // Space for everything including conservation
              headerHeight = availableComponentHeight;
            } else {
              // Only space for basic header (no conservation)
              headerHeight = minHeaderHeight;
            }

            scroller.headerHeight = headerHeight;

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

            const titleX = cellBounds.x + cellBounds.width / 2;
            const titleY = cellBounds.y + headerTitleSpace / 2;

            ctx.fillText(seqCol.name ?? '', titleX, titleY);

            ctx.restore();
          } else {
            // Not enough space for title - just use available height
            scroller.headerHeight = Math.max(minHeaderHeight, cellBounds.height);
          }

          const titleShift = Math.max(0, availableTitleSpace ?? 0) > headerTitleSpace ? headerTitleSpace : 0;
          const font = getFontSize();
          scroller.positionWidth = font + (isMSA ? 8 : 0);
          const start = getStart();
          const startPadding = isMSA ? 0 : 4;
          scroller.draw(cellBounds.x + startPadding, cellBounds.y + titleShift, cellBounds.width - startPadding, cellBounds.height - titleShift, getCurrent(), start, e);
        }));


        grid.root.addEventListener('click', (e) => {
          // Check if grid is still valid
          if (grid.isDetached) return;

          // Get click coordinates relative to the canvas
          const rect = grid.canvas.getBoundingClientRect();
          const x = e.clientX - rect.left;
          const y = e.clientY - rect.top;

          // Get cell at click position

          // const cell = grid.cellAtCoordinates(x, y);
          const cell = grid.hitTest(x, y);
          if (!cell || cell.isColHeader) return;

          // Check if click is in our sequence column
          const clickedColName = cell?.gridColumn?.name;
          if (clickedColName !== seqCol.name) return;

          // Get position in the sequence
          // Get row index
          const rowIdx = cell.tableRowIndex;
          if (rowIdx === null) return;

          // Get cell bounds
          const bounds = cell.bounds;
          const positionX = x - bounds.x;

          // Calculate clicked position similar to MonomerPlacer's logic
          const start = getStart();
          const font = getFontSize();
          const charWidth = font + 8; // Same as positionWidth

          // Calculate character index
          const charIdx = Math.floor(positionX / charWidth);
          const clickedPosition = start + charIdx;

          console.log(`Clicked position ${clickedPosition} (charIdx=${charIdx}, start=${start})`);

          // Update the selected position if valid
          if (clickedPosition >= 1 && clickedPosition <= maxSeqLen) {
            seqCol.setTag(bioTAGS.selectedPosition, clickedPosition.toString());

            // Center view if needed
            const visiblePositions = Math.floor(bounds.width / charWidth);
            if (clickedPosition < start || clickedPosition > start + visiblePositions - 1) {
              // Center view on clicked position if possible
              const newStart = Math.max(1, clickedPosition - Math.floor(visiblePositions / 2));
              seqCol.setTag(bioTAGS.positionShift, (newStart - 1).toString());
            }
          }
        });

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

