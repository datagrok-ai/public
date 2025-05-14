/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {MSAScrollingHeader} from '@datagrok-libraries/bio/src/utils/sequence-position-scroller';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {ALPHABET, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_package} from '../package';

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
        if (sh.isHelm() || sh.alphabet === ALPHABET.UN || !sh.isMsa())
          continue;

        const gCol = grid.col(seqCol.name);
        if (!gCol)
          continue;

        let positionStatsViewerAddedOnce = false;

        const ifNan = (a: number, els: number) => (Number.isNaN(a) ? els : a);
        const getStart = () => ifNan(Math.max(Number.parseInt(seqCol.getTag(bioTAGS.positionShift) ?? '0'), 0), 0) + 1;
        const getCurrent = () => ifNan(Number.parseInt(seqCol.getTag(bioTAGS.selectedPosition) ?? '-2'), -2);
        const getFontSize = () => MonomerPlacer.getFontSettings(seqCol).fontWidth;
        // get the maximum length of sequences by randomly taking 10 sequences;
        let maxSeqLen = 0;
        for (let i = 0; i < Math.min(10, df.rowCount); i++) {
          const row = Math.floor(Math.random() * df.rowCount - 1);
          const seq = sh.getSplitted(row);
          if (seq)
            maxSeqLen = Math.max(maxSeqLen, seq.length);
        }

        // makes no sense to have scroller if we have shorter than 50 positions
        if (maxSeqLen < 50)
          continue;

        const defaultHeaderHeight = 40;
        const scroller = new MSAScrollingHeader({
          canvas: grid.overlay,
          headerHeight: defaultHeaderHeight,
          totalPositions: maxSeqLen + 1,
          onPositionChange: (scrollerCur, scrollerRange) => {
            setTimeout(() => {
              const start = getStart();
              const cur = getCurrent();
              if (start !== scrollerRange.start)
                seqCol.setTag(bioTAGS.positionShift, (scrollerRange.start - 1).toString());
              if (scrollerCur >= 0 && cur !== scrollerCur) {
                seqCol.setTag(bioTAGS.selectedPosition, (scrollerCur).toString());
                if (!positionStatsViewerAddedOnce && grid.tableView) {
                  positionStatsViewerAddedOnce = true;
                  const v = grid.tableView.addViewer('Sequence Position Statistics', {sequenceColumnName: seqCol.name});
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
          if (!cellBounds || cellBounds.height <= 50)
            return;

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
          scroller.positionWidth = font + 8;
          const start = getStart();
          //const positionXShift = start > 1 ? font + 8 : 0;
          // pass in the event to the scroller so it can internally preventDefault if all is well
          scroller.draw(cellBounds.x, cellBounds.y + titleShift, cellBounds.width, cellBounds.height - titleShift, getCurrent(), start, e);
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
