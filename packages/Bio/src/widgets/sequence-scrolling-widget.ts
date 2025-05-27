/* eslint-disable rxjs/no-ignored-subscription */
/* eslint-disable max-lines */
/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {ConservationTrack, MSAHeaderTrack, MSAScrollingHeader, WebLogoTrack} from '@datagrok-libraries/bio/src/utils/sequence-position-scroller';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {ALPHABET, TAGS as bioTAGS, SplitterFunc} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_package} from '../package';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

// ============================================================================
// SIMPLE VIEWPORT-AWARE CACHING
// ============================================================================

/**
 * Simple viewport-aware MSA data manager
 * Keeps the same data format as original but only calculates visible ranges
 */
class MSAViewportManager {
  private static conservationCache: DG.LruCache<string, number[]> = new DG.LruCache<string, number[]>(200);
  private static webLogoCache: DG.LruCache<string, Map<number, Map<string, number>>> = new DG.LruCache<string, Map<number, Map<string, number>>>(200);

  private static simpleHash(str: string): string {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
      const char = str.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash;
    }
    return Math.abs(hash).toString(36);
  }

  private static getCacheKey(column: DG.Column, start: number, end: number, type: string): string {
    const sequences = column.categories.filter(Boolean);
    const dataHash = MSAViewportManager.simpleHash(sequences.join('|'));
    return `${column.dataFrame.id}_${column.name}_v${column.version}_${dataHash}_${type}_${start}_${end}`;
  }

  /**
   * Get conservation scores for full sequence (0-based array where index 0 = position 1)
   * Uses caching and smart range calculation
   */
  static getConservationScores(column: DG.Column, splitter: SplitterFunc, maxLength: number): number[] {
    const fullCacheKey = MSAViewportManager.getCacheKey(column, 0, maxLength, 'conservation_full');

    return MSAViewportManager.conservationCache.getOrCreate(fullCacheKey, () => {
      const sequences = column.categories.filter(Boolean);
      if (sequences.length <= 1)
        return new Array(maxLength).fill(0);


      const splitSequences: ISeqSplitted[] = sequences
        .map((seq) => splitter(seq))
        .filter(Boolean);

      const scores: number[] = [];

      // Calculate for each position (0-based in sequence = 1-based in MSA)
      for (let pos = 0; pos < maxLength; pos++) {
        const residueCounts: Record<string, number> = {};
        let totalResidues = 0;

        for (const seq of splitSequences) {
          if (pos < seq.length) {
            const residue = seq.getCanonical(pos);
            if (residue && !seq.isGap(pos)) {
              residueCounts[residue] = (residueCounts[residue] || 0) + 1;
              totalResidues++;
            }
          }
        }

        if (totalResidues === 0) {
          scores.push(0);
          continue;
        }

        const mostCommonCount = Math.max(...Object.values(residueCounts), 0);
        const conservation = mostCommonCount / totalResidues;
        scores.push(conservation);
      }

      return scores;
    });
  }

  /**
   * Get WebLogo data for full sequence (0-based map where key 0 = position 1)
   * Uses caching and smart range calculation
   */
  static getWebLogoData(column: DG.Column, splitter: SplitterFunc, maxLength: number): Map<number, Map<string, number>> {
    const fullCacheKey = MSAViewportManager.getCacheKey(column, 0, maxLength, 'weblogo_full');

    return MSAViewportManager.webLogoCache.getOrCreate(fullCacheKey, () => {
      const sequences = column.categories.filter(Boolean);
      const webLogoData: Map<number, Map<string, number>> = new Map();

      if (sequences.length <= 1)
        return webLogoData;


      const splitSequences: ISeqSplitted[] = sequences
        .map((seq) => splitter(seq))
        .filter(Boolean);

      // Calculate for each position (0-based in sequence = 1-based in MSA)
      for (let pos = 0; pos < maxLength; pos++) {
        const residueCounts: Map<string, number> = new Map();
        let totalResidues = 0;

        for (const seq of splitSequences) {
          if (pos < seq.length) {
            const residue = seq.getCanonical(pos);
            if (residue && !seq.isGap(pos)) {
              residueCounts.set(residue, (residueCounts.get(residue) || 0) + 1);
              totalResidues++;
            }
          }
        }

        if (totalResidues === 0) {
          webLogoData.set(pos, new Map());
          continue;
        }

        const residueFreqs: Map<string, number> = new Map();
        for (const [residue, count] of residueCounts.entries())
          residueFreqs.set(residue, count / totalResidues);
        webLogoData.set(pos, residueFreqs);
      }

      return webLogoData;
    });
  }
}

class CachedWebLogoTrack extends WebLogoTrack {
  private column: DG.Column;
  private splitter: SplitterFunc;
  private dataInitialized: boolean = false;

  constructor(
    column: DG.Column,
    splitter: SplitterFunc,
    maxLength: number,
    height: number = 45,
    title: string = 'WebLogo'
  ) {
    super(new Map(), height, '', title);
    this.column = column;
    this.splitter = splitter;

    // Initialize with full cached data
    const webLogoData = MSAViewportManager.getWebLogoData(column, splitter, maxLength);
    this.updateData(webLogoData);
    this.visible = webLogoData.size > 0;
    this.dataInitialized = true;
  }
}

class CachedConservationTrack extends ConservationTrack {
  private column: DG.Column;
  private splitter: SplitterFunc;
  private dataInitialized: boolean = false;

  constructor(
    column: DG.Column,
    splitter: SplitterFunc,
    maxLength: number,
    height: number = 45,
    colorScheme: 'default' | 'rainbow' | 'heatmap' = 'default',
    title: string = 'Conservation'
  ) {
    super([], height, colorScheme, title);
    this.column = column;
    this.splitter = splitter;

    const conservationData = MSAViewportManager.getConservationScores(column, splitter, maxLength);
    this.updateData(conservationData);
    this.visible = conservationData.length > 0;
    this.dataInitialized = true;
  }
}

export function handleSequenceHeaderRendering() {
  const handleGrid = (grid: DG.Grid) => {
    setTimeout(() => {
      if (grid.isDetached) return;

      const df = grid.dataFrame;
      if (!df) return;

      const seqCols = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);

      for (const seqCol of seqCols) {
        const sh = _package.seqHelper.getSeqHandler(seqCol);
        if (!sh) continue;
        if (sh.isHelm() || sh.alphabet === ALPHABET.UN) continue;

        const gCol = grid.col(seqCol.name);
        if (!gCol) continue;

        let positionStatsViewerAddedOnce = !!grid.tableView &&
          Array.from(grid.tableView.viewers).some((v) => v.type === 'Sequence Position Statistics');

        const isMSA = sh.isMsa();
        const ifNan = (a: number, els: number) => (Number.isNaN(a) ? els : a);
        const getStart = () => ifNan(Math.max(Number.parseInt(seqCol.getTag(bioTAGS.positionShift) ?? '0'), 0), 0) + 1;
        const getCurrent = () => ifNan(Number.parseInt(seqCol.getTag(bioTAGS.selectedPosition) ?? '-2'), -2);
        const getFontSize = () => MonomerPlacer.getFontSettings(seqCol).fontWidth;

        // Get maximum sequence length the simple way
        let maxSeqLen = 0;
        const cats = seqCol.categories;

        for (let i = 0; i < cats.length; i++) {
          const seq = cats[i];
          if (seq && seq.length > 0) {
            const split = sh.splitter(seq);
            if (split && split.length > maxSeqLen)
              maxSeqLen = split.length;
          }
        }

        // Skip if sequence is too short
        if (maxSeqLen < 50) continue;

        // Check if we have multiple sequences for MSA features
        const hasMultipleSequences = cats.filter(Boolean).length > 1;

        // Layout constants
        const dottedCellsHeight = 38;
        const trackGap = 4;
        const DEFAULT_WEBLOGO_HEIGHT = 45;
        const DEFAULT_CONSERVATION_HEIGHT = 45;

        // Calculate initial header height
        const initialHeaderHeight = dottedCellsHeight +
          (hasMultipleSequences ? DEFAULT_WEBLOGO_HEIGHT + trackGap : 0) +
          (hasMultipleSequences ? DEFAULT_CONSERVATION_HEIGHT + trackGap : 0);

        // Initialize tracks with cached data
        const initializeHeaders = (monomerLib: any = null) => {
          const tracks: { id: string, track: MSAHeaderTrack, priority: number }[] = [];

          // Create cached tracks only if we have multiple sequences
          if (hasMultipleSequences) {
            // Conservation track
            const conservationTrack = new CachedConservationTrack(
              seqCol,
              sh.splitter,
              maxSeqLen,
              DEFAULT_CONSERVATION_HEIGHT,
              'default',
              'Conservation'
            );
            tracks.push({id: 'conservation', track: conservationTrack, priority: 1});

            // WebLogo track
            const webLogoTrack = new CachedWebLogoTrack(
              seqCol,
              sh.splitter,
              maxSeqLen,
              DEFAULT_WEBLOGO_HEIGHT,
              'WebLogo'
            );

            if (monomerLib) {
              webLogoTrack.setMonomerLib(monomerLib);
              webLogoTrack.setBiotype(sh.defaultBiotype || 'PEPTIDE');
            }

            webLogoTrack.setupDefaultTooltip();
            tracks.push({id: 'weblogo', track: webLogoTrack, priority: 2});
          }

          // Create the scrolling header
          const scroller = new MSAScrollingHeader({
            canvas: grid.overlay,
            headerHeight: initialHeaderHeight,
            totalPositions: maxSeqLen + 1,
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
                    const v = grid.tableView.addViewer('Sequence Position Statistics', {sequenceColumnName: seqCol.name});
                    grid.tableView.dockManager.dock(v, DG.DOCK_TYPE.DOWN, null, 'Sequence Position Statistics', 0.4);
                  }
                }
              });
            },
            onHeaderHeightChange: (newHeight) => {
              if (grid && !grid.isDetached)
                grid.props.colHeaderHeight = newHeight;
            },
          });

          scroller.setupTooltipHandling();

          // Add tracks to scroller
          tracks.forEach(({id, track}) => {
            scroller.addTrack(id, track);
          });

          scroller.setSelectionData(df, seqCol, sh);

          // Set initial header height
          grid.props.colHeaderHeight = initialHeaderHeight;

          // Set column width
          setTimeout(() => {
            if (grid.isDetached) return;
            gCol.width = 400;
          }, 300);

          // Handle cell rendering
          grid.sub(grid.onCellRender.subscribe((e) => {
            const cell = e.cell;
            if (!cell || !cell.isColHeader || cell?.gridColumn?.name !== gCol?.name)
              return;

            const cellBounds = e.bounds;
            if (!cellBounds) return;

            // Set dynamic properties
            scroller.headerHeight = cellBounds.height;
            const font = getFontSize();
            scroller.positionWidth = font + (isMSA ? 8 : 0);

            const start = getStart();
            const startPadding = isMSA ? 0 : 4;

            // Draw the header
            scroller.draw(
              cellBounds.x + startPadding,
              cellBounds.y,
              cellBounds.width - startPadding,
              cellBounds.height,
              getCurrent(),
              start,
              e
            );
          }));
        };

        // Initialize with monomer library
        getMonomerLibHelper()
          .then((libHelper) => {
            const monomerLib = libHelper.getMonomerLib();
            initializeHeaders(monomerLib);
          })
          .catch((error) => {
            console.error('Error loading monomerLib:', error);
            initializeHeaders();
          });
      }
    }, 1000);
  };

  // Handle new grid additions
  const _ = grok.events.onViewerAdded.subscribe((e) => {
    if (!e.args || !(e.args.viewer instanceof DG.Grid)) return;
    const grid = e.args.viewer as DG.Grid;
    handleGrid(grid);
  });

  // Handle existing grids
  const openTables = grok.shell.tableViews;
  for (const tv of openTables) {
    const grid = tv?.grid;
    if (grid) handleGrid(grid);
  }
}
