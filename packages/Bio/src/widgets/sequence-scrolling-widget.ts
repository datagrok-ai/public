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
// OPTIMIZED VIEWPORT-AWARE CACHING WITH FORCE UPDATE SUPPORT
// ============================================================================
class MSAViewportManager {
  private static conservationCache: DG.LruCache<string, number[]> = new DG.LruCache<string, number[]>(100);
  private static webLogoCache: DG.LruCache<string, Map<number, Map<string, number>>> = new DG.LruCache<string, Map<number, Map<string, number>>>(100);

  // Track when data was last invalidated for force updates
  private static lastInvalidationTime: number = 0;

  // Cache chunks of 200 positions at a time
  private static readonly CHUNK_SIZE = 200;

  /**
   * Clear all caches and update invalidation timestamp for force updates
   */
  static clearAllCaches(): void {
    MSAViewportManager.conservationCache = new DG.LruCache<string, number[]>(100);
    MSAViewportManager.webLogoCache = new DG.LruCache<string, Map<number, Map<string, number>>>(100);
    MSAViewportManager.lastInvalidationTime = Date.now();
  }

  /**
   * Get the last invalidation time (used by tracks to detect forced updates)
   */
  static getLastInvalidationTime(): number {
    return MSAViewportManager.lastInvalidationTime;
  }

  private static getChunkCacheKey(column: DG.Column, chunkStart: number, chunkEnd: number, type: string): string {
    const dataFrame = column.dataFrame;
    // Use framework's built-in versioning instead of complex hash
    return `${dataFrame.id}_${column.name}_v${column.version}_f${dataFrame.filter.version}_${type}_${chunkStart}_${chunkEnd}`;
  }

  private static getConservationChunk(seqHandler: any, start: number, end: number, cacheKey: string): number[] {
    return MSAViewportManager.conservationCache.getOrCreate(cacheKey, () => {
      const dataFrame = seqHandler.column.dataFrame;
      const chunkSize = end - start;

      // Pre-allocate result array
      const scores = new Array(chunkSize);

      const visibleCount = dataFrame.filter.trueCount;
      if (visibleCount <= 1) {
        scores.fill(0);
        return scores;
      }

      // Pre-allocate sequences array with known size
      const sequences = new Array(visibleCount);
      let seqIndex = 0;

      const filter = dataFrame.filter;
      for (let i = -1; (i = filter.findNext(i, true)) !== -1;)
        sequences[seqIndex++] = seqHandler.getSplitted(i);


      // Calculate conservation for the chunk
      for (let pos = start; pos < end; pos++) {
        const residueCounts: Record<string, number> = {};
        let totalResidues = 0;
        const relativePos = pos - start;

        for (let seqIdx = 0; seqIdx < sequences.length; seqIdx++) {
          const seq = sequences[seqIdx];
          if (pos < seq.length) {
            const residue = seq.getCanonical(pos);
            if (residue && !seq.isGap(pos)) {
              residueCounts[residue] = (residueCounts[residue] || 0) + 1;
              totalResidues++;
            }
          }
        }

        if (totalResidues === 0) {
          scores[relativePos] = 0;
          continue;
        }

        const mostCommonCount = Math.max(...Object.values(residueCounts), 0);
        scores[relativePos] = mostCommonCount / totalResidues;
      }

      return scores;
    });
  }

  private static getWebLogoChunk(seqHandler: any, start: number, end: number, cacheKey: string): Map<number, Map<string, number>> {
    return MSAViewportManager.webLogoCache.getOrCreate(cacheKey, () => {
      const dataFrame = seqHandler.column.dataFrame;
      const webLogoData: Map<number, Map<string, number>> = new Map();

      // OPTIMIZED: Use filter iterator
      const visibleCount = dataFrame.filter.trueCount;
      if (visibleCount <= 1)
        return webLogoData;


      // Pre-allocate sequences array
      const sequences = new Array(visibleCount);
      let seqIndex = 0;

      const filter = dataFrame.filter;
      for (let i = -1; (i = filter.findNext(i, true)) !== -1;) {
        // OPTIMIZED: Use seqHandler.getSplitted which has internal caching
        sequences[seqIndex++] = seqHandler.getSplitted(i);
      }

      // Calculate WebLogo data for the chunk
      for (let pos = start; pos < end; pos++) {
        const residueCounts: Map<string, number> = new Map();
        let totalResidues = 0;

        for (let seqIdx = 0; seqIdx < sequences.length; seqIdx++) {
          const seq = sequences[seqIdx];
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

  static getConservationForViewport(seqHandler: any, viewportStart: number, viewportEnd: number, maxLength: number): number[] {
    // Create a full-sized array filled with zeros
    const result: number[] = new Array(maxLength).fill(0);

    // Calculate which chunks we need
    const startChunk = Math.floor(viewportStart / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;
    const endChunk = Math.ceil(viewportEnd / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;

    // Get all needed chunks and place them in the correct positions
    for (let chunkStart = startChunk; chunkStart < endChunk; chunkStart += MSAViewportManager.CHUNK_SIZE) {
      const chunkEnd = Math.min(chunkStart + MSAViewportManager.CHUNK_SIZE, maxLength);
      const cacheKey = MSAViewportManager.getChunkCacheKey(seqHandler.column, chunkStart, chunkEnd, 'conservation');
      const chunkData = MSAViewportManager.getConservationChunk(seqHandler, chunkStart, chunkEnd, cacheKey);

      // Copy chunk data to the correct positions in the result array
      for (let i = 0; i < chunkData.length && (chunkStart + i) < maxLength; i++)
        result[chunkStart + i] = chunkData[i];
    }

    return result;
  }

  static getWebLogoForViewport(seqHandler: any, viewportStart: number, viewportEnd: number, maxLength: number): Map<number, Map<string, number>> {
    const result: Map<number, Map<string, number>> = new Map();

    // Calculate which chunks we need
    const startChunk = Math.floor(viewportStart / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;
    const endChunk = Math.ceil(viewportEnd / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;

    // Get all needed chunks and place them in the correct positions
    for (let chunkStart = startChunk; chunkStart < endChunk; chunkStart += MSAViewportManager.CHUNK_SIZE) {
      const chunkEnd = Math.min(chunkStart + MSAViewportManager.CHUNK_SIZE, maxLength);
      const cacheKey = MSAViewportManager.getChunkCacheKey(seqHandler.column, chunkStart, chunkEnd, 'weblogo');
      const chunkData = MSAViewportManager.getWebLogoChunk(seqHandler, chunkStart, chunkEnd, cacheKey);

      // Copy chunk data to the result map (positions are already correct)
      for (const [pos, data] of chunkData.entries()) {
        if (pos < maxLength)
          result.set(pos, data);
      }
    }

    return result;
  }
}

// ============================================================================
//  SMART TRACKS
// ============================================================================

class LazyWebLogoTrack extends WebLogoTrack {
  private seqHandler: any;
  private maxLength: number;
  private lastViewportStart: number = -1;
  private lastViewportEnd: number = -1;
  private lastInvalidationTime: number = 0;
  private forceNextUpdate: boolean = false;

  constructor(
    seqHandler: any,
    maxLength: number,
    height: number = 45,
    title: string = 'WebLogo'
  ) {
    super(new Map(), height, '', title);
    this.seqHandler = seqHandler;
    this.maxLength = maxLength;

    this.visible = seqHandler.column.dataFrame.filter.trueCount > 1;
  }

  public forceUpdate(): void {
    this.forceNextUpdate = true;
  }

  public resetViewportTracking(): void {
    this.lastViewportStart = -1;
    this.lastViewportEnd = -1;
  }

  private updateForViewport(viewportStart: number, viewportEnd: number): void {
    const currentInvalidationTime = MSAViewportManager.getLastInvalidationTime();
    const cacheWasInvalidated = currentInvalidationTime > this.lastInvalidationTime;

    const viewportChanged = (
      Math.abs(this.lastViewportStart - viewportStart) >= 10 ||
      Math.abs(this.lastViewportEnd - viewportEnd) >= 10
    );

    // Only update if viewport changed significantly OR cache was invalidated OR force update is set
    if (!viewportChanged && !cacheWasInvalidated && !this.forceNextUpdate)
      return;


    // Reset flags and update tracking
    this.forceNextUpdate = false;
    this.lastViewportStart = viewportStart;
    this.lastViewportEnd = viewportEnd;
    this.lastInvalidationTime = currentInvalidationTime;

    // Add buffer to reduce cache misses on small scrolls
    const bufferedStart = Math.max(0, viewportStart - 50);
    const bufferedEnd = Math.min(this.maxLength, viewportEnd + 50);

    const webLogoData = MSAViewportManager.getWebLogoForViewport(
      this.seqHandler,
      bufferedStart,
      bufferedEnd,
      this.maxLength
    );

    this.updateData(webLogoData);
  }

  draw(x: number, y: number, width: number, height: number, windowStart: number,
    positionWidth: number, totalPositions: number, currentPosition: number): void {
    // Calculate what positions are actually visible
    const visiblePositions = Math.ceil(width / positionWidth) + 2;
    const viewportEnd = Math.min(windowStart + visiblePositions, totalPositions);

    // Update data for just the visible range
    this.updateForViewport(windowStart - 1, viewportEnd - 1); // Convert to 0-based

    // Call parent draw method
    super.draw(x, y, width, height, windowStart, positionWidth, totalPositions, currentPosition);
  }
}

class LazyConservationTrack extends ConservationTrack {
  private seqHandler: any;
  private maxLength: number;
  private lastViewportStart: number = -1;
  private lastViewportEnd: number = -1;
  private lastInvalidationTime: number = 0;
  private forceNextUpdate: boolean = false;

  constructor(
    seqHandler: any,
    maxLength: number,
    height: number = 45,
    colorScheme: 'default' | 'rainbow' | 'heatmap' = 'default',
    title: string = 'Conservation'
  ) {
    super([], height, colorScheme, title);
    this.seqHandler = seqHandler;
    this.maxLength = maxLength;

    // OPTIMIZED: Use dataFrame.filter.trueCount directly
    this.visible = seqHandler.column.dataFrame.filter.trueCount > 1;
  }

  public forceUpdate(): void {
    this.forceNextUpdate = true;
  }

  public resetViewportTracking(): void {
    this.lastViewportStart = -1;
    this.lastViewportEnd = -1;
  }

  private updateForViewport(viewportStart: number, viewportEnd: number): void {
    const currentInvalidationTime = MSAViewportManager.getLastInvalidationTime();
    const cacheWasInvalidated = currentInvalidationTime > this.lastInvalidationTime;

    const viewportChanged = (
      Math.abs(this.lastViewportStart - viewportStart) >= 10 ||
      Math.abs(this.lastViewportEnd - viewportEnd) >= 10
    );

    // Only update if viewport changed significantly OR cache was invalidated OR force update is set
    if (!viewportChanged && !cacheWasInvalidated && !this.forceNextUpdate)
      return;


    // Reset flags and update tracking
    this.forceNextUpdate = false;
    this.lastViewportStart = viewportStart;
    this.lastViewportEnd = viewportEnd;
    this.lastInvalidationTime = currentInvalidationTime;

    // Add buffer to reduce cache misses on small scrolls
    const bufferedStart = Math.max(0, viewportStart - 50);
    const bufferedEnd = Math.min(this.maxLength, viewportEnd + 50);

    const conservationScores = MSAViewportManager.getConservationForViewport(
      this.seqHandler,
      bufferedStart,
      bufferedEnd,
      this.maxLength
    );

    this.updateData(conservationScores);
  }

  draw(x: number, y: number, width: number, height: number, windowStart: number,
    positionWidth: number, totalPositions: number, currentPosition: number): void {
    // Calculate what positions are actually visible
    const visiblePositions = Math.ceil(width / positionWidth) + 2;
    const viewportEnd = Math.min(windowStart + visiblePositions, totalPositions);

    this.updateForViewport(windowStart - 1, viewportEnd - 1); // Convert to 0-based

    // Call parent draw method
    super.draw(x, y, width, height, windowStart, positionWidth, totalPositions, currentPosition);
  }
}

// ============================================================================
//  MAIN HANDLER
// ============================================================================

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

        // Get maximum sequence length
        let maxSeqLen = 0;
        if (sh.maxLength)
          maxSeqLen = sh.maxLength;
        else {
          // Fallback to manual calculation if needed
          const cats = seqCol.categories;
          for (let i = 0; i < cats.length; i++) {
            const seq = cats[i];
            if (seq && seq.length > 0) {
              const split = sh.splitter(seq);
              if (split && split.length > maxSeqLen)
                maxSeqLen = split.length;
            }
          }
        }

        // Skip if sequence is too short
        if (maxSeqLen < 50) continue;

        const getHasMultipleSequences = () => df.filter.trueCount > 1;
        let hasMultipleSequences = getHasMultipleSequences();

        const STRICT_THRESHOLDS = {
          BASE: 38, // DOTTED_CELL_HEIGHT(30) + SLIDER_HEIGHT(8)
          WITH_TITLE: 58, // BASE + TITLE_HEIGHT(16) + TRACK_GAP(4)
          WITH_WEBLOGO: 107, // WITH_TITLE + DEFAULT_TRACK_HEIGHT(45) + TRACK_GAP(4)
          WITH_BOTH: 156 // WITH_WEBLOGO + DEFAULT_TRACK_HEIGHT(45) + TRACK_GAP(4)
        };

        let initialHeaderHeight: number;
        if (!hasMultipleSequences) {
          // Single sequence: just dotted cells
          initialHeaderHeight = STRICT_THRESHOLDS.BASE;
        } else {
          // Multiple sequences: show both tracks by default
          initialHeaderHeight = STRICT_THRESHOLDS.WITH_BOTH;
        }

        let webLogoTrackRef: LazyWebLogoTrack | null = null;
        let conservationTrackRef: LazyConservationTrack | null = null;
        const filterChangeSub = DG.debounce(df.onFilterChanged, 100).subscribe(() => {
          MSAViewportManager.clearAllCaches();

          if (webLogoTrackRef) {
            webLogoTrackRef.resetViewportTracking();
            webLogoTrackRef.forceUpdate();
          }
          if (conservationTrackRef) {
            conservationTrackRef.resetViewportTracking();
            conservationTrackRef.forceUpdate();
          }
          const newHasMultipleSequences = getHasMultipleSequences();
          if (newHasMultipleSequences !== hasMultipleSequences) {
            hasMultipleSequences = newHasMultipleSequences;
            let newHeaderHeight: number;
            if (!hasMultipleSequences)
              newHeaderHeight = STRICT_THRESHOLDS.BASE;
            else
              newHeaderHeight = STRICT_THRESHOLDS.WITH_BOTH;


            if (newHeaderHeight !== initialHeaderHeight) {
              initialHeaderHeight = newHeaderHeight;
              grid.props.colHeaderHeight = newHeaderHeight;
            }
          }
          grid.invalidate();
          setTimeout(() => {
            if (!grid.isDetached)
              grid.invalidate();
          }, 50);
        });

        grid.sub(filterChangeSub);

        const initializeHeaders = (monomerLib: any = null) => {
          const tracks: { id: string, track: MSAHeaderTrack, priority: number }[] = [];

          // Create lazy tracks only if we have multiple sequences
          if (hasMultipleSequences) {
            // OPTIMIZED: Pass seqHandler directly instead of column/splitter
            const conservationTrack = new LazyConservationTrack(
              sh,
              maxSeqLen,
              45, // DEFAULT_TRACK_HEIGHT
              'default',
              'Conservation'
            );
            conservationTrackRef = conservationTrack; // Store reference
            tracks.push({id: 'conservation', track: conservationTrack, priority: 1});

            // OPTIMIZED: Pass seqHandler directly
            const webLogoTrack = new LazyWebLogoTrack(
              sh,
              maxSeqLen,
              45, // DEFAULT_TRACK_HEIGHT
              'WebLogo'
            );
            webLogoTrackRef = webLogoTrack; // Store reference

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
              if (grid && !grid.isDetached) {
                const validHeight = Math.max(STRICT_THRESHOLDS.BASE, newHeight);
                grid.props.colHeaderHeight = validHeight;
              }
            },
          });

          scroller.setupTooltipHandling();

          // Add tracks to scroller
          tracks.forEach(({id, track}) => {
            scroller.addTrack(id, track);
          });

          scroller.setSelectionData(df, seqCol, sh);

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

            scroller.draw(
              cellBounds.x + startPadding,
              cellBounds.y,
              cellBounds.width - startPadding,
              cellBounds.height,
              getCurrent(),
              start,
              e,
              seqCol.name
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

  const _ = grok.events.onViewerAdded.subscribe((e) => {
    if (!e.args || !(e.args.viewer instanceof DG.Grid)) return;
    const grid = e.args.viewer as DG.Grid;
    handleGrid(grid);
  });

  const openTables = grok.shell.tableViews;
  for (const tv of openTables) {
    const grid = tv?.grid;
    if (grid) handleGrid(grid);
  }
}
