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
// VIEWPORT-AWARE CACHING WITH FORCE UPDATE SUPPORT
// ============================================================================

/**
 * Lazy viewport-aware MSA data manager with force update capability
 */
class MSAViewportManager {
  private static conservationCache: DG.LruCache<string, number[]> = new DG.LruCache<string, number[]>(100);
  private static webLogoCache: DG.LruCache<string, Map<number, Map<string, number>>> = new DG.LruCache<string, Map<number, Map<string, number>>>(100);

  // Track when data was last invalidated for force updates
  private static lastInvalidationTime: number = 0;

  // Cache chunks of 200 positions at a time
  private static readonly CHUNK_SIZE = 200;

  private static simpleHash(str: string): string {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
      const char = str.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash;
    }
    return Math.abs(hash).toString(36);
  }

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

  /**
   * Modified cache key generation that includes filter state
   */
  private static getChunkCacheKey(column: DG.Column, chunkStart: number, chunkEnd: number, type: string): string {
    // Get only the visible/filtered sequences
    const dataFrame = column.dataFrame;
    const visibleSequences: string[] = [];

    if (dataFrame.filter.trueCount === dataFrame.rowCount) {
      // No filters applied, use all sequences
      visibleSequences.push(...column.categories.filter(Boolean));
    } else {
      // Filters applied, get only visible sequences
      for (let i = 0; i < dataFrame.rowCount; i++) {
        if (dataFrame.filter.get(i)) {
          const seq = column.get(i);
          if (seq) visibleSequences.push(seq);
        }
      }
    }

    const dataHash = MSAViewportManager.simpleHash(visibleSequences.join('|'));
    const filterHash = dataFrame.filter.trueCount; // Simple filter state indicator

    return `${column.dataFrame.id}_${column.name}_v${column.version}_${dataHash}_f${filterHash}_${type}_chunk_${chunkStart}_${chunkEnd}`;
  }

  /**
   * Conservation chunk method that works with filtered data
   */
  private static getConservationChunk(column: DG.Column, splitter: SplitterFunc, start: number, end: number): number[] {
    const cacheKey = MSAViewportManager.getChunkCacheKey(column, start, end, 'conservation');

    return MSAViewportManager.conservationCache.getOrCreate(cacheKey, () => {
      // Get only visible/filtered sequences
      const dataFrame = column.dataFrame;
      const sequences: string[] = [];

      if (dataFrame.filter.trueCount === dataFrame.rowCount)
        sequences.push(...column.categories.filter(Boolean));
      else {
        for (let i = 0; i < dataFrame.rowCount; i++) {
          if (dataFrame.filter.get(i)) {
            const seq = column.get(i);
            if (seq) sequences.push(seq);
          }
        }
      }

      if (sequences.length <= 1)
        return new Array(end - start).fill(0);

      const splitSequences: ISeqSplitted[] = sequences
        .map((seq) => splitter(seq))
        .filter(Boolean);

      const scores: number[] = [];

      // Only calculate for the requested range
      for (let pos = start; pos < end; pos++) {
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
   * WebLogo chunk method that works with filtered data
   */
  private static getWebLogoChunk(column: DG.Column, splitter: SplitterFunc, start: number, end: number): Map<number, Map<string, number>> {
    const cacheKey = MSAViewportManager.getChunkCacheKey(column, start, end, 'weblogo');

    return MSAViewportManager.webLogoCache.getOrCreate(cacheKey, () => {
      // Get only visible/filtered sequences
      const dataFrame = column.dataFrame;
      const sequences: string[] = [];

      if (dataFrame.filter.trueCount === dataFrame.rowCount)
        sequences.push(...column.categories.filter(Boolean));
      else {
        for (let i = 0; i < dataFrame.rowCount; i++) {
          if (dataFrame.filter.get(i)) {
            const seq = column.get(i);
            if (seq) sequences.push(seq);
          }
        }
      }

      const webLogoData: Map<number, Map<string, number>> = new Map();

      if (sequences.length <= 1)
        return webLogoData;

      const splitSequences: ISeqSplitted[] = sequences
        .map((seq) => splitter(seq))
        .filter(Boolean);

      // Only calculate for the requested range
      for (let pos = start; pos < end; pos++) {
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

  /**
   * Get conservation scores for viewport (returns properly aligned sparse array)
   */
  static getConservationForViewport(column: DG.Column, splitter: SplitterFunc, viewportStart: number, viewportEnd: number, maxLength: number): number[] {
    // Create a full-sized array filled with zeros
    const result: number[] = new Array(maxLength).fill(0);

    // Calculate which chunks we need
    const startChunk = Math.floor(viewportStart / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;
    const endChunk = Math.ceil(viewportEnd / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;

    // Get all needed chunks and place them in the correct positions
    for (let chunkStart = startChunk; chunkStart < endChunk; chunkStart += MSAViewportManager.CHUNK_SIZE) {
      const chunkEnd = Math.min(chunkStart + MSAViewportManager.CHUNK_SIZE, maxLength);
      const chunkData = MSAViewportManager.getConservationChunk(column, splitter, chunkStart, chunkEnd);

      // Copy chunk data to the correct positions in the result array
      for (let i = 0; i < chunkData.length && (chunkStart + i) < maxLength; i++)
        result[chunkStart + i] = chunkData[i];
    }

    return result;
  }

  /**
   * Get WebLogo data for viewport (returns properly aligned sparse map)
   */
  static getWebLogoForViewport(column: DG.Column, splitter: SplitterFunc, viewportStart: number, viewportEnd: number, maxLength: number): Map<number, Map<string, number>> {
    const result: Map<number, Map<string, number>> = new Map();

    // Calculate which chunks we need
    const startChunk = Math.floor(viewportStart / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;
    const endChunk = Math.ceil(viewportEnd / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;

    // Get all needed chunks and place them in the correct positions
    for (let chunkStart = startChunk; chunkStart < endChunk; chunkStart += MSAViewportManager.CHUNK_SIZE) {
      const chunkEnd = Math.min(chunkStart + MSAViewportManager.CHUNK_SIZE, maxLength);
      const chunkData = MSAViewportManager.getWebLogoChunk(column, splitter, chunkStart, chunkEnd);

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
// SMART TRACKS WITH FORCE UPDATE CAPABILITY
// ============================================================================

/**
 * WebLogoTrack that loads data on-demand with force update capability
 */
class LazyWebLogoTrack extends WebLogoTrack {
  private column: DG.Column;
  private splitter: SplitterFunc;
  private maxLength: number;
  private lastViewportStart: number = -1;
  private lastViewportEnd: number = -1;
  private lastInvalidationTime: number = 0;
  private forceNextUpdate: boolean = false;

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
    this.maxLength = maxLength;

    // Check if we have data to show (but don't calculate it yet!)
    const sequences = column.categories.filter(Boolean);
    this.visible = sequences.length > 1;
  }

  /**
   * Force the next update regardless of viewport
   */
  public forceUpdate(): void {
    this.forceNextUpdate = true;
  }

  /**
   * Reset viewport tracking to force update on next draw
   */
  public resetViewportTracking(): void {
    this.lastViewportStart = -1;
    this.lastViewportEnd = -1;
  }

  /**
   * Update data only for the positions that are actually being rendered
   */
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

    // Get properly aligned data for the tracks
    const webLogoData = MSAViewportManager.getWebLogoForViewport(
      this.column,
      this.splitter,
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

/**
 * ConservationTrack that loads data on-demand with force update capability
 */
class LazyConservationTrack extends ConservationTrack {
  private column: DG.Column;
  private splitter: SplitterFunc;
  private maxLength: number;
  private lastViewportStart: number = -1;
  private lastViewportEnd: number = -1;
  private lastInvalidationTime: number = 0;
  private forceNextUpdate: boolean = false;

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
    this.maxLength = maxLength;

    // Check if we have data to show (but don't calculate it yet!)
    const sequences = column.categories.filter(Boolean);
    this.visible = sequences.length > 1;
  }

  /**
   * Force the next update regardless of viewport
   */
  public forceUpdate(): void {
    this.forceNextUpdate = true;
  }

  /**
   * Reset viewport tracking to force update on next draw
   */
  public resetViewportTracking(): void {
    this.lastViewportStart = -1;
    this.lastViewportEnd = -1;
  }

  /**
   * Update data only for the positions that are actually being rendered
   */
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

    // Get properly aligned data for the tracks
    const conservationScores = MSAViewportManager.getConservationForViewport(
      this.column,
      this.splitter,
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

    // Update data for just the visible range
    this.updateForViewport(windowStart - 1, viewportEnd - 1); // Convert to 0-based

    // Call parent draw method
    super.draw(x, y, width, height, windowStart, positionWidth, totalPositions, currentPosition);
  }
}

// ============================================================================
// MAIN HANDLER
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
        const getHasMultipleSequences = () => getVisibleSequenceCount(seqCol) > 1;
        let hasMultipleSequences = getHasMultipleSequences();

        const STRICT_THRESHOLDS = {
          BASE: 38, // DOTTED_CELL_HEIGHT(30) + SLIDER_HEIGHT(8)
          WITH_TITLE: 58, // BASE + TITLE_HEIGHT(16) + TRACK_GAP(4)
          WITH_WEBLOGO: 107, // WITH_TITLE + DEFAULT_TRACK_HEIGHT(45) + TRACK_GAP(4)
          WITH_BOTH: 156 // WITH_WEBLOGO + DEFAULT_TRACK_HEIGHT(45) + TRACK_GAP(4)
        };

        // Calculate initial header height based on available data
        let initialHeaderHeight: number;
        if (!hasMultipleSequences) {
          // Single sequence: just dotted cells
          initialHeaderHeight = STRICT_THRESHOLDS.BASE;
        } else {
          // Multiple sequences: show both tracks by default
          initialHeaderHeight = STRICT_THRESHOLDS.WITH_BOTH;
        }

        // Store track references for force updates
        let webLogoTrackRef: LazyWebLogoTrack | null = null;
        let conservationTrackRef: LazyConservationTrack | null = null;

        // ADD FILTER CHANGE SUBSCRIPTION WITH FORCE UPDATE CAPABILITY
        const filterChangeSub = DG.debounce(df.onFilterChanged, 100).subscribe(() => {
          // 1. Clear all caches (this also updates invalidation timestamp)
          MSAViewportManager.clearAllCaches();

          // 2. Force track updates by resetting viewport tracking
          if (webLogoTrackRef) {
            webLogoTrackRef.resetViewportTracking();
            webLogoTrackRef.forceUpdate();
          }
          if (conservationTrackRef) {
            conservationTrackRef.resetViewportTracking();
            conservationTrackRef.forceUpdate();
          }

          // 3. Update hasMultipleSequences based on filtered data
          const newHasMultipleSequences = getHasMultipleSequences();
          if (newHasMultipleSequences !== hasMultipleSequences) {
            hasMultipleSequences = newHasMultipleSequences;

            // Potentially adjust header height if sequence count changed significantly
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

          // 4. Force grid redraw
          grid.invalidate();

          // 5. Additional redraw attempts to ensure update
          setTimeout(() => {
            if (!grid.isDetached)
              grid.invalidate();
          }, 50);
        });

        // Store the subscription for cleanup if needed
        if (!grid.temp) grid.temp = {};
        if (!grid.temp.filterSubs) grid.temp.filterSubs = [];
        grid.temp.filterSubs.push(filterChangeSub);

        // Initialize tracks with cached data
        const initializeHeaders = (monomerLib: any = null) => {
          const tracks: { id: string, track: MSAHeaderTrack, priority: number }[] = [];

          // Create lazy tracks only if we have multiple sequences
          if (hasMultipleSequences) {
            // Conservation track
            const conservationTrack = new LazyConservationTrack(
              seqCol,
              sh.splitter,
              maxSeqLen,
              45, // DEFAULT_TRACK_HEIGHT
              'default',
              'Conservation'
            );
            conservationTrackRef = conservationTrack; // Store reference
            tracks.push({id: 'conservation', track: conservationTrack, priority: 1});

            // WebLogo track
            const webLogoTrack = new LazyWebLogoTrack(
              seqCol,
              sh.splitter,
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

// Helper function to get visible sequence count
function getVisibleSequenceCount(seqCol: DG.Column): number {
  const dataFrame = seqCol.dataFrame;
  if (dataFrame.filter.trueCount === dataFrame.rowCount)
    return seqCol.categories.filter(Boolean).length;

  const visibleSequences = new Set<string>();
  for (let i = 0; i < dataFrame.rowCount; i++) {
    if (dataFrame.filter.get(i)) {
      const seq = seqCol.get(i);
      if (seq) visibleSequences.add(seq);
    }
  }
  return visibleSequences.size;
}

// Optional: Add cleanup method for when grids are destroyed
export function cleanupMSAHeaderSubscriptions(grid: DG.Grid) {
  if (grid.temp && grid.temp.filterSubs) {
    grid.temp.filterSubs.forEach((sub: any) => sub.unsubscribe());
    grid.temp.filterSubs = [];
  }
}
