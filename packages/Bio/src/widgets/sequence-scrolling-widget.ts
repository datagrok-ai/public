/* eslint-disable rxjs/no-ignored-subscription */
/* eslint-disable max-lines */
/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {ConservationTrack, MSAHeaderTrack, MSAScrollingHeader, WebLogoTrack} from '@datagrok-libraries/bio/src/utils/sequence-position-scroller';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {ALPHABET, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_package} from '../package';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import * as RxJs from 'rxjs';
import {filter} from 'rxjs/operators';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';

// ============================================================================
// OPTIMIZED VIEWPORT-AWARE CACHING WITH FORCE UPDATE SUPPORT
// ============================================================================
class MSAViewportManager {
  private static conservationCache: DG.LruCache<string, number[]> = new DG.LruCache<string, number[]>(100);
  private static webLogoCache: DG.LruCache<string, Map<number, Map<string, number>>> = new DG.LruCache<string, Map<number, Map<string, number>>>(100);

  // Track when data was last invalidated for force updates
  private static lastInvalidationTime: number = 0;

  // Cache chunks of 200 positions at a time
  private static readonly CHUNK_SIZE = 50;

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

  private static getChunkCacheKey(column: DG.Column, chunkStart: number, chunkEnd: number): string {
    const dataFrame = column.dataFrame;
    // Use framework's built-in versioning instead of complex hash
    return `${dataFrame.id}_${column.name}_f${dataFrame.filter.version}_${chunkStart}_${chunkEnd}`;
  }

  private static getConservationChunk(seqHandler: ISeqHandler, start: number, end: number, cacheKey: string): number[] {
    return MSAViewportManager.conservationCache.getOrCreate(cacheKey, () => {
      // conservation track is only shown together with weblogo track, so we can reuse weblogo for this caluclation
      const webLogoChunk = MSAViewportManager.getWebLogoChunk(seqHandler, start, end, cacheKey);
      const chunkSize = end - start;
      const scores = new Array(chunkSize).fill(0);
      for (let i = start; i < end; i++) {
        const entries = webLogoChunk.get(i);
        if (!entries || entries.size === 0)
          continue;
        let totalResidues = 0;
        let mostCommonCount = 0;
        const gapS = seqHandler.defaultGapOriginal;
        // Single pass to find both total and max
        for (const [mon, count] of entries.entries()) {
          if (!mon || mon === gapS) continue; // Skip gaps
          totalResidues += count;
          if (count > mostCommonCount)
            mostCommonCount = count;
        }
        // Calculate conservation score
        scores[i - start] = totalResidues > 0 ? mostCommonCount / totalResidues : 0;
      }
      return scores;
    });
  }

  private static getWebLogoChunk(seqHandler: ISeqHandler, start: number, end: number, cacheKey: string): Map<number, Map<string, number>> {
    return MSAViewportManager.webLogoCache.getOrCreate(cacheKey, () => {
      const dataFrame = seqHandler.column.dataFrame;
      const webLogoData: Map<number, Map<string, number>> = new Map();
      // OPTIMIZED: Use filter iterator
      const visibleCount = dataFrame.filter.trueCount;
      if (visibleCount <= 1)
        return webLogoData;

      // pre-alocate residue counts map
      for (let i = start; i < end; i++)
        webLogoData.set(i, new Map());

      const filter = dataFrame.filter;

      const oneOverCount = 1 / visibleCount; // Pre-calculate for performance

      for (let i = -1; (i = filter.findNext(i, true)) !== -1;) {
        // OPTIMIZED: Use seqHandler.getSplitted which has internal caching
        const seqChunk = seqHandler.getSplitted(i).getOriginalRegion(start, end);
        if (seqChunk.length === 0) continue;
        for (let pos = 0; pos < seqChunk.length; pos++) {
          const residue = seqChunk[pos];
          // in weblogo, do not skip gaps, they are important
          const residueMap = webLogoData.get(start + pos)!;
          // add oneOverCount here to avoid division later
          residueMap.set(residue, (residueMap.get(residue) || 0) + oneOverCount);
        }
      }
      return webLogoData;
    });
  }

  static getConservationForViewport(seqHandler: ISeqHandler, viewportStart: number, viewportEnd: number, maxLength: number): number[] {
    // Create a full-sized array filled with zeros
    const result: number[] = new Array(maxLength).fill(0);

    // Calculate which chunks we need
    const startChunk = Math.floor(viewportStart / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;
    const endChunk = Math.ceil(viewportEnd / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;

    // Get all needed chunks and place them in the correct positions
    for (let chunkStart = startChunk; chunkStart < endChunk; chunkStart += MSAViewportManager.CHUNK_SIZE) {
      const chunkEnd = Math.min(chunkStart + MSAViewportManager.CHUNK_SIZE, maxLength);
      const cacheKey = MSAViewportManager.getChunkCacheKey(seqHandler.column, chunkStart, chunkEnd);
      const chunkData = MSAViewportManager.getConservationChunk(seqHandler, chunkStart, chunkEnd, cacheKey);

      // Copy chunk data to the correct positions in the result array
      for (let i = 0; i < chunkData.length && (chunkStart + i) < maxLength; i++)
        result[chunkStart + i] = chunkData[i];
    }

    return result;
  }

  static getWebLogoForViewport(seqHandler: ISeqHandler, viewportStart: number, viewportEnd: number, maxLength: number): Map<number, Map<string, number>> {
    const result: Map<number, Map<string, number>> = new Map();

    // Calculate which chunks we need
    const startChunk = Math.floor(viewportStart / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;
    const endChunk = Math.ceil(viewportEnd / MSAViewportManager.CHUNK_SIZE) * MSAViewportManager.CHUNK_SIZE;

    // Get all needed chunks and place them in the correct positions
    for (let chunkStart = startChunk; chunkStart < endChunk; chunkStart += MSAViewportManager.CHUNK_SIZE) {
      const chunkEnd = Math.min(chunkStart + MSAViewportManager.CHUNK_SIZE, maxLength);
      const cacheKey = MSAViewportManager.getChunkCacheKey(seqHandler.column, chunkStart, chunkEnd);
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
  private seqHandler: ISeqHandler;
  private maxLength: number;
  private lastViewportStart: number = -1;
  private lastViewportEnd: number = -1;
  private lastInvalidationTime: number = 0;
  private forceNextUpdate: boolean = false;

  constructor(
    seqHandler: ISeqHandler,
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
    const bufferedStart = Math.max(0, viewportStart - 20);
    const bufferedEnd = Math.min(this.maxLength, viewportEnd + 20);

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
  private seqHandler: ISeqHandler;
  private maxLength: number;
  private lastViewportStart: number = -1;
  private lastViewportEnd: number = -1;
  private lastInvalidationTime: number = 0;
  private forceNextUpdate: boolean = false;

  constructor(
    seqHandler: ISeqHandler,
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
    const bufferedStart = Math.max(0, viewportStart - 20);
    const bufferedEnd = Math.min(this.maxLength, viewportEnd + 20);

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
        let sh: ISeqHandler | null = null;

        try {
          sh = _package.seqHelper.getSeqHandler(seqCol);
        } catch (_e) {
          continue;
        }
        if (!sh) continue;
        if (sh.isHelm() || sh.alphabet === ALPHABET.UN) continue;

        const gCol = grid.col(seqCol.name);
        if (!gCol) continue;

        let positionStatsViewerAddedOnce = !!grid.tableView &&
          Array.from(grid.tableView.viewers).some((v) => v.type === 'Sequence Position Statistics');

        const isMSA = sh.isMsa();

        if (isMSA) {
          const ifNan = (a: number, els: number) => (Number.isNaN(a) ? els : a);
          const getStart = () => ifNan(Math.max(Number.parseInt(seqCol.getTag(bioTAGS.positionShift) ?? '0'), 0), 0) + 1;
          const getCurrent = () => ifNan(Number.parseInt(seqCol.getTag(bioTAGS.selectedPosition) ?? '-2'), -2);
          const getFontSize = () => MonomerPlacer.getFontSettings(seqCol).fontWidth;

          // Get maximum sequence length. since this scroller is only applicable to Single character monomeric sequences,
          // we do not need to check every single sequence and split it, instead, max length will coorelate with length of the longest string
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

          // Do not Skip if sequences are too short, rather, just don't render the tracks by default

          const STRICT_THRESHOLDS = {
            WITH_TITLE: 58, // BASE + TITLE_HEIGHT(16) + TRACK_GAP(4)
            WITH_WEBLOGO: 107, // WITH_TITLE + DEFAULT_TRACK_HEIGHT(45) + TRACK_GAP(4)
            WITH_BOTH: 156 // WITH_WEBLOGO + DEFAULT_TRACK_HEIGHT(45) + TRACK_GAP(4)
          };

          let initialHeaderHeight: number;
          if (seqCol.length > 100_000 || maxSeqLen < 50) {
            // Single sequence: just dotted cells
            initialHeaderHeight = STRICT_THRESHOLDS.WITH_TITLE;
          } else {
            if (seqCol.length > 50_000)
              initialHeaderHeight = STRICT_THRESHOLDS.WITH_WEBLOGO;
            else
              initialHeaderHeight = STRICT_THRESHOLDS.WITH_BOTH;
          }

          let webLogoTrackRef: LazyWebLogoTrack | null = null;
          let conservationTrackRef: LazyConservationTrack | null = null;
          const filterChangeSub = DG.debounce(
            RxJs.merge(df.onFilterChanged, df.onDataChanged.pipe(filter((a) => a?.args?.column === seqCol))), 100
          ).subscribe(() => {
            MSAViewportManager.clearAllCaches();

            if (webLogoTrackRef) {
              webLogoTrackRef.resetViewportTracking();
              webLogoTrackRef.forceUpdate();
            }
            if (conservationTrackRef) {
              conservationTrackRef.resetViewportTracking();
              conservationTrackRef.forceUpdate();
            }
            setTimeout(() => {
              if (!grid.isDetached)
                grid.invalidate();
            }, 50);
          });

          grid.sub(filterChangeSub);

          const initializeHeaders = (monomerLib: IMonomerLib) => {
            const tracks: { id: string, track: MSAHeaderTrack, priority: number }[] = [];

            // Create lazy tracks only for MSA sequences

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
              webLogoTrack.setBiotype(sh.defaultBiotype || 'HELM_AA');
            }

            webLogoTrack.setupDefaultTooltip();
            tracks.push({id: 'weblogo', track: webLogoTrack, priority: 2});

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
              onHeaderHeightChange: (_newHeight) => {
                // Update grid header height
                if (grid.isDetached || _newHeight < STRICT_THRESHOLDS.WITH_TITLE) return;
                setTimeout(() => grid.props.colHeaderHeight = _newHeight);
              },
            });

            scroller.setupTooltipHandling();

            // Add tracks to scroller
            tracks.forEach(({id, track}) => {
              scroller.addTrack(id, track);
            });

            scroller.setSelectionData(df, seqCol, sh);

            if (maxSeqLen > 50) {
              grid.props.colHeaderHeight = initialHeaderHeight;

              // Set column width
              setTimeout(() => {
                if (grid.isDetached) return;
                gCol.width = 400;
              }, 300);
            }

            // Handle cell rendering for MSA
            grid.sub(grid.onCellRender.subscribe((e) => {
              const cell = e.cell;
              if (!cell || !cell.isColHeader || cell?.gridColumn?.name !== gCol?.name)
                return;

              const cellBounds = e.bounds;
              if (!cellBounds) return;

              // Set dynamic properties
              scroller.headerHeight = cellBounds.height;
              const font = getFontSize();
              scroller.positionWidth = font + 8; // MSA always has padding

              const start = getStart();
              const startPadding = 0; // MSA doesn't need extra padding

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

          // Initialize with monomer library for MSA sequences
          getMonomerLibHelper()
            .then((libHelper) => {
              const monomerLib = libHelper.getMonomerLib();
              initializeHeaders(monomerLib);
            })
            .catch((error) => {
              grok.shell.warning(`Failed to initialize monomer library`);
              console.error('Failed to initialize monomer library:', error);
            });
        } else {
          // For non-MSA sequences, just use standard sequence rendering.
        }
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
