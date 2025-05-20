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


/**
 * Calculate WebLogo data for each position in a multiple sequence alignment
 * @param {string[]} sequences - Array of sequence strings
 * @param {SplitterFunc} splitter - Function to split sequences into monomers
 * @return {Map<number, Map<string, number>>} Map of position to monomer frequencies
 */
function calculateWebLogoData(sequences: string[], splitter: SplitterFunc): Map<number, Map<string, number>> {
  if (!sequences || sequences.length <= 1) return new Map();

  // Split sequences using the specialized splitter function
  const splitSequences: ISeqSplitted[] = sequences
    .map((seq) => splitter(seq))
    .filter(Boolean);

  if (splitSequences.length <= 1) return new Map();

  // Find the maximum sequence length
  const maxLen = Math.max(...splitSequences.map((seq) => seq.length));

  // Create a map to store residue frequencies at each position
  const webLogoData: Map<number, Map<string, number>> = new Map();

  // For each position in the alignment
  for (let pos = 0; pos < maxLen; pos++) {
    const residueCounts: Map<string, number> = new Map();
    let totalResidues = 0;

    // Count the residues at this position
    for (const seq of splitSequences) {
      if (pos < seq.length) {
        // Use getCanonical method to get standardized residue representation
        const residue = seq.getCanonical(pos);

        // Skip gaps - using the interface's method to check
        if (residue && !seq.isGap(pos)) {
          residueCounts.set(residue, (residueCounts.get(residue) || 0) + 1);
          totalResidues++;
        }
      }
    }

    if (totalResidues === 0) {
      webLogoData.set(pos, new Map()); // Empty map for positions with no data
      continue;
    }

    // Convert counts to frequencies
    const residueFreqs: Map<string, number> = new Map();
    for (const [residue, count] of residueCounts.entries())
      residueFreqs.set(residue, count / totalResidues);


    webLogoData.set(pos, residueFreqs);
  }

  return webLogoData;
}
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

        // Get maximum sequence length
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

        // Skip if sequence is too short
        if (maxSeqLen < 50)
          continue;

        // Calculate conservation data if we have multiple sequences
        let conservationData: number[] = [];
        let webLogoData: Map<number, Map<string, number>> = new Map();

        if (cats.length > 1) {
          // Calculate both conservation and WebLogo data
          conservationData = calculateConservation(cats.filter(Boolean), sh.splitter);
          webLogoData = calculateWebLogoData(cats.filter(Boolean), sh.splitter);
        }

        // Determine initial header height based on available data
        const hasConservation = conservationData.length > 0;
        const hasWebLogo = webLogoData.size > 0;

        // Define track heights, layout parameters, and thresholds
        const dottedCellsHeight = 38; // Fixed height for dotted cells + slider
        const trackGap = 4; // Gap between tracks

        // Default track heights with title space included
        const DEFAULT_WEBLOGO_HEIGHT = 45; // Increased from 40 to accommodate title
        const DEFAULT_CONSERVATION_HEIGHT = 45; // Increased from 40 to accommodate title

        // Calculate total required height with all tracks
        const initialHeaderHeight = dottedCellsHeight +
          (hasWebLogo ? DEFAULT_WEBLOGO_HEIGHT + trackGap : 0) +
          (hasConservation ? DEFAULT_CONSERVATION_HEIGHT + trackGap : 0);

        // Create tracks in the correct vertical order (Conservation on top, WebLogo in middle)
        const tracks: { id: string, track: MSAHeaderTrack, priority: number }[] = [];

        // Function to initialize headers with monomerLib
        const initializeHeaders = (monomerLib: any = null) => {
          // 1. Conservation track first (top position)
          if (hasConservation) {
            const conservationTrack = new ConservationTrack(
              conservationData,
              DEFAULT_CONSERVATION_HEIGHT,
              'default', // color scheme
              'Conservation' // Track title
            );
            tracks.push({id: 'conservation', track: conservationTrack, priority: 1}); // Lower visibility priority (1)
          }

          // 2. WebLogo track second (middle position)
          if (hasWebLogo) {
            const webLogoTrack = new WebLogoTrack(
              webLogoData,
              DEFAULT_WEBLOGO_HEIGHT,
              'default', // This parameter is ignored in our simplified implementation
              'WebLogo' // Track title
            );

            // Set monomerLib for color lookup if available
            if (monomerLib) {
              webLogoTrack.setMonomerLib(monomerLib);
              webLogoTrack.setBiotype(sh.defaultBiotype || 'PEPTIDE');
              console.log(`WebLogo track: using monomerLib with biotype ${sh.defaultBiotype || 'PEPTIDE'}`);
            } else
              console.log('WebLogo track: no monomerLib available, using fallback colors');


            webLogoTrack.setupDefaultTooltip(); // Enable tooltips
            tracks.push({id: 'weblogo', track: webLogoTrack, priority: 2}); // Higher visibility priority (2)
          }

          // Create MSAScrollingHeader with initial height
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
          });

          scroller.setupTooltipHandling();
          // Add tracks to scroller in the defined vertical order (top to bottom)
          tracks.forEach(({id, track}) => {
            scroller.addTrack(id, track);
          });

          // Set initial header height in grid
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

            // Always set header height based on available space
            scroller.headerHeight = cellBounds.height;

            // Set position width based on font size
            const font = getFontSize();
            scroller.positionWidth = font + (isMSA ? 8 : 0);

            // Get current state
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

        // Get monomerLib helper to access the monomer library
        getMonomerLibHelper()
          .then((libHelper) => {
            const monomerLib = libHelper.getMonomerLib();
            initializeHeaders(monomerLib);
          })
          .catch((error) => {
            console.error('Error loading monomerLib:', error);
            initializeHeaders(); // Initialize with fallback colors
          });
      }
    }, 1000);
  };

  // Handle new grid additions
  const _ = grok.events.onViewerAdded.subscribe((e) => {
    if (!e.args || !(e.args.viewer instanceof DG.Grid))
      return;
    const grid = e.args.viewer as DG.Grid;
    handleGrid(grid);
  });

  // Handle existing grids
  const openTables = grok.shell.tableViews;
  for (const tv of openTables) {
    const grid = tv?.grid;
    if (grid)
      handleGrid(grid);
  }
}
