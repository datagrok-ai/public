/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {getSeqHelper, ISeqHelper} from './seq-helper';
import {ALPHABET, MonomerToShortFunc, NOTATION, SplitterFunc, TAGS as bioTAGS} from './macromolecule';
import {ISeqGraphInfoHelm, ISeqSplitted} from './macromolecule/types';
import {CellRendererBackBase} from './cell-renderer-back-base';
import {ILogger} from './logger';
import {getMonomerLibHelper} from '../types/monomer-library';
import {errInfo} from './err-info';
import {DrawStyle, getPrintColumnTemps, printLeftOrCentered, TAGS as mmcrTAGS, PrintOptions} from './cell-renderer';
import {MmcrTemps, rendererSettingsChangedState, tempTAGS} from './cell-renderer-consts';
import {AnnotationRenderer} from './cell-renderer-annotations';
import {IMonomerLibBase} from '../types/monomer-library';
import {HelmTypes} from '../helm/consts';
import {ISeqMonomer} from '../helm/types';
import {execMonomerHoverLinks} from '../monomer-works/monomer-hover';
import {
  DEFAULT_MACROMOLECULE_HIGHLIGHT_FILL,
  MACROMOLECULE_HIGHLIGHT_EVENT_ID, MACROMOLECULE_HIGHLIGHT_TEMP,
  MacromoleculeHighlightEntry, MacromoleculeHighlightEventArgs, macromoleculeHighlightColorToCss,
} from './macromolecule-highlight';
import {fromEvent, Unsubscribable} from 'rxjs';
import * as operators from 'rxjs/operators';
import {
  CONNECTION_GREY, IConnectionsGeometry, IConnectionsLayout, IDrawnConnection, assignLanes, bondNumberFontSize, connectionsGeometry, darkerColor,
  drawBondNumbers, drawConnections, getDrawnConnections, layoutConnections, measureConnectionsFont, sideHeight,
  stackedLineGap,
} from './cell-renderer-connections';

type MonomerPlacerProps = {
  separatorWidth: number,
  monomerToShort: MonomerToShortFunc,
  font: string,
  fontCharWidth: number,
  /** The size of {@link font}; when not given, the column's font size temp or 12 */
  fontSize?: number,
  /** The max monomer length when the column sets none (the package default); when not given, the constructor's */
  maxMonomerLength?: number,
  /** Whether rows with connections draw them as lines when the column has no `.mm.cellRenderer.showConnections` tag
   * (the package default); off when not given */
  showConnections?: boolean,
};

export const undefinedColor = 'rgb(100,100,100)';

export const shiftedLeftPaddingText = '...';

/** Be ware, this can return -1 meaning that hovering/clicking happened on the three dots
 * and not on the monomer itself */
export function hitBounds(bounds: number[], x: number, positionShiftPixels?: number): number | null {
  if ((positionShiftPixels ?? 0) > 0 && x < (bounds[0] ?? 0) + positionShiftPixels!)
    return -1;
  x -= positionShiftPixels ?? 0;

  let iterationCount: number = 100;
  let leftI: number = 0;
  let rightI = bounds.length - 1;
  let midI;
  while (leftI <= rightI) {
    midI = Math.floor((rightI + leftI) / 2);
    if (bounds[midI] <= x && x < bounds[midI + 1])
      return midI;
    else if (x < bounds[midI])
      rightI = midI - 1;
    else /* if (bounds[midI + 1] <= x) */
      leftI = midI + 1;

    if (--iterationCount <= 0) {
      throw new Error(`Get position for pointer x = ${x} searching has not converged ` +
        `on ${JSON.stringify(bounds)}. `);
    }
  }
  return null;
}

interface IMonomerLayoutData {
  lineIdx: number;
  monomerIdx: number;
  bounds: DG.Rect;
  sequencePosition: number;
}

/** What a row drawn one line per chain needs besides its cell */
interface IStackedRow {
  sh: { isSeparator(): boolean };
  seqSS: ISeqSplitted;
  layout: IConnectionsLayout;
  maxLengthOfMonomer: number;
  separator: string;
  monomerColor: (cm: string) => string;
  referenceSequence: string[];
  selectedPosition: number;
  highlightEntry: MacromoleculeHighlightEntry | null;
  annotRenderer: AnnotationRenderer | null;
  /** The table row, null for a render outside a grid */
  rowIdx: number | null;
  /** Hit tests count positions from it (a stacked row is laid out unshifted) */
  positionShift: number;
}

/** A stacked row's horizontal layout: per flat position its x from its chain's start, its chain and its label width;
 * per chain its cumulative x (printLeftOrCentered's maxWord) and its shift */
interface IStackedGeometry {
  sepW: number;
  localX: number[];
  chainOf: number[];
  labelW: number[];
  maxWords: number[][];
  offsets: number[];
}

const SOLID: number[] = [];
const DASHED = [2, 2];

interface IMultilineLayoutElement {
  posIdx?: number;
  x: number;
  width: number;
  om: string; // original monomer or separator text
  isSeparator: boolean;
}

/** {@link MonomerPlacer.hoverPartners} of a cell without a hovered connected monomer */
const NO_HOVER: { pos: number | null, partners: Set<number> } = Object.freeze({pos: null, partners: new Set<number>()});

export class MonomerPlacer extends CellRendererBackBase<string> {
  private colWidth: number = 0;
  private _monomerLengthList: number[][] | null = null;

  // width of separator symbol
  private get separatorWidth() { return this.props?.fontCharWidth ? this.props?.fontCharWidth : 5; }
  public props: MonomerPlacerProps;
  private _processedRows: DG.BitSet; // rows for which monomer lengths were processed
  private _processedMaxVisibleSeqLength: number = 0;

  public _monomerLengthMap: { [key: string]: TextMetrics } = {}; // caches the lengths to save time on g.measureText
  public _monomerStructureMap: { [key: string]: HTMLElement } = {}; // caches the atomic structures of monomers

  private _ellipsisBounds: DG.Rect | undefined = undefined;
  private _totalLinesNeeded: number = 0;
  private _lineHeight: number = 20;

  private _cellBounds: Map<number, IMonomerLayoutData[]> = new Map();

  /** Rows with connections draw them as lines: the column's tag, else the package default; resolved on reset */
  private _showConnections: boolean = false;
  /** The rows' connection layouts per cell height (a BILN renderer serves every grid of its column), null for a row
   * printed as before; emptied on reset */
  private _connLayouts = new Map<number, Map<number, IConnectionsLayout | null>>();
  /** Stacked rows' hit bounds by cell height, then row: a BILN renderer serves grids of other row heights, each hit-tested
   * against its own drawing (multiline rows keep theirs in _cellBounds) */
  private _stackedBounds = new Map<number, Map<number, IMonomerLayoutData[]>>();
  /** Stacked rows' geometries, by their layout (dropped with the layouts) */
  private _stackedGeometries = new WeakMap<IConnectionsLayout, IStackedGeometry>();
  /** The cell height the monomer widths are laid out for (the last drawn or hovered cell's): a row that draws its
   * connections prints canonical symbols */
  private _widthsH: number = -1;
  /** Per row: whether its cached monomer widths are of canonical symbols */
  private _monomerLengthCanonical: boolean[] = [];
  /** The cell height the MSA widths (per position, over the rows) were laid out for */
  private _msaWidthsH: number = -1;
  /** Monomer colours one shade darker, for connection lines */
  private _darkerColors = new Map<string, string>();
  /** The column's font size temp this renderer has applied: the first renderer to render clears the changed flag */
  private _appliedFontSizeTemp: any = undefined;
  /** When the props were last compared with the provider's (the package defaults change without an event) */
  private _propsCheckedAt: number = 0;
  /** Whether this render pass has compared the settings that change without an event (see checkSettingsOncePerPass) */
  private _passChecked: boolean = false;
  /** The connected monomer under the mouse (a table row and a flat position, in a grid): its partners are marked */
  private _hoverConn: { row: number, pos: number, grid: DG.Grid | null } | null = null;
  /** While a connected monomer is hovered: its grid, whose overlay is watched for the mouse leaving the cells */
  private _hoverLeave: { grid: DG.Grid, subs: Unsubscribable[] } | null = null;

  private seqHelper!: ISeqHelper;

  private sysMonomerLib: IMonomerLibBase | null = null;

  /** Lazily created annotation renderer — null cost when no annotations exist. */
  private _annotationRenderer: AnnotationRenderer | null = null;

  private getAnnotationRenderer(): AnnotationRenderer | null {
    if (!this._annotationRenderer)
      this._annotationRenderer = new AnnotationRenderer(this.tableCol);
    return this._annotationRenderer.hasAnnotations() ? this._annotationRenderer : null;
  }

  /** View is required to subscribe and handle for data frame changes */
  constructor(
    gridCol: DG.GridColumn | null,
    tableCol: DG.Column<string>,
    logger: ILogger,
    public monomerLengthLimit: number,
    private readonly propsProvider: () => MonomerPlacerProps,
  ) {
    super(gridCol, tableCol, logger);
    this.props = this.propsProvider();
    this._processedRows = DG.BitSet.create(this.tableCol!.length);
    this.subs.push({unsubscribe: () => this.watchHoverLeave(null)});

    if (this.gridCol) {
      this.subs.push(this.gridCol.grid.onAfterDrawContent.subscribe(() => {
        this._onRendered.next();
      }));
    }

    if (this.tableCol && this.gridCol) {
      this.subs.push(this.tableCol.dataFrame.onCurrentRowChanged.subscribe(() => {
        const df = this.tableCol.dataFrame;
        if (df.currentRowIdx === -1) {
          this.tableCol.temp[tempTAGS.referenceSequence] = null;
          this.tableCol.temp[tempTAGS.currentWord] = null;
          this.invalidateGrid();
        }
      }));
      const resetTriggerTags = [
        bioTAGS.positionShift,
        'renderMultiline',
        mmcrTAGS.showConnections,
      ];

      this.subs.push(DG.debounce(this.tableCol.dataFrame.onMetadataChanged.pipe(
        operators.filter((a) =>
          a.args.source === this.tableCol && resetTriggerTags.includes(a.args.key)
        )
      ), 200).subscribe((_) => {
        this.reset();
      }));

      this.subs.push(DG.debounce(this.tableCol.dataFrame.onMetadataChanged.pipe(
        operators.filter((a) => a.args.source === this.tableCol && a.args.key === bioTAGS.positionShift)
      ), 200).subscribe((_) => {
        this.reset();
      }));
    }
  }


  private calculateFontBasedSpacing(g: CanvasRenderingContext2D): {
    lineHeight: number;
    monomerSpacing: number;
  } {
    const metrics = g.measureText('M');

    // Get font size directly from the column's temp properties for safety.
    // This avoids parsing the font string and breaking the props interface.
    let fontSize = this.props?.fontSize ?? 12; // the package default, when the renderer's creator passes it
    if (this.tableCol?.temp[MmcrTemps.fontSize]) {
      const sizeFromCol = this.tableCol.temp[MmcrTemps.fontSize];
      if (typeof sizeFromCol === 'number' && !isNaN(sizeFromCol))
        fontSize = Math.max(sizeFromCol, 1);
    }

    // Line height is calculated based on the safe font size.
    const lineHeight = Math.max(fontSize * 1.4, metrics.fontBoundingBoxAscent + metrics.fontBoundingBoxDescent + 4);

    // Monomer spacing is proportional to character width.
    const monomerSpacing = Math.max(2, this.props.fontCharWidth * 0.2);

    return {lineHeight, monomerSpacing};
  }

  private shouldUseMultilineRendering(tableCol: DG.Column): boolean {
    const renderMultiline = tableCol.getTag('renderMultiline');
    return renderMultiline === 'true';
  }

  private calculateMultiLineLayoutDynamic(
    g: CanvasRenderingContext2D,
    w: number,
    h: number,
    subParts: ISeqSplitted,
    positionShift: number,
    maxLengthOfMonomer: number,
    canonical: boolean = false,
    extraLineH: number = 0,
  ): {
    lineLayouts: Array<{ lineIdx: number; elements: IMultilineLayoutElement[]; }>;
    lineHeight: number;
    } {
    if (this.dirty) {
      try {
        this.reset();
      } catch (err) {
        const [errMsg, errStack] = errInfo(err);
        this.logger.error(errMsg, undefined, errStack);
      }
    }
    // --- 1. Setup ---
    const {lineHeight: baseLineHeight, monomerSpacing} = this.calculateFontBasedSpacing(g);
    // room above each line for the connections' brackets or bond numbers
    const lineHeight = baseLineHeight + extraLineH;
    const availableWidth = w - (this.padding * 2);

    // --- 2. Find the widest monomer in the sequence to set a uniform column width ---
    let maxMonomerWidth = 0;
    const monomers: {text: string, posIdx: number, isSeparator: boolean}[] = [];
    if (subParts.length > 0) {
      // a '|' slot between chains, as the single line draws it
      const chainStarts = subParts.graphInfo?.disjointSeqStarts;
      let nextChain = 1;
      for (let i = positionShift; i < subParts.length; i++) {
        for (; chainStarts && nextChain < chainStarts.length && chainStarts[nextChain] <= i; nextChain++) {
          if (chainStarts[nextChain] === i && i > positionShift)
            monomers.push({text: '|', posIdx: i, isSeparator: true});
        }
        const om = canonical ? subParts.getCanonical(i) : subParts.getOriginal(i);
        const shortMon = this.props.monomerToShort(om, maxLengthOfMonomer);
        monomers.push({text: shortMon, posIdx: i, isSeparator: false});
        maxMonomerWidth = Math.max(maxMonomerWidth, (this._monomerLengthMap[shortMon] ??= g.measureText(shortMon)).width);
      }
    }

    if (monomers.length === 0)
      return {lineLayouts: [], lineHeight: lineHeight};

    // --- 3. Calculate how many uniform columns can fit ---
    const uniformColumnWidth = maxMonomerWidth;
    let colsPerLine = Math.floor((availableWidth + monomerSpacing) / (uniformColumnWidth + monomerSpacing));
    colsPerLine = Math.max(1, colsPerLine);

    // --- 4. Generate the final grid layout ---
    const availableHeightForLines = h - (this.padding * 2);
    const linesToRenderCount = Math.max(0, Math.floor(availableHeightForLines / lineHeight));
    const lineLayouts: Array<{ lineIdx: number, elements: IMultilineLayoutElement[] }> = [];

    let monomerIdx = 0;
    for (let lineIdx = 0; lineIdx < linesToRenderCount && monomerIdx < monomers.length; lineIdx++) {
      const elementsForLine: IMultilineLayoutElement[] = [];
      for (let colIdx = 0; colIdx < colsPerLine && monomerIdx < monomers.length; colIdx++) {
        const monomer = monomers[monomerIdx];
        const xPos = this.padding + colIdx * (uniformColumnWidth + monomerSpacing);

        elementsForLine.push({
          posIdx: monomer.isSeparator ? undefined : monomer.posIdx,
          x: xPos,
          width: uniformColumnWidth, // Use the uniform width for all slots
          om: monomer.text,
          isSeparator: monomer.isSeparator,
        });
        monomerIdx++;
      }
      lineLayouts.push({lineIdx: lineIdx, elements: elementsForLine});
    }

    return {lineLayouts, lineHeight};
  }

  public async init(): Promise<void> {
    await Promise.all([
      (async () => {
        this.seqHelper = await getSeqHelper();
        this.invalidateGrid();
      })(),
      (async () => {
        const libHelper = await getMonomerLibHelper();
        this.sysMonomerLib = libHelper.getMonomerLib();
      })(),
    ]);

    this.subs.push(this.sysMonomerLib!.onChanged.subscribe(() => {
      this.reset();
    }));

    this.subs.push(grok.events.onCustomEvent(MACROMOLECULE_HIGHLIGHT_EVENT_ID)
      .subscribe((args: MacromoleculeHighlightEventArgs) => this.handleHighlightEvent(args)));

    this.reset();
  }

  private handleHighlightEvent(args: MacromoleculeHighlightEventArgs): void {
    if (!args || !this.tableCol || !this.tableCol.dataFrame) return;
    if (args.columnName !== this.tableCol.name) return;
    const df = this.tableCol.dataFrame;
    if (args.tableId && df.id !== args.tableId) return;
    if (args.tableName && df.name !== args.tableName) return;

    const map = (this.tableCol.temp[MACROMOLECULE_HIGHLIGHT_TEMP] ??=
      new Map<number, MacromoleculeHighlightEntry>()) as Map<number, MacromoleculeHighlightEntry>;
    if (args.monomers == null || args.monomers.length === 0)
      map.delete(args.rowIdx);
    else {
      map.set(args.rowIdx, {
        monomers: args.monomers.slice(),
        fillColor: args.fillColor,
        strokeColor: args.strokeColor,
      });
    }
    this.invalidateGrid();
  }

  private getHighlightForRow(rowIdx: number | null): MacromoleculeHighlightEntry | null {
    if (rowIdx == null) return null;
    const map = this.tableCol.temp[MACROMOLECULE_HIGHLIGHT_TEMP] as
      Map<number, MacromoleculeHighlightEntry> | undefined;
    return map?.get(rowIdx) ?? null;
  }

  // Paints a translucent background behind the monomer at `posIdx` if this row has a highlight
  // entry that includes that position. Uses the entry's `fillColor` (alpha forced) or the
  // default fill color. Called from both single-line and multi-line render paths.
  private drawHighlightBackground(
    g: CanvasRenderingContext2D, entry: MacromoleculeHighlightEntry | null, posIdx: number,
    rectX: number, rectY: number, rectW: number, rectH: number,
  ): void {
    if (!entry || !entry.monomers || entry.monomers.length === 0) return;
    if (!entry.monomers.includes(posIdx)) return;
    const fillCss = macromoleculeHighlightColorToCss(entry.fillColor ?? DEFAULT_MACROMOLECULE_HIGHLIGHT_FILL);
    g.save();
    try {
      g.fillStyle = fillCss;
      g.fillRect(rectX, rectY, rectW, rectH);
    } finally {
      g.restore();
    }
  }

  /** The monomers' font: the column's font size temp, else {@link defaultFontSize} (the package default) */
  public static getFontSettings(tableCol?: DG.Column, defaultFontSize: number = 12) {
    const fontFamily = 'monospace';
    let fontSize = typeof defaultFontSize === 'number' && defaultFontSize > 0 ? defaultFontSize : 12;
    if (tableCol && tableCol.temp[MmcrTemps.fontSize] &&
        typeof tableCol.temp[MmcrTemps.fontSize] === 'number' && !isNaN(tableCol.temp[MmcrTemps.fontSize]))
      fontSize = Math.max(tableCol.temp[MmcrTemps.fontSize], 1);
    // monospace ratio
    const fontWidth = fontSize * 0.6;
    return {font: `${fontSize}px ${fontFamily}`, fontWidth, fontSize};
  }

  override toLog(): string {
    return `MonomerPlacer<${this.viewerId}>`;
  }

  protected getMonomerLib(): IMonomerLibBase | null {
    return this.tableCol.temp[MmcrTemps.overriddenLibrary] ?? this.sysMonomerLib;
  }

  protected override reset(): void {
    if (this.propsProvider)
      this.props = this.propsProvider();
    this._processedRows = DG.BitSet.create(this.tableCol!.length);
    this._monomerLengthList = null;
    this._monomerLengthMap = {};
    this._monomerStructureMap = {};
    this._cellBounds.clear();
    this._connLayouts.clear();
    this._stackedBounds.clear();
    this._monomerLengthCanonical = [];
    this._darkerColors.clear();
    this._showConnections = this.resolveShowConnections() && this.columnHasGraphInfo();
    super.reset();
    this.invalidateGrid();
  }

  protected invalidateGrid(): void {
    if (this.gridCol && this.gridCol.dart) this.gridCol.grid?.invalidate();
  }

  /** Returns monomers lengths of the {@link rowIdx} and cumulative sums for borders, monomer places */
  public getCellMonomerLengths(rowIdx: number, newWidth: number): [number[], number[]] {
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    if (this.colWidth < newWidth) {
      this.colWidth = newWidth;
      this.dirty = true;
    }
    if (this.dirty) {
      try { this.reset(); } catch (err) {
        const [errMsg, errStack] = errInfo(err);
        this.logger.error(errMsg, undefined, errStack);
      }
    }

    const res: number[] = sh.isMsa() ? this.getCellMonomerLengthsForSeqMsa() :
      this.getCellMonomerLengthsForSeq(rowIdx);

    return [res, this.getSummedMonomerLengths(res)];
  }

  private getSummedMonomerLengths(res: number[]): number[] {
    const resSum: number[] = new Array<number>(res.length + 1);
    resSum[0] = this.padding; // padding
    for (let pos: number = 1; pos < resSum.length; pos++)
      resSum[pos] = resSum[pos - 1] + res[pos - 1];
    // due to implementation specifics and performance, resSum can have NaN s at the end for stuff that is not visible
    let lastNumValue = resSum[0];
    for (let pos = 1; pos < resSum.length; pos++) {
      if (!resSum[pos])// 0 values will def not be here, so ! check is good enaugh
        resSum[pos] = lastNumValue;
      else
        lastNumValue = resSum[pos];
    }
    return resSum;
  }

  private getCellMonomerLengthsForSeqValue(value: string, width: number, canonical: boolean = false): number[] {
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
    const visibleSeqStart = this.positionShift;
    const maxVisibleSeqLength: number = Math.ceil(width / minMonWidth) + visibleSeqStart;
    const seqSS: ISeqSplitted = sh.splitter(value);
    const visibleSeqEnd: number = Math.min(maxVisibleSeqLength, seqSS.length);
    const res = new Array<number>(visibleSeqEnd - visibleSeqStart);
    let seqWidth: number = 0;
    for (let seqMonI = visibleSeqStart; seqMonI < visibleSeqEnd; ++seqMonI) {
      const seqMonLabel = this.monomerLabel(seqSS, seqMonI, canonical);
      const shortMon: string = this.props.monomerToShort(seqMonLabel, this.monomerLengthLimit);
      const separatorWidth = sh.isSeparator() ? this.separatorWidth : this.props.separatorWidth;
      const seqMonWidth: number = separatorWidth + shortMon.length * this.props.fontCharWidth;
      res[seqMonI - visibleSeqStart] = seqMonWidth;
      seqWidth += seqMonWidth;
      if (seqWidth > width) break;
    }
    return res;
  }

  private getCellMonomerLengthsForSeq(rowIdx: number): number[] {
    const logPrefix = `${this.toLog()}.getCellMonomerLengthsForSeq()`;
    // this.logger.debug(`${logPrefix}, start`);

    if (this._monomerLengthList == null)
      this._monomerLengthList = new Array(this.tableCol.length).fill(null);

    const visibleSeqStart = this.positionShift;
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
    const maxVisibleSeqLength: number = Math.ceil(this.colWidth / minMonWidth) + visibleSeqStart;
    const seqSS: ISeqSplitted = sh.getSplitted(rowIdx);
    const visibleSeqEnd: number = Math.min(maxVisibleSeqLength, seqSS.length);

    let res: number[] = this._monomerLengthList[rowIdx];
    const canonical = this._showConnections &&
      this.connectionsLayout(rowIdx, seqSS, this._widthsH, visibleSeqStart) !== null;
    if (res == null || res.length != visibleSeqEnd - visibleSeqStart || this._monomerLengthCanonical[rowIdx] !== canonical) {
      res = this._monomerLengthList[rowIdx] = new Array<number>(seqSS.length);
      this._monomerLengthCanonical[rowIdx] = canonical;

      let seqWidth: number = 0;
      for (let seqMonI = visibleSeqStart; seqMonI < visibleSeqEnd; ++seqMonI) {
        const seqMonLabel = this.monomerLabel(seqSS, seqMonI, canonical);
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.monomerLengthLimit);
        const separatorWidth = sh.isSeparator() ? this.separatorWidth : this.props.separatorWidth;
        const seqMonWidth: number = separatorWidth + shortMon.length * this.props.fontCharWidth;
        res[seqMonI - visibleSeqStart] = seqMonWidth;
        seqWidth += seqMonWidth;
        if (seqWidth > this.colWidth) break;
      }
    }
    return res;
  }

  private getCellMonomerLengthsForSeqMsa(): number[] {
    const logPrefix = `${this.toLog()}.getCellMonomerLengthsForSeqMsa()`;
    // this.logger.debug(`${logPrefix}, start`);
    if (this._showConnections && this._msaWidthsH !== this._widthsH) {
      this._msaWidthsH = this._widthsH;
      this._monomerLengthList = null;
      this._processedRows = DG.BitSet.create(this.tableCol.length);
      this._processedMaxVisibleSeqLength = 0;
    }
    if (this._monomerLengthList == null)
      this._monomerLengthList = new Array(1).fill(null);
    this._monomerLengthList[0] ??= new Array(0);
    const res = this._monomerLengthList[0];
    // const startIdx = Math.max(Math.floor((this.grid?.vertScroll.min ?? 0) - 10), 0);
    // const endIdx = Math.min(Math.ceil((this.grid?.vertScroll.max ?? 0) + 10), this.col.length);

    const {startIdx, endIdx} = (() => {
      try {
        const grid: DG.Grid | null = this.gridCol && this.gridCol.dart ? this.gridCol.grid : null;
        if (grid && grid.dart) {
          return {
            startIdx: Math.max(Math.floor((grid?.vertScroll.min ?? 0) - 10), 0),
            endIdx: Math.min(Math.ceil((grid?.vertScroll.max ?? 0) + 10), this.tableCol.length)
          };
        } else return {startIdx: 0, endIdx: Math.min(this.tableCol.length, 10)};
      } catch (_e) {
        return {startIdx: 0, endIdx: Math.min(this.tableCol.length, 10)};
      }
    })();

    const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
    const visibleSequenceStart = this.positionShift;
    const maxVisibleSeqLength: number = Math.ceil(this.colWidth / minMonWidth) + visibleSequenceStart;
    for (let seqIdx = startIdx; seqIdx < endIdx; seqIdx++) {
      if (this._processedRows.get(seqIdx) && maxVisibleSeqLength <= this._processedMaxVisibleSeqLength)
        continue;
      const sh = this.seqHelper.getSeqHandler(this.tableCol);
      const seqSS: ISeqSplitted = sh.getSplitted(seqIdx, maxVisibleSeqLength);
      // a row that draws its connections prints canonical symbols; its graph info needs the whole sequence
      const canonical = this._showConnections &&
        this.connectionsLayout(seqIdx, sh.getSplitted(seqIdx), this._widthsH, visibleSequenceStart) !== null;

      const visibleSeqEnd: number = Math.min(maxVisibleSeqLength, seqSS.length);
      if (visibleSeqEnd - visibleSequenceStart > res.length)
        res.push(...new Array<number>(visibleSeqEnd - visibleSequenceStart - res.length).fill(minMonWidth));

      let seqWidth: number = 0;
      for (let seqMonI = visibleSequenceStart; seqMonI < visibleSeqEnd; ++seqMonI) {
        const seqMonLabel: string = this.monomerLabel(seqSS, seqMonI, canonical);
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.monomerLengthLimit);
        const seqMonWidth: number = this.props.separatorWidth + shortMon.length * this.props.fontCharWidth;
        res[seqMonI - visibleSequenceStart] = Math.max(res[seqMonI - visibleSequenceStart] ?? 0, seqMonWidth);
        seqWidth += seqMonWidth;
        if (seqWidth >= this.colWidth) break;
      }
      this._processedMaxVisibleSeqLength = Math.max(this._processedMaxVisibleSeqLength, maxVisibleSeqLength);
      this._processedRows.set(seqIdx, true);
    }
    return res; // first (and single) row of data
  }

  /** Returns seq position for pointer x */
  public getPosition(rowIdx: number, x: number, width: number, positionShiftPadding?: number): number | null {
    const [_monomerMaxLengthList, monomerMaxLengthSumList]: [number[], number[]] =
      this.getCellMonomerLengths(rowIdx, width);

    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const seqSS = sh.getSplitted(rowIdx);
    if (seqSS.length === 0) return null;

    return hitBounds(monomerMaxLengthSumList, x, positionShiftPadding);
  }

  /** Where the monomer at sequence position {@link pos} (0-based) of the cell's row sits in the cell, in px from its
   * top left, as {@link onMouseMove} hit-tests it: the point the pointer rests on for its tooltip. Null when the cell
   * does not draw it. */
  public getMonomerPoint(gridCell: DG.GridCell, pos: number): {x: number, y: number} | null {
    const rowIdx = gridCell.tableRowIndex;
    const idx = pos - this.positionShift;
    if (!this.seqHelper || rowIdx == null || idx < 0)
      return null;
    const {width: cellWidth, height: cellHeight} = gridCell.bounds;
    const bounds = this._cellBounds.get(rowIdx) ?? this._stackedBounds.get(cellHeight)?.get(rowIdx);
    if (bounds) {
      const b = bounds.find((bound) => bound.monomerIdx === idx)?.bounds;
      return b ? {x: b.x + b.width / 2, y: b.y + b.height / 2} : null;
    }
    if (this._showConnections)
      this._widthsH = cellHeight;
    const [lengths, sums] = this.getCellMonomerLengths(rowIdx, cellWidth);
    const leftPadding = this.shouldRenderShiftedThreeDots(this.positionShift) ? this._leftThreeDotsPadding ?? 0 : 0;
    const x = leftPadding + sums[idx] + lengths[idx] / 2;
    return idx < lengths.length && lengths[idx] > 0 && x < cellWidth ? {x, y: cellHeight / 2} : null;
  }

  public setMonomerLengthLimit(value: number): void {
    if (this.monomerLengthLimit != value) {
      this.monomerLengthLimit = value;
      this.dirty = true;
    }
  }

  public setSeparatorWidth(value: number): void {
    if (this.separatorWidth != value) {
      this.props.separatorWidth = value;
      this.dirty = true;
    }
  }

  private padding: number = 5;

  get positionShift(): number {
    const parsed = Number.parseInt(this.tableCol?.tags[bioTAGS.positionShift] ?? '0') ?? 0;
    return isNaN(parsed) ? 0 : Math.max(parsed, 0);
  }

  // --- Connections drawn as lines (cell-renderer-connections.ts) ---

  /** The column's `.mm.cellRenderer.showConnections` tag ('true' / 'false'), else the package default from the props */
  private resolveShowConnections(): boolean {
    const tag = this.tableCol?.getTag(mmcrTAGS.showConnections);
    if (tag === 'true') return true;
    if (tag === 'false') return false;
    return this.props?.showConnections ?? false;
  }

  /** Whether the column's rows carry graph info (BILN, HELM, a custom notation's splitter): only then can they draw
   * connections. The other columns take none of the connection steps */
  private columnHasGraphInfo(): boolean {
    if (!this.seqHelper || !this.tableCol || this.tableCol.length === 0) return false;
    try {
      return this.seqHelper.getSeqHandler(this.tableCol).getSplitted(0).graphInfo !== undefined;
    } catch {
      return false;
    }
  }

  /** The label a position prints: the canonical symbol in a row that draws its connections, else the original */
  private monomerLabel(seqSS: ISeqSplitted, posIdx: number, canonical: boolean): string {
    return canonical ? seqSS.getCanonical(posIdx) : seqSS.getOriginal(posIdx);
  }

  private fontSizeOf(): number {
    if (this.props?.fontSize) return this.props.fontSize;
    const v = this.tableCol?.temp[MmcrTemps.fontSize];
    return typeof v === 'number' && !isNaN(v) ? Math.max(v, 1) : 12;
  }

  /** How the row's connections fit a cell of height {@link h}, or null (the setting is off, the row has none, they do
   * not fit): cached per row and height for a grid row ({@link rowIdx} not null) */
  private connectionsLayout(rowIdx: number | null, seqSS: ISeqSplitted, h: number, positionShift: number,
  ): IConnectionsLayout | null {
    if (!this._showConnections || h <= 0) return null;
    const gi = seqSS.graphInfo;
    if (!gi || (!gi.connections?.length && !(gi as ISeqGraphInfoHelm).hbonds?.length)) return null;
    let byRow: Map<number, IConnectionsLayout | null> | undefined;
    if (rowIdx !== null) {
      byRow = this._connLayouts.get(h);
      if (!byRow || byRow.size >= 20000) {
        if (!byRow && this._connLayouts.size >= 4) this._connLayouts.clear();
        this._connLayouts.set(h, byRow = new Map());
      }
      const cached = byRow.get(rowIdx);
      if (cached !== undefined) return cached;
    }
    const fontSize = this.fontSizeOf();
    // one line per chain only for a sequence shown from its start
    const res = layoutConnections(seqSS, h, fontSize, measureConnectionsFont(this.props.font, fontSize), positionShift === 0);
    byRow?.set(rowIdx!, res);
    return res;
  }

  /** A connection's colour: its monomers' colour one shade darker when both ends are the same monomer (a Cys-Cys
   * bridge), grey otherwise and for rings and combs */
  private connectionColor(seqSS: ISeqSplitted, monomerColor: (cm: string) => string, colorCode: boolean,
  ): (c: IDrawnConnection) => string {
    const darker = this._darkerColors;
    return (c: IDrawnConnection): string => {
      if (!colorCode || c.headToTail || c.comb) return CONNECTION_GREY;
      const cm = seqSS.getCanonical(c.from);
      if (cm !== seqSS.getCanonical(c.to)) return CONNECTION_GREY;
      let color = darker.get(cm);
      if (color === undefined) darker.set(cm, color = darkerColor(monomerColor(cm)));
      return color;
    };
  }

  /** The positions bonded to the hovered monomer of this row (empty when none is hovered or the setting is off). The
   * grid of {@link gridCell} (null off the grid) is read only for the row with the hovered monomer. */
  private hoverPartners(rowIdx: number | null, seqSS: ISeqSplitted, gridCell: DG.GridCell | null,
  ): { pos: number | null, partners: Set<number> } {
    const hover = this._hoverConn;
    if (!this._showConnections || !hover || hover.row !== rowIdx || hover.grid?.dart !== gridCell?.grid?.dart)
      return NO_HOVER;
    const partners = new Set<number>();
    for (const c of getDrawnConnections(seqSS)) {
      if (c.from === hover.pos) partners.add(c.to);
      else if (c.to === hover.pos) partners.add(c.from);
    }
    return {pos: partners.size > 0 ? hover.pos : null, partners};
  }

  /** Compares the settings that change without an event, once per render pass (the grid draws its cells in one task,
   * so its first cell does it): the column's font size, which a grid's panel sets for every grid that shows the column
   * while the first renderer to draw clears the changed flag, and the package defaults */
  private checkSettingsOncePerPass(): void {
    if (this._passChecked) return;
    this._passChecked = true;
    queueMicrotask(() => this._passChecked = false);
    const fontSizeTemp = this.tableCol.temp[MmcrTemps.fontSize];
    if (!Object.is(fontSizeTemp, this._appliedFontSizeTemp)) {
      this._appliedFontSizeTemp = fontSizeTemp;
      this.dirty = true;
    }
    this.checkPackageDefaults();
  }

  /** The package defaults (font size, max monomer length, connections) the props provider reads change without an
   * event: compared with the props once a second, a change lays the cells out anew */
  private checkPackageDefaults(): void {
    const now = performance.now();
    if (now - this._propsCheckedAt < 1000 || !this.propsProvider) return;
    this._propsCheckedAt = now;
    try {
      const p = this.propsProvider();
      const c = this.props;
      if (p.font !== c.font || p.maxMonomerLength !== c.maxMonomerLength || p.showConnections !== c.showConnections)
        this.dirty = true;
    } catch (err) {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(errMsg, undefined, errStack);
    }
  }

  /** A soft box behind a monomer bonded to the hovered one */
  private drawPartnerBox(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number): void {
    g.save();
    try {
      g.fillStyle = 'rgba(255, 196, 0, 0.28)';
      g.fillRect(x, y, w, h);
    } finally {
      g.restore();
    }
  }

  private _leftThreeDotsPadding: number = 0;

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle,
  ) {
    // for cases when we render it on somewhere other than grid, gridRow might be null or incorrect (set to 0).
    //for this case we can just recalculate split sequence without caching
    const isRenderedOnGrid = gridCell.grid?.dart && gridCell.grid?.canvas === g.canvas;

    if (!this.seqHelper) return;
    const tableCol = this.tableCol;
    const positionShift = this.positionShift;

    g.save();
    try {
      const sh = this.seqHelper.getSeqHandler(tableCol);
      // the column's max monomer length below, else the package default (when the renderer's creator passes it)
      let maxLengthOfMonomer: number = this.props.maxMonomerLength ?? this.monomerLengthLimit;
      if (mmcrTAGS.maxMonomerLength in tableCol.tags) {
        const v = parseInt(tableCol.getTag(mmcrTAGS.maxMonomerLength));
        maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
      }

      if (MmcrTemps.maxMonomerLength in tableCol.temp) {
        const v = tableCol.temp[MmcrTemps.maxMonomerLength];
        const vn = typeof v === 'number' ? v : parseInt(v);
        maxLengthOfMonomer = !isNaN(vn) && vn ? vn : 50;
      }

      if (
        tableCol.temp[MmcrTemps.rendererSettingsChanged] === rendererSettingsChangedState.true ||
        this.monomerLengthLimit != maxLengthOfMonomer
      ) {
        // this if means that the mm renderer settings have changed,
        // particularly monomer representation and max width.
        let gapLength = 0;
        const msaGapLength = 8;
        gapLength = tableCol.temp[MmcrTemps.gapLength] as number ?? gapLength;
        this.setMonomerLengthLimit(maxLengthOfMonomer);
        this.setSeparatorWidth(sh.isMsa() ? msaGapLength : gapLength);
        tableCol.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.false;
        this.dirty = true;
      }
      this.checkSettingsOncePerPass();
      // The reset reads the props anew (the font among them) and empties the monomers' text metrics: before the font is
      // set, so this cell prints and measures in the new font. After it, this cell printed in the old font and its
      // monomers kept the old metrics, by which every later cell placed them above or below the others.
      if (this.dirty) {
        try { this.reset(); } catch (err) {
          const [errMsg, errStack] = errInfo(err);
          this.logger.error(errMsg, undefined, errStack);
        }
      }


      const rowIdx = gridCell.cell.rowIndex;
      const value: any = gridCell.cell.value;
      // The cell's own width: the multiline layout wraps at the cell's bounds, so its lines stay put while the grid
      // scrolls horizontally. Below, `w` is cut to the grid's visible edge, which is enough for the single line.
      const cellW = w;
      if (isRenderedOnGrid)
        w = getUpdatedWidth(gridCell.grid, g, x, w, window.devicePixelRatio);

      g.beginPath();
      g.rect(x, y, w, h);
      g.clip();
      g.font = this.props?.font ?? '12px monospace';
      g.textBaseline = 'top';

      const units = tableCol.meta.units;
      const aligned: string = tableCol.getTag(bioTAGS.aligned);
      const separator = tableCol.getTag(bioTAGS.separator) ?? '';
      const subParts: ISeqSplitted = isRenderedOnGrid ? sh.getSplitted(rowIdx) : sh.splitter(value);

      let drawStyle = DrawStyle.classic;
      if (aligned?.includes('MSA') && units === NOTATION.SEPARATOR)
        drawStyle = DrawStyle.MSA;

      const tempReferenceSequence: string | null = tableCol.temp[tempTAGS.referenceSequence];
      const tempCurrentWord: string | null = this.tableCol.temp[tempTAGS.currentWord];
      const referenceSequence: string[] = (() => {
        const splitterFunc: SplitterFunc = sh.splitter;
        const seqSS = splitterFunc(
          ((tempReferenceSequence != null) && (tempReferenceSequence !== '')) ?
            tempReferenceSequence : tempCurrentWord ?? '');
        return wu.count(0).take(seqSS.length).slice(positionShift).map((posIdx) => seqSS.getCanonical(posIdx)).toArray();
      })();
      const selectedPosition = Number.parseInt(tableCol.getTag(bioTAGS.selectedPosition) ?? '-200');
      // Read once per cell, not per monomer (each is a call into the platform): the row, its highlight entry and, at
      // the first monomer that asks, whether it is the current row (they do not change mid-render)
      const tableRowIdx = gridCell.tableRowIndex;
      const highlightEntry = this.getHighlightForRow(tableRowIdx);
      let notCurrentRow: boolean | undefined = undefined;
      // A monomer's text colour, once per monomer symbol in this cell (the library does not change mid-render)
      const textColors = new Map<string, string>();
      const monomerTextColor = (cm: string): string => {
        let color = textColors.get(cm);
        if (color === undefined) {
          color = undefinedColor;
          const lib = this.getMonomerLib();
          if (lib)
            color = lib.getMonomerTextColor(sh.defaultBiotype, cm);
          textColors.set(cm, color);
        }
        return color;
      };


      const shouldUseMultiLine = this.shouldUseMultilineRendering(tableCol);

      const annotRenderer = this.getAnnotationRenderer();

      if (shouldUseMultiLine) {
        const currentCellBounds: IMonomerLayoutData[] = [];
        // Connections in a wrapped sequence: canonical symbols; a bond whose two ends share a line is bracketed over it,
        // a bond a line break cuts gets the same small number at both ends
        const mlConns = this._showConnections ? getDrawnConnections(subParts) : [];
        let ml = mlConns.length > 0;
        let layout = this.calculateMultiLineLayoutDynamic(g, cellW, h, subParts, positionShift, maxLengthOfMonomer, ml);
        let extraH = 0;
        let mlFontSize = 0;
        let mlGeo: IConnectionsGeometry | null = null;
        let lineOf: Map<number, number> | null = null;
        let bracketsByLine: Map<number, IDrawnConnection[]> | null = null;
        const numbered: IDrawnConnection[] = [];
        if (ml) {
          lineOf = new Map<number, number>();
          for (const ll of layout.lineLayouts)
            for (const el of ll.elements) if (!el.isSeparator) lineOf.set(el.posIdx!, ll.lineIdx);
          bracketsByLine = new Map<number, IDrawnConnection[]>();
          for (const c of mlConns) {
            const line = lineOf.get(c.from);
            if (line !== undefined && line === lineOf.get(c.to)) {
              let cs = bracketsByLine.get(line);
              if (!cs) bracketsByLine.set(line, cs = []);
              cs.push(c);
            } else
              numbered.push(c);
          }
          mlFontSize = this.fontSizeOf();
          mlGeo = connectionsGeometry(mlFontSize, false);
          let lanes = 0;
          for (const cs of bracketsByLine.values()) {
            lanes = Math.max(lanes, assignLanes(cs, () => 1).above);
            cs.sort((a, b) => b.lane - a.lane); // outer lanes first
          }
          // the room over every line, as much as the line that needs the most takes (the lines keep one height); a
          // monomer's line does not depend on it, only how many lines the cell has room for
          extraH = Math.max(numbered.length > 0 ? bondNumberFontSize(mlFontSize) + 1 : 0, sideHeight(lanes, mlGeo));
          if (Math.floor((h - this.padding * 2) / (layout.lineHeight + extraH)) >= layout.lineLayouts.length)
            layout.lineHeight += extraH;
          else {
            // the room would cost lines of the sequence: printed as before, with the marks
            ml = false;
            extraH = 0;
            layout = this.calculateMultiLineLayoutDynamic(g, cellW, h, subParts, positionShift, maxLengthOfMonomer);
          }
        }
        const {pos: hoverPos, partners} = this.hoverPartners(tableRowIdx, subParts, isRenderedOnGrid ? gridCell : null);
        const centerX = ml ? new Map<number, number>() : null;
        const lineYs = ml ? new Map<number, number>() : null;

        // --- NEW: Vertical Centering Logic for Single Lines ---
        // Default to top-aligned layout
        let yBase = y + this.padding;

        // If there's exactly one line of content, calculate a new base Y to center it vertically.
        if (layout.lineLayouts.length === 1)
          yBase = y + (h - layout.lineHeight) / 2;

        // --- End of New Logic ---

        for (const lineLayout of layout.lineLayouts) {
          // The Y position for each line is now based on our calculated `yBase`.
          const lineY = yBase + (lineLayout.lineIdx * layout.lineHeight);
          lineYs?.set(lineLayout.lineIdx, lineY);

          for (const element of lineLayout.elements) {
            const elementX = x + element.x;
            if (element.isSeparator) { // a chain break
              printLeftOrCentered(g, element.om, elementX, lineY + extraH, element.width, layout.lineHeight - extraH,
                {color: 'rgb(128,128,128)', isMultiLineContext: true, wordIdx: -1});
              continue;
            }
            const monomer = element;
            const monomerIndex = monomer.posIdx!;
            const cm = subParts.getCanonical(monomerIndex);

            const color = monomerTextColor(cm);

            let transparencyRate = 0.0;
            if ((notCurrentRow ??= tableRowIdx !== tableCol.dataFrame.currentRowIdx) && referenceSequence.length > 0) {
              const refIndex = monomerIndex - positionShift;
              if (refIndex >= 0 && refIndex < referenceSequence.length) {
                const currentMonomerCanonical = cm;
                const refMonomerCanonical = referenceSequence[refIndex];
                if (currentMonomerCanonical === refMonomerCanonical)
                  transparencyRate = 0.7;
              }
            }
            currentCellBounds.push({
              lineIdx: lineLayout.lineIdx,
              monomerIdx: monomerIndex - positionShift,
              bounds: new DG.Rect(element.x, (lineY - y), element.width, layout.lineHeight),
              sequencePosition: monomerIndex,
            });

            if (annotRenderer)
              annotRenderer.drawPositionBackground(g, elementX, lineY, element.width + 2, layout.lineHeight, monomerIndex, tableRowIdx!);

            this.drawHighlightBackground(
              g, highlightEntry, monomerIndex,
              elementX, lineY, element.width + 2, layout.lineHeight);
            if (hoverPos !== null && (partners.has(monomerIndex) || monomerIndex === hoverPos))
              this.drawPartnerBox(g, elementX, lineY + extraH, element.width + 2, layout.lineHeight - extraH);
            centerX?.set(monomerIndex, elementX + element.width / 2);

            printLeftOrCentered(g, monomer.om, elementX, lineY + extraH, element.width, layout.lineHeight - extraH, {
              color: color,
              isMultiLineContext: true,
              transparencyRate: transparencyRate,
              selectedPosition: isNaN(selectedPosition) || selectedPosition < 1 ? undefined : selectedPosition,
              wordIdx: monomerIndex,
            });
          }
        }
        if (ml && lineYs && centerX && lineOf && bracketsByLine && mlGeo) {
          const met = measureConnectionsFont(this.props.font, mlFontSize);
          const capH = met.capBotOff - met.capTopOff;
          // the letters are centred in the band under the extra room
          const capTopOf = (line: number): number =>
            lineYs.get(line)! + extraH + (layout.lineHeight - extraH) / 2 - capH / 2 - 0.5;
          const colorOf = this.connectionColor(subParts, monomerTextColor, getPrintColumnTemps(gridCell).colorCode);
          for (const [line, cs] of bracketsByLine) {
            if (lineYs.has(line)) // a line the cell has room for
              drawConnections(g, cs, (p: number) => centerX.get(p) ?? x - 20, capTopOf(line) - 1.5, 0, mlGeo, colorOf, hoverPos);
          }
          // the cut bonds numbered in the order of their first monomer
          const numbersAt = new Map<number, number[]>();
          numbered.sort((a, b) => a.from - b.from).forEach((c, i) => {
            for (const p of [c.from, c.to]) {
              let l = numbersAt.get(p);
              if (!l) numbersAt.set(p, l = []);
              l.push(i + 1);
            }
          });
          const items: { x: number, y: number, numbers: number[] }[] = [];
          for (const [p, numbers] of numbersAt) {
            const line = lineOf.get(p);
            if (line !== undefined && lineYs.has(line) && centerX.has(p))
              items.push({x: centerX.get(p)!, y: capTopOf(line) - 1, numbers});
          }
          drawBondNumbers(g, items, mlFontSize);
        }
        if (isRenderedOnGrid && gridCell.tableRowIndex !== null)
          this._cellBounds.set(gridCell.tableRowIndex, currentCellBounds);
      } else {
        // --- Single-line rendering path ---
        // the current row is the one the others are compared with: it is drawn whole, as the multiline path draws it
        const rowReference = referenceSequence.length > 0 && !(notCurrentRow ??= tableRowIdx !== tableCol.dataFrame.currentRowIdx) ?
          [] : referenceSequence;
        // The row's connections drawn as lines, when the setting is on and they fit the cell: such a row prints the
        // canonical symbols, in the widths and hit tests too, which follow the height the cell was drawn at
        // the height a row's widths are decided at; each row's cached widths follow its own canonical state, so another
        // height (a heatmap, a second grid of this BILN column) does not lay out the others anew
        if (this._showConnections && isRenderedOnGrid)
          this._widthsH = h;
        const conn = this.connectionsLayout(isRenderedOnGrid ? tableRowIdx : null, subParts, h, positionShift);
        if (conn !== null && conn.chains.length > 0) {
          this.renderStackedChains(g, x, y, w, h, gridCell, {sh, seqSS: subParts, layout: conn, maxLengthOfMonomer,
            separator, monomerColor: monomerTextColor, referenceSequence: rowReference, selectedPosition, highlightEntry, annotRenderer,
            rowIdx: isRenderedOnGrid ? tableRowIdx : null, positionShift});
          return;
        }
        this._leftThreeDotsPadding = this.shouldRenderShiftedThreeDots(positionShift) ? g.measureText(shiftedLeftPaddingText).width : 0;
        let [, maxLengthWordsSum]: [number[], number[]] = this.getCellMonomerLengths(gridCell.tableRowIndex!, w);
        if (!isRenderedOnGrid)
          maxLengthWordsSum = this.getSummedMonomerLengths(this.getCellMonomerLengthsForSeqValue(value, w, conn !== null));
        const minMonomerWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
        const visibleSeqLength = Math.min(subParts.length, Math.ceil(w / (minMonomerWidth)) + positionShift);

        // The monomers' band: the whole cell, or the part the connections' lanes leave above and below
        let bandY = y;
        let bandH = h;
        let aboveY0 = 0;
        let belowY0 = 0;
        if (conn) {
          const met = measureConnectionsFont(this.props.font, this.fontSizeOf());
          const aboveH = sideHeight(conn.above, conn.geo);
          const belowH = sideHeight(conn.below, conn.geo);
          const capH = met.capBotOff - met.capTopOff;
          const capTop = y + (h - (aboveH + capH + Math.max(belowH + 2, met.descBelow))) / 2 + aboveH;
          const drawY = capTop - met.capTopOff;
          // printLeftOrCentered prints 1 px under the top of a band of the em height
          bandY = drawY - 1;
          bandH = met.emH;
          aboveY0 = capTop - 1.5;
          belowY0 = drawY + met.capBotOff + 2;
        }

        // Pre-compute text bottom y for annotation underlines (text is vertically centered)
        let textBottomY: number | undefined;
        // Compute a monomer-sized highlight rect (a little taller than the glyph itself) so the
        // background doesn't span the whole cell height.
        // the column temps every printed monomer of this cell reads, once
        const columnTemps = getPrintColumnTemps(gridCell);
        const monomerMetrics = g.measureText('M');
        const monomerTextHeight = monomerMetrics.fontBoundingBoxAscent + monomerMetrics.fontBoundingBoxDescent;
        const hlPadding = 4;
        const hlRectHeight = Math.min(h, monomerTextHeight + hlPadding);
        const hlRectY = bandY + (bandH - hlRectHeight) / 2;
        if (annotRenderer) {
          const dy = bandH / 2 - monomerTextHeight / 2 + 1;
          textBottomY = bandY + dy + monomerTextHeight;
        }
        const {pos: hoverPos, partners} = this.hoverPartners(tableRowIdx, subParts, isRenderedOnGrid ? gridCell : null);

        for (let posIdx: number = positionShift; posIdx < visibleSeqLength; ++posIdx) {
          const cm: string = posIdx < subParts.length ? subParts.getCanonical(posIdx) : sh.defaultGapOriginal;
          const om: string = conn ? cm : posIdx < subParts.length ? subParts.getOriginal(posIdx) : sh.defaultGapOriginal;

          const color = monomerTextColor(cm);

          const last = posIdx === subParts.length - 1;
          // here, the first index of disjointSeqStarts is always 0, so we check from position 1, and also, next index
          // will show where the next disjoint sequence starts, so we check posIdx + 1
          const drawnSeparator = (subParts?.graphInfo?.disjointSeqStarts?.indexOf(posIdx + 1) ?? 0) > 0 ? '|' : separator;

          // Draw annotation background before the monomer text
          const wordIdx = posIdx - positionShift;
          const placeX = (maxLengthWordsSum[wordIdx] ?? 0) - (maxLengthWordsSum[0] ?? 0);
          const nextPlaceX = wordIdx + 1 < maxLengthWordsSum.length ?
            (maxLengthWordsSum[wordIdx + 1] ?? 0) - (maxLengthWordsSum[0] ?? 0) : placeX + this.props.fontCharWidth;
          const monomerW = nextPlaceX - placeX;
          const monomerLeftX = x + this.padding + this._leftThreeDotsPadding + placeX;
          if (annotRenderer && tableRowIdx != null) {
            annotRenderer.drawPositionBackground(
              g, monomerLeftX, y, monomerW, h,
              posIdx, tableRowIdx, textBottomY,
            );
          }

          this.drawHighlightBackground(g, highlightEntry, posIdx,
            monomerLeftX, hlRectY, monomerW, hlRectHeight);
          if (hoverPos !== null && (partners.has(posIdx) || posIdx === hoverPos)) {
            const short = this.props.monomerToShort(om, maxLengthOfMonomer);
            const labelW = (this._monomerLengthMap[short] ??= g.measureText(short)).width;
            this.drawPartnerBox(g, monomerLeftX - 1, hlRectY, labelW + 2, hlRectHeight);
          }

          const opts: Partial<PrintOptions> = {
            color: color, pivot: 0, left: true, transparencyRate: 0.0,
            separator: drawnSeparator, last: last,
            drawStyle: drawStyle, maxWord: maxLengthWordsSum, wordIdx: posIdx - positionShift, gridCell: gridCell,
            referenceSequence: rowReference, maxLengthOfMonomer: maxLengthOfMonomer,
            monomerTextSizeMap: this._monomerLengthMap, logger: this.logger,
            selectedPosition: isNaN(selectedPosition) || selectedPosition < 1 ? undefined : selectedPosition - positionShift,
            columnTemps,
          };
          printLeftOrCentered(g, om, x + this.padding + this._leftThreeDotsPadding, bandY, w, bandH, opts);
        }
        if (this.shouldRenderShiftedThreeDots(positionShift)) {
          const opts = {
            color: undefinedColor, pivot: 0, left: true, transparencyRate: 0, separator: separator, last: false,
            drawStyle: drawStyle, maxWord: maxLengthWordsSum, wordIdx: 0, gridCell: gridCell,
            maxLengthOfMonomer: maxLengthOfMonomer,
            monomerTextSizeMap: this._monomerLengthMap, logger: this.logger, columnTemps,
          };
          printLeftOrCentered(g, shiftedLeftPaddingText, x + this.padding, bandY, w, bandH, opts);
        }
        if (conn) {
          // a monomer's x for a bond's R-group: R1 at its left edge, R2 at its right edge, a side chain at its centre;
          // monomers scrolled out or past the cell's edge lie beyond the clip, so their brackets run off the cell
          const x0 = x + this.padding + this._leftThreeDotsPadding;
          const s0 = maxLengthWordsSum[0] ?? 0;
          const xOf = (pos: number, r?: number): number => {
            const wordIdx = pos - positionShift;
            if (wordIdx < 0) return x - 20;
            const left = x0 + (maxLengthWordsSum[wordIdx] ?? 0) - s0;
            if (pos >= visibleSeqLength || left >= x + w) return x + w + 20;
            const short = this.props.monomerToShort(subParts.getCanonical(pos), maxLengthOfMonomer);
            const labelW = (this._monomerLengthMap[short] ??= g.measureText(short)).width;
            return r === 1 ? left + 1 : r === 2 ? left + labelW - 1 : left + labelW / 2;
          };
          drawConnections(g, conn.conns, xOf, aboveY0, belowY0, conn.geo,
            this.connectionColor(subParts, monomerTextColor, columnTemps.colorCode), hoverPos);
        }
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(errMsg, undefined, errStack);
      this.errors.push(err);
    } finally {
      g.restore();
    }
  }
  /** One line per chain, for a multi-chain row whose connections allow it in a tall cell: a later chain moves under its
   * first partner in an earlier chain, a one-monomer chain (a hub) under the middle of its partners; the chains' own
   * bonds are bracketed above the first line and below the last, inter-chain bonds run straight between the lines.
   * Hit tests use the drawn monomers' bounds, as in multiline mode. */
  private renderStackedChains(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, row: IStackedRow,
  ): void {
    const {seqSS, layout, maxLengthOfMonomer, separator, monomerColor, referenceSequence, selectedPosition,
      highlightEntry, annotRenderer, rowIdx, positionShift} = row;
    const fontSize = this.fontSizeOf();
    const met = measureConnectionsFont(this.props.font, fontSize);
    const chains = layout.chains;
    const n = chains.length;
    const capH = met.capBotOff - met.capTopOff;
    const aboveH = sideHeight(layout.above, layout.geo);
    const belowH = sideHeight(layout.below, layout.geo);
    // the room a taller row leaves goes between the lines, up to most of an em: longer, clearer inter-chain lines
    const fixedH = aboveH + n * capH + Math.max(belowH + 2, met.descBelow);
    let gap = stackedLineGap(fontSize);
    const room = h - fixedH - (n - 1) * gap - 4;
    if (room > 0) gap += Math.min(room / (n - 1), fontSize * 0.7);
    const top = y + (h - (fixedH + (n - 1) * gap)) / 2 + aboveH;
    const capTops = chains.map((_, k) => top + k * (capH + gap));

    const geo = this.stackedGeometry(g, row);
    const x0 = x + this.padding;
    const xOf = (p: number, r?: number): number => x0 + geo.offsets[geo.chainOf[p]] + geo.localX[p] +
      (r === 1 ? 1 : r === 2 ? geo.labelW[p] - 1 : geo.labelW[p] / 2);

    const columnTemps = getPrintColumnTemps(gridCell);
    const {pos: hoverPos, partners} = this.hoverPartners(rowIdx, seqSS, gridCell);
    const cellBounds: IMonomerLayoutData[] = [];
    const boxH = capH + 6;
    for (let k = 0; k < n; k++) {
      const s = chains[k];
      const e = k + 1 < n ? chains[k + 1] : seqSS.length;
      const lineX = x0 + geo.offsets[k];
      const drawY = capTops[k] - met.capTopOff;
      const boxY = capTops[k] - 3;
      // the reference and the selected position count from this chain's start
      const chainReference = referenceSequence.length > 0 ? referenceSequence.slice(s, e) : referenceSequence;
      for (let p = s; p < e; p++) {
        const left = lineX + geo.localX[p];
        if (left >= x + w) break; // the rest of the line is past the cell's edge
        const cm = seqSS.getCanonical(p);
        const labelW = geo.labelW[p];
        if (annotRenderer && rowIdx !== null)
          annotRenderer.drawPositionBackground(g, left, boxY, labelW + geo.sepW, boxH, p, rowIdx, boxY + boxH);
        this.drawHighlightBackground(g, highlightEntry, p, left, boxY, labelW, boxH);
        if (hoverPos !== null && (partners.has(p) || p === hoverPos))
          this.drawPartnerBox(g, left - 1, boxY, labelW + 2, boxH);
        printLeftOrCentered(g, cm, lineX, drawY - 1, w, met.emH, {
          color: monomerColor(cm), pivot: 0, left: true, transparencyRate: 0.0, separator, last: p === e - 1,
          drawStyle: DrawStyle.classic, maxWord: geo.maxWords[k], wordIdx: p - s, gridCell,
          referenceSequence: chainReference, maxLengthOfMonomer, monomerTextSizeMap: this._monomerLengthMap,
          logger: this.logger,
          selectedPosition: isNaN(selectedPosition) || selectedPosition < 1 ? undefined : selectedPosition - s,
          columnTemps,
        });
        cellBounds.push({lineIdx: k, monomerIdx: p - positionShift, sequencePosition: p,
          bounds: new DG.Rect(left - x, capTops[k] - gap / 2 - y, labelW + geo.sepW, capH + gap)});
      }
    }
    if (rowIdx !== null) { // a grid row
      let byRow = this._stackedBounds.get(h);
      if (!byRow || byRow.size >= 20000) {
        if (!byRow && this._stackedBounds.size >= 4) this._stackedBounds.clear();
        this._stackedBounds.set(h, byRow = new Map());
      }
      byRow.set(rowIdx, cellBounds);
    }

    const colorOf = this.connectionColor(seqSS, monomerColor, columnTemps.colorCode);
    drawConnections(g, layout.conns.filter((c) => c.chainFrom === c.chainTo), xOf, capTops[0] - 1.5,
      capTops[n - 1] + capH + 2, layout.geo, colorOf, hoverPos);
    g.save();
    try {
      g.lineCap = 'round';
      for (const c of layout.conns) {
        if (c.chainFrom === c.chainTo) continue;
        const hovered = hoverPos === c.from || hoverPos === c.to;
        g.strokeStyle = colorOf(c);
        g.globalAlpha = hoverPos !== null && !hovered ? 0.3 : 1;
        g.lineWidth = hovered ? layout.geo.lineWidth * 2 : layout.geo.lineWidth;
        g.setLineDash(c.hbond ? DASHED : SOLID);
        g.beginPath();
        g.moveTo(xOf(c.from), capTops[geo.chainOf[c.from]] + capH + 2);
        g.lineTo(xOf(c.to), capTops[geo.chainOf[c.to]] - 1.5);
        g.stroke();
      }
    } finally {
      g.restore();
    }
  }

  /** A stacked row's horizontal geometry, once per row layout (layouts are dropped on reset, with the font and the
   * max monomer length they were laid out with): each chain laid out like a single line from its own start, a later
   * chain shifted under its first partner, a one-monomer chain under the middle of its partners */
  private stackedGeometry(g: CanvasRenderingContext2D, row: IStackedRow): IStackedGeometry {
    let geo = this._stackedGeometries.get(row.layout);
    if (geo) return geo;
    const {sh, seqSS, layout, maxLengthOfMonomer} = row;
    const chains = layout.chains;
    const n = chains.length;
    const len = seqSS.length;
    const sepW = sh.isSeparator() ? this.separatorWidth : this.props.separatorWidth;
    const cw = this.props.fontCharWidth;
    const chainEnd = (k: number): number => k + 1 < n ? chains[k + 1] : len;
    const localX = new Array<number>(len);
    const chainOf = new Array<number>(len);
    const labelW = new Array<number>(len);
    const maxWords: number[][] = [];
    for (let k = 0; k < n; k++) {
      let acc = 0;
      const maxWord: number[] = [];
      for (let p = chains[k]; p < chainEnd(k); p++) {
        const label = this.props.monomerToShort(seqSS.getCanonical(p), maxLengthOfMonomer);
        labelW[p] = (this._monomerLengthMap[label] ??= g.measureText(label)).width;
        localX[p] = acc;
        chainOf[p] = k;
        maxWord.push(acc);
        acc += sepW + label.length * cw;
      }
      maxWord.push(acc);
      maxWords.push(maxWord);
    }
    const offsets = new Array<number>(n).fill(0);
    for (let k = 1; k < n; k++) {
      const toEarlier = layout.conns.filter((c) => c.chainTo === k && c.chainFrom < k);
      if (toEarlier.length === 0) continue;
      const used = chainEnd(k) - chains[k] === 1 ? toEarlier : toEarlier.slice(0, 1);
      let sum = 0;
      for (const c of used)
        sum += offsets[c.chainFrom] + localX[c.from] + labelW[c.from] / 2 - localX[c.to] - labelW[c.to] / 2;
      offsets[k] = sum / used.length;
    }
    const minOffset = Math.min(...offsets);
    for (let k = 0; k < n; k++) offsets[k] -= minOffset;
    geo = {sepW, localX, chainOf, labelW, maxWords, offsets};
    this._stackedGeometries.set(layout, geo);
    return geo;
  }

  /** Sets the connected monomer under the mouse, repainting the grids it leaves and enters (a BILN renderer serves
   * every grid of its column) */
  private setHoverConn(hover: { row: number, pos: number, grid: DG.Grid | null } | null): void {
    const cur = this._hoverConn;
    if (cur?.row === hover?.row && cur?.pos === hover?.pos && cur?.grid?.dart === hover?.grid?.dart) return;
    this._hoverConn = hover;
    this.watchHoverLeave(hover?.grid ?? null);
    for (const grid of [cur?.grid, hover?.grid]) {
      if (grid?.dart) grid.invalidate();
      else if (grid === null) this.invalidateGrid();
    }
  }

  /** Clears the hover when the mouse leaves the hovered cell for a spot with no cell (the empty area below the rows,
   * outside the grid), which the grid does not report to renderers. Listens on the grid's overlay while hovering. */
  private watchHoverLeave(grid: DG.Grid | null): void {
    if (this._hoverLeave?.grid.dart === grid?.dart) return;
    for (const sub of this._hoverLeave?.subs ?? []) sub.unsubscribe();
    this._hoverLeave = null;
    const overlay = grid?.dart ? grid.overlay : null;
    if (!grid || !overlay) return;
    this._hoverLeave = {grid, subs: [
      fromEvent(overlay, 'mouseleave').subscribe(() => this.setHoverConn(null)),
      fromEvent<MouseEvent>(overlay, 'mousemove').subscribe((e) => {
        const hover = this._hoverConn;
        const cell = hover ? grid.hitTest(e.offsetX, e.offsetY) : null;
        const onHoveredCell = cell?.isTableCell && cell.tableRowIndex === hover?.row &&
          cell.tableColumn?.name === this.tableCol.name;
        if (hover && !onHoveredCell)
          this.setHoverConn(null);
      }),
    ]};
  }

  /** The hovered monomer's bonds for its tooltip, e.g. `Bond: C8 ↔ C20 (R3–R3)`, chains named when there are several */
  private connectionsTooltip(seqSS: ISeqSplitted, pos: number, conns: IDrawnConnection[]): HTMLElement {
    const starts = seqSS.graphInfo?.disjointSeqStarts?.length ? seqSS.graphInfo.disjointSeqStarts : [0];
    const where = (p: number): string => {
      let chain = starts.length - 1;
      while (chain > 0 && starts[chain] > p) chain--;
      return `${seqSS.getCanonical(p)}${p - starts[chain] + 1}` + (starts.length > 1 ? ` (chain ${chain + 1})` : '');
    };
    const lines = conns.map((c) => {
      const own = c.from === pos;
      const other = own ? c.to : c.from;
      const rGroups = c.hbond ? '' : ` (R${own ? c.rFrom : c.rTo}–R${own ? c.rTo : c.rFrom})`;
      const kind = c.hbond ? 'Pair' : c.headToTail ? 'Head-to-tail' : 'Bond';
      return ui.divText(`${kind}: ${where(pos)} ↔ ${where(other)}${rGroups}`, {style: {fontSize: '12px'}});
    });
    return ui.divV(lines, {style: {borderTop: '1px solid #ccc', marginTop: '4px', paddingTop: '2px'}});
  }

  private shouldRenderShiftedThreeDots(positionShift: number): boolean {
    return positionShift > 0 && (!this.gridCol || !this.gridCol.dart || !this.gridCol.grid || !this.gridCol.grid.dart || (this.gridCol.grid.props.colHeaderHeight ?? 0) <= 50);
  }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const logPrefix = `${this.toLog()}.onMouseMove()`;
    if (!this.seqHelper || gridCell.tableRowIndex == null) return;

    const positionShift = this.positionShift;
    const gridCellBounds: DG.Rect = gridCell.bounds;

    const argsX = e.offsetX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCellBounds.x);
    const argsY = e.offsetY - gridCellBounds.y;

    // Reset cursor to default, as ellipsis is gone
    if (this.gridCol?.grid?.canvas)
      this.gridCol.grid.canvas.style.cursor = 'default';

    let left: number | null = null;
    const boundsForCell = this._cellBounds.get(gridCell.tableRowIndex!) ??
      this._stackedBounds.get(gridCellBounds.height)?.get(gridCell.tableRowIndex!);

    if (boundsForCell) {
      for (const bound of boundsForCell) {
        if (bound.bounds.contains(argsX, argsY)) {
          left = bound.monomerIdx;
          break;
        }
      }
    } else {
      // Single-line hit detection, by the widths of the labels this cell's height prints
      if (this._showConnections)
        this._widthsH = gridCellBounds.height;
      const leftPadding = this.shouldRenderShiftedThreeDots(positionShift) && (this._leftThreeDotsPadding ?? 0) > 0 ? this._leftThreeDotsPadding : 0;
      left = this.getPosition(gridCell.tableRowIndex!, argsX, gridCellBounds.width, leftPadding);
    }

    this.logger.debug(`${logPrefix}, argsX: ${argsX}, argsY: ${argsY}, left: ${left}`);

    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const seqSS = sh.getSplitted(gridCell.tableRowIndex!);

    // a connected monomer under the mouse marks its partners, and its bonds join its tooltip
    const hoverPos = left !== null && left >= 0 ? left + positionShift : -1;
    const hoverConns = this._showConnections && hoverPos >= 0 ?
      getDrawnConnections(seqSS).filter((c) => c.from === hoverPos || c.to === hoverPos) : [];
    this.setHoverConn(hoverConns.length > 0 ? {row: gridCell.tableRowIndex!, pos: hoverPos, grid: gridCell.grid ?? null} : null);

    if (left !== null && left >= 0 && left + positionShift < seqSS.length) {
      const alphabet = sh.alphabet ?? ALPHABET.UN;
      const seqMonomer = {
        position: left,
        biotype: alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA,
        symbol: seqSS.getCanonical(left + positionShift),
      } as ISeqMonomer;

      const tooltipElements: HTMLElement[] = [];
      let monomerDiv = this._monomerStructureMap[seqMonomer.symbol];
      if (!monomerDiv) {
        const monomerLib = this.getMonomerLib();
        monomerDiv = this._monomerStructureMap[seqMonomer.symbol] = (() => {
          return monomerLib ? monomerLib.getTooltip(seqMonomer.biotype, seqMonomer.symbol) :
            ui.divText('Monomer library is not available');
        })();
      }
      tooltipElements.push(monomerDiv);
      if (hoverConns.length > 0)
        tooltipElements.push(this.connectionsTooltip(seqSS, hoverPos, hoverConns));

      // Append annotation info to tooltip
      const annotR = this.getAnnotationRenderer();
      if (annotR) {
        const annotInfo = annotR.getTooltipInfo(left + positionShift, gridCell.tableRowIndex!);
        if (annotInfo.length > 0) {
          const sep = document.createElement('hr');
          sep.style.margin = '4px 0';
          sep.style.border = 'none';
          sep.style.borderTop = '1px solid #ccc';
          tooltipElements.push(sep);
          for (const line of annotInfo)
            tooltipElements.push(ui.divText(line, {style: {fontSize: '12px', marginBottom: '2px'}}));
        }
      }

      ui.tooltip.show(ui.divV(tooltipElements), e.x + 16, e.y + 16);

      execMonomerHoverLinks(gridCell, seqMonomer);
    } else {
      if (left === -1)
        ui.tooltip.show(ui.divText(`${Math.min(positionShift, seqSS.length)} hidden monomers`), e.x + 16, e.y + 16);
      else
        ui.tooltip.hide();

      execMonomerHoverLinks(gridCell, null);
    }
  }

  override onMouseLeave(_gridCell: DG.GridCell, _e: MouseEvent): void {
    this.setHoverConn(null);
  }
}

export function getUpdatedWidth(
  grid: DG.Grid | null | undefined, g: CanvasRenderingContext2D, x: number, w: number, dpr: number
): number {
  return grid?.dart ? Math.max(Math.min(grid.canvas.width / dpr - x, w)) : Math.max(g.canvas.width / dpr - x, 0);
}
