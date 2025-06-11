/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {getSeqHelper, ISeqHelper} from './seq-helper';
import {ALPHABET, MonomerToShortFunc, NOTATION, SplitterFunc, TAGS as bioTAGS} from './macromolecule';
import {ISeqSplitted} from './macromolecule/types';
import {CellRendererBackBase, getGridCellColTemp} from './cell-renderer-back-base';
import {ILogger} from './logger';
import {getMonomerLibHelper} from '../monomer-works/monomer-utils';
import {errInfo} from './err-info';
import {DrawStyle, printLeftOrCentered, TAGS as mmcrTAGS, PrintOptions} from './cell-renderer';
import {MmcrTemps, rendererSettingsChangedState, tempTAGS} from './cell-renderer-consts';
import {IMonomerLibBase} from '../types/index';
import {HelmTypes} from '../helm/consts';
import {ISeqMonomer} from '../helm/types';
import {execMonomerHoverLinks} from '../monomer-works/monomer-hover';
import * as operators from 'rxjs/operators';

type MonomerPlacerProps = {
  separatorWidth: number,
  monomerToShort: MonomerToShortFunc,
  font: string, fontCharWidth: number,
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

  private seqHelper!: ISeqHelper;

  private sysMonomerLib: IMonomerLibBase | null = null;

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

    if (this.gridCol) {
      this.subs.push(this.gridCol.grid.onAfterDrawContent.subscribe(() => {
        this._onRendered.next();
      }));
    }

    if (this.tableCol && this.gridCol) {
      this.subs.push(this.tableCol.dataFrame.onCurrentRowChanged.subscribe(() => {
        const df = this.tableCol.dataFrame;
        const grid = this.gridCol!.grid;
        if (df.currentRowIdx === -1) {
          this.tableCol.temp[tempTAGS.referenceSequence] = null;
          this.tableCol.temp[tempTAGS.currentWord] = null;
          this.invalidateGrid();
        }
      }));

      this.subs.push(DG.debounce(this.tableCol.dataFrame.onMetadataChanged.pipe(
        operators.filter((a) => a.args.source === this.tableCol && a.args.key === bioTAGS.positionShift)
      ), 200).subscribe((_) => {
        this.reset();
      }));
    }
  }

  private _multiLineMonomerBounds: Array<{
  lineIdx: number;
  monomerIdx: number;
  bounds: DG.Rect;
  sequencePosition: number;
}> = [];


  private calculateFontBasedSpacing(g: CanvasRenderingContext2D): {
  lineHeight: number;
  monomerSpacing: number;
} {
  // Get font metrics to derive spacing
    const metrics = g.measureText('M'); // Use 'M' as it's typically the widest character
    const fontSize = parseInt(this.props.font.match(/(\d+)px/)?.[1] || '12');

    // Line height should be font size + some padding (typically 1.2-1.5x font size)
    const lineHeight = Math.max(fontSize * 1.4, metrics.fontBoundingBoxAscent + metrics.fontBoundingBoxDescent + 4);

    // Monomer spacing should be proportional to character width
    const monomerSpacing = Math.max(2, this.props.fontCharWidth * 0.2);

    return {lineHeight, monomerSpacing};
  }


  private shouldUseMultilineRendering(tableCol: DG.Column): boolean {
    const renderMultiline = tableCol.getTag('renderMultiline');
    return renderMultiline === 'true';
  }


  // This stuff doesn't seem to trigger for whatever reason on my ellipsis element. Leaving it out for now.
  override onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!this._ellipsisBounds || !this.gridCol?.grid) return;

    const grid = this.gridCol.grid;
    const gridCellBounds: DG.Rect = gridCell.bounds;
    // Calculate mouse coordinates relative to the cell's top-left corner
    const argsX = e.offsetX - gridCellBounds.x;
    const argsY = e.offsetY - gridCellBounds.y;

    if (this._ellipsisBounds.contains(argsX, argsY)) {
      (grid.props as any).allowDynamicRowHeight = true;

      const requiredHeight = (this._totalLinesNeeded * this._lineHeight) + (this.padding * 2);

      (grid as any).setRowHeight(gridCell.gridRow, requiredHeight);

      grid.invalidate();
    }
  }


  private calculateMultiLineLayoutDynamic(
    g: CanvasRenderingContext2D,
    x: number,
    y: number,
    w: number,
    h: number,
    subParts: ISeqSplitted,
    positionShift: number,
    maxLengthOfMonomer: number
  ): {
    lineLayouts: Array<{
      lineIdx: number;
      elements: Array<{
        posIdx?: number;
        x: number;
        width: number;
        om: string;
        isSeparator: boolean;
      }>;
    }>;
    lineHeight: number;
  } {
    // --- 1. Setup ---
    const {lineHeight, monomerSpacing} = this.calculateFontBasedSpacing(g);
    const availableWidth = w - (this.padding * 2);
    const separator = this.tableCol.getTag(bioTAGS.separator) ?? '';
    const separatorWidth = separator ? g.measureText(separator).width : 0;

    // --- 2. Create a unified list of elements (monomers + separators) ---
    const elementsToLayout: { text: string, width: number, isSeparator: boolean, posIdx?: number }[] = [];
    for (let i = positionShift; i < subParts.length; i++) {
      const om = subParts.getOriginal(i);
      const shortMon = this.props.monomerToShort(om, maxLengthOfMonomer);
      const monomerWidth = g.measureText(shortMon).width;
      elementsToLayout.push({text: shortMon, width: monomerWidth, isSeparator: false, posIdx: i});

      if (separator && i < subParts.length - 1)
        elementsToLayout.push({text: separator, width: separatorWidth, isSeparator: true, posIdx: undefined});
    }

    // --- 3. Virtual layout to create all possible lines ---
    const allLineElements: Array<Array<typeof elementsToLayout[0]>> = [];
    if (elementsToLayout.length > 0) {
      let currentElementIdx = 0;
      while (currentElementIdx < elementsToLayout.length) {
        const lineElements: Array<typeof elementsToLayout[0]> = [];
        let currentLineWidth = 0;
        while (currentElementIdx < elementsToLayout.length) {
          const element = elementsToLayout[currentElementIdx];
          const elementWidthWithSpacing = element.width + (lineElements.length > 0 ? monomerSpacing : 0);
          if (currentLineWidth + elementWidthWithSpacing > availableWidth && lineElements.length > 0)
            break;

          lineElements.push(element);
          currentLineWidth += elementWidthWithSpacing;
          currentElementIdx++;
        }
        allLineElements.push(lineElements);
      }
    }

    // --- 4. Determine how many lines to actually render ---
    const availableHeightForLines = h - (this.padding * 2);
    const linesToRenderCount = Math.max(0, Math.floor(availableHeightForLines / lineHeight));

    // --- 5. Generate final layout for the visible lines ---
    const lineLayouts: Array<{ lineIdx: number, elements: any[] }> = [];
    for (let i = 0; i < linesToRenderCount && i < allLineElements.length; i++) {
      const lineData = allLineElements[i];
      if (!lineData) continue;

      const elementsForLine = [];
      let currentX = x + this.padding;
      for (const element of lineData) {
        elementsForLine.push({
          posIdx: element.posIdx,
          x: currentX,
          width: element.width,
          om: element.text,
          isSeparator: element.isSeparator,
        });
        currentX += element.width + monomerSpacing;
      }
      lineLayouts.push({lineIdx: i, elements: elementsForLine});
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

    this.reset();
  }

  public static getFontSettings(tableCol?: DG.Column) {
    const fontFamily = 'monospace';
    let fontSize = 12;
    if (tableCol && tableCol.temp[MmcrTemps.fontSize] &&
        typeof tableCol.temp[MmcrTemps.fontSize] === 'number' && !isNaN(tableCol.temp[MmcrTemps.fontSize]))
      fontSize = Math.max(tableCol.temp[MmcrTemps.fontSize], 1);
    // monospace ratio
    const fontWidth = fontSize * 0.6;
    return {font: `${fontSize}px ${fontFamily}`, fontWidth};
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

  private getCellMonomerLengthsForSeq(rowIdx: number): number[] {
    const logPrefix = `${this.toLog()}.getCellMonomerLengthsForSeq()`;
    // this.logger.debug(`${logPrefix}, start`);

    if (this._monomerLengthList === null)
      this._monomerLengthList = new Array(this.tableCol.length).fill(null);

    const visibleSeqStart = this.positionShift;
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
    const maxVisibleSeqLength: number = Math.ceil(this.colWidth / minMonWidth) + visibleSeqStart;
    const seqSS: ISeqSplitted = sh.getSplitted(rowIdx);
    const visibleSeqEnd: number = Math.min(maxVisibleSeqLength, seqSS.length);

    let res: number[] = this._monomerLengthList[rowIdx];
    if (res === null || res.length != visibleSeqEnd - visibleSeqStart) {
      res = this._monomerLengthList[rowIdx] = new Array<number>(seqSS.length);

      let seqWidth: number = 0;
      for (let seqMonI = visibleSeqStart; seqMonI < visibleSeqEnd; ++seqMonI) {
        const seqMonLabel = seqSS.getOriginal(seqMonI);
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

  private getCellMonomerLengthsForSeqValue(value: string, width: number): number[] {
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
    const visibleSeqStart = this.positionShift;
    const maxVisibleSeqLength: number = Math.ceil(width / minMonWidth) + visibleSeqStart;
    const seqSS: ISeqSplitted = sh.splitter(value);
    const visibleSeqEnd: number = Math.min(maxVisibleSeqLength, seqSS.length);
    const res = new Array<number>(visibleSeqEnd - visibleSeqStart);
    let seqWidth: number = 0;
    for (let seqMonI = visibleSeqStart; seqMonI < visibleSeqEnd; ++seqMonI) {
      const seqMonLabel = seqSS.getOriginal(seqMonI);
      const shortMon: string = this.props.monomerToShort(seqMonLabel, this.monomerLengthLimit);
      const separatorWidth = sh.isSeparator() ? this.separatorWidth : this.props.separatorWidth;
      const seqMonWidth: number = separatorWidth + shortMon.length * this.props.fontCharWidth;
      res[seqMonI - visibleSeqStart] = seqMonWidth;
      seqWidth += seqMonWidth;
      if (seqWidth > width) break;
    }
    return res;
  }

  private getCellMonomerLengthsForSeqMsa(): number[] {
    const logPrefix = `${this.toLog()}.getCellMonomerLengthsForSeqMsa()`;
    // this.logger.debug(`${logPrefix}, start`);
    if (this._monomerLengthList === null)
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

      const visibleSeqEnd: number = Math.min(maxVisibleSeqLength, seqSS.length);
      if (visibleSeqEnd - visibleSequenceStart > res.length)
        res.push(...new Array<number>(visibleSeqEnd - visibleSequenceStart - res.length).fill(minMonWidth));

      let seqWidth: number = 0;
      for (let seqMonI = visibleSequenceStart; seqMonI < visibleSeqEnd; ++seqMonI) {
        const seqMonLabel: string = seqSS.getOriginal(seqMonI);
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
  // Check if we have multi-line bounds stored
    if (this._multiLineMonomerBounds.length > 0) {
    // Use multi-line hit detection
      return this.getPositionMultiLine(x, 0); // y is handled by bounds checking
    }

    // Fall back to single-line detection
    const [_monomerMaxLengthList, monomerMaxLengthSumList]: [number[], number[]] =
    this.getCellMonomerLengths(rowIdx, width);
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const seqSS = sh.getSplitted(rowIdx);
    if (seqSS.length === 0) return null;
    return hitBounds(monomerMaxLengthSumList, x, positionShiftPadding);
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

  private _leftThreeDotsPadding: number = 0;

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle
  ) {
    // for cases when we render it on somewhere other than grid, gridRow might be null or incorrect (set to 0).
    //for this case we can just recalculate split sequence without caching
    const isRenderedOnGrid = gridCell.grid?.canvas === g.canvas;

    if (!this.seqHelper) return;
    const gridCol = this.gridCol;
    const tableCol = this.tableCol;
    const dpr = window.devicePixelRatio;
    const positionShift = this.positionShift;


    const logPrefix = `${this.toLog()}.render()`;
    this.logger.debug(`${logPrefix}, start`);

    // Cell renderer settings
    let maxLengthOfMonomer: number = this.monomerLengthLimit;
    if (mmcrTAGS.maxMonomerLength in tableCol.tags) {
      const v = parseInt(tableCol.getTag(mmcrTAGS.maxMonomerLength));
      maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
    }
    if (MmcrTemps.maxMonomerLength in tableCol.temp) {
      const v = tableCol.temp[MmcrTemps.maxMonomerLength];
      maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
    }
    //TODO: this is hacky, but i wanted to leave a parameter for max visible lines initally. Perhaps this should just be calculated dynamically.
    const maxVisibleLines = parseInt(tableCol.getTag('maxVisibleLines') ?? '999');


    g.save();
    try {
      const sh = this.seqHelper.getSeqHandler(tableCol);

      if (
        tableCol.temp[MmcrTemps.rendererSettingsChanged] === rendererSettingsChangedState.true ||
        this.monomerLengthLimit != maxLengthOfMonomer
      ) {
        let gapLength = 0;
        const msaGapLength = 8;
        gapLength = tableCol.temp[MmcrTemps.gapLength] as number ?? gapLength;
        // this event means that the mm renderer settings have changed,
        // particularly monomer representation and max width.
        this.setMonomerLengthLimit(maxLengthOfMonomer);
        this.setSeparatorWidth(sh.isMsa() ? msaGapLength : gapLength);
        tableCol.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.false;
      }

      let [maxLengthWords, maxLengthWordsSum]: [number[], number[]] =
        this.getCellMonomerLengths(gridCell.tableRowIndex!, w);
      const _maxIndex = maxLengthWords.length;

      const value: any = gridCell.cell.value;
      const rowIdx = gridCell.cell.rowIndex;
      const minDistanceRenderer = 50;
      if (isRenderedOnGrid)
        w = getUpdatedWidth(gridCol?.grid, g, x, w, dpr);

      g.beginPath();
      g.rect(x + this.padding, y, w - this.padding - 1, h);
      g.clip();
      g.font = this.props?.font ?? '12px monospace';
      g.textBaseline = 'top';


      //TODO: can this be replaced/merged with splitSequence?
      const units = tableCol.meta.units;
      const aligned: string = tableCol.getTag(bioTAGS.aligned);

      const separator = tableCol.getTag(bioTAGS.separator) ?? '';
      const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
      const splitLimit = Math.ceil(w / minMonWidth) + positionShift;

      const tempReferenceSequence: string | null = tableCol.temp[tempTAGS.referenceSequence];
      const tempCurrentWord: string | null = this.tableCol.temp[tempTAGS.currentWord];
      if (tempCurrentWord && tableCol?.dataFrame?.currentRowIdx === -1)
        this.tableCol.temp[tempTAGS.currentWord] = null;

      const referenceSequence: string[] = (() => {
        // @ts-ignore
        const splitterFunc: SplitterFunc = sh.getSplitter(splitLimit);
        const seqSS = splitterFunc(
          ((tempReferenceSequence != null) && (tempReferenceSequence != '')) ?
            tempReferenceSequence : tempCurrentWord ?? '');
        return wu.count(0).take(seqSS.length).slice(positionShift).map((posIdx) => seqSS.getOriginal(posIdx)).toArray();
      })();

      const subParts: ISeqSplitted = isRenderedOnGrid ? sh.getSplitted(rowIdx) : sh.splitter(value);

      if (!isRenderedOnGrid)
        maxLengthWordsSum = this.getSummedMonomerLengths(this.getCellMonomerLengthsForSeqValue(value, w));

      let drawStyle = DrawStyle.classic;
      if (aligned && aligned.includes('MSA') && units == NOTATION.SEPARATOR)
        drawStyle = DrawStyle.MSA;

      const shouldUseMultiLine = this.shouldUseMultilineRendering(tableCol) && drawStyle !== DrawStyle.MSA;

      // if the sequence is rendered in shifted mode, we will also render three dots at start, indicating the shift
      this._leftThreeDotsPadding = this.shouldRenderShiftedThreeDots(positionShift) ? g.measureText(shiftedLeftPaddingText).width : 0;
      const selectedPosition = Number.parseInt(tableCol.getTag(bioTAGS.selectedPosition) ?? '-200');
      const visibleSeqLength = Math.min(subParts.length, splitLimit);


      if (shouldUseMultiLine) {
        // Clear previous bounds for hit detection
        this._multiLineMonomerBounds = [];
        this._ellipsisBounds = undefined; // Clear any leftover state

        // Call the updated layout function
        const layout = this.calculateMultiLineLayoutDynamic(
          g, x, y, w, h, subParts, positionShift, maxLengthOfMonomer
        );


        for (const lineLayout of layout.lineLayouts) {
          const lineY = y + (lineLayout.lineIdx * layout.lineHeight);

          for (const element of lineLayout.elements) {
          // Handle separators differently from monomers
            if (element.isSeparator) {
              printLeftOrCentered(g, element.om, element.x, lineY, element.width + 4, layout.lineHeight, {
                color: undefinedColor,
                isMultiLineContext: true,
              });
              continue;
            }

            // --- Existing monomer rendering logic ---
            const monomer = element;
            const cm: string = monomer.posIdx as number < subParts.length ? subParts.getCanonical(monomer.posIdx as number) : sh.defaultGapOriginal;
            let color = undefinedColor;
            const monomerLib = this.getMonomerLib();
            if (monomerLib) {
              const biotype = sh.defaultBiotype;
              color = monomerLib.getMonomerTextColor(biotype, cm);
            }

            const last = monomer.posIdx === subParts.length - 1;
            let transparencyRate = 0.0;
            // ... (rest of transparency logic is unchanged)
            if (referenceSequence && referenceSequence.length > 0) {
              const refIndex = monomer.posIdx as number - positionShift;
              if (refIndex >= 0 && refIndex < referenceSequence.length) {
                const currentMonomerCanonical = referenceSequence[refIndex];
                const compareWithCurrent = gridCell?.cell?.column?.temp?.['compare-with-current'] ?? true;
                const highlightDifference = gridCell?.cell?.column?.temp?.['highlight-difference'] ?? 'difference';

                if (compareWithCurrent && highlightDifference === 'difference')
                  transparencyRate = (monomer.om === currentMonomerCanonical) ? 0.7 : 0.0;
                else if (compareWithCurrent && highlightDifference === 'equal')
                  transparencyRate = (monomer.om !== currentMonomerCanonical) ? 0.7 : 0.0;
              }
            }

            const opts: Partial<PrintOptions> = {
              color: color,
              pivot: 0,
              left: true,
              transparencyRate: transparencyRate,
              separator: '', // No separator in multiline mode
              last: last,
              drawStyle: drawStyle,
              maxWord: maxLengthWordsSum,
              wordIdx: monomer.posIdx as number - positionShift,
              gridCell: gridCell,
              referenceSequence: referenceSequence,
              maxLengthOfMonomer: maxLengthOfMonomer,
              monomerTextSizeMap: this._monomerLengthMap,
              logger: this.logger,
              selectedPosition: isNaN(selectedPosition) || selectedPosition < 1 ? undefined : selectedPosition - positionShift,
              isMultiLineContext: true,
              lineNumber: lineLayout.lineIdx,
            };

            // Store bounds for hit detection
            this._multiLineMonomerBounds.push({
              lineIdx: lineLayout.lineIdx,
              monomerIdx: monomer.posIdx as number - positionShift,
              bounds: new DG.Rect(
                monomer.x - x, // Make relative to cell
                lineY - y, // Make relative to cell
                monomer.width,
                layout.lineHeight
              ),
              sequencePosition: monomer.posIdx as number
            });

            printLeftOrCentered(g, monomer.om, monomer.x, lineY, monomer.width + 4, layout.lineHeight, opts);
          }
        }
      } else {
        this._multiLineMonomerBounds = [];
        this._ellipsisBounds = undefined;


        // Existing single-line rendering (unchanged)
        for (let posIdx: number = positionShift; posIdx < visibleSeqLength; ++posIdx) {
          const om: string = posIdx < subParts.length ? subParts.getOriginal(posIdx) : sh.defaultGapOriginal;
          const cm: string = posIdx < subParts.length ? subParts.getCanonical(posIdx) : sh.defaultGapOriginal;

          let color = undefinedColor;
          const monomerLib = this.getMonomerLib();
          if (monomerLib) {
            const biotype = sh.defaultBiotype;
            color = monomerLib.getMonomerTextColor(biotype, cm);
          }
          g.fillStyle = undefinedColor;
          const last = posIdx === subParts.length - 1;
          const opts: Partial<PrintOptions> = {
            color: color, pivot: 0, left: true, transparencyRate: 0.0,
            separator: separator, last: last,
            drawStyle: drawStyle, maxWord: maxLengthWordsSum, wordIdx: posIdx - positionShift, gridCell: gridCell,
            referenceSequence: referenceSequence, maxLengthOfMonomer: maxLengthOfMonomer,
            monomerTextSizeMap: this._monomerLengthMap, logger: this.logger,
            selectedPosition: isNaN(selectedPosition) || selectedPosition < 1 ? undefined : selectedPosition - positionShift,
          };
          printLeftOrCentered(g, om, x + this.padding + this._leftThreeDotsPadding, y, w, h, opts);
          if (minDistanceRenderer > w) break;
        }
      }


      // Handle position shift dots (only for single-line)
      if (this.shouldRenderShiftedThreeDots(positionShift) && !shouldUseMultiLine) {
        const opts: Partial<PrintOptions> = {
          color: undefinedColor, pivot: 0, left: true, transparencyRate: 1.0, separator: separator, last: false,
          drawStyle: drawStyle, maxWord: maxLengthWordsSum, wordIdx: 0, gridCell: gridCell,
          maxLengthOfMonomer: maxLengthOfMonomer,
          monomerTextSizeMap: this._monomerLengthMap, logger: this.logger,
        };
        printLeftOrCentered(g, shiftedLeftPaddingText, x + this.padding, y, w, h, opts);
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(errMsg, undefined, errStack);
      this.errors.push(err);
      //throw err; // Do not throw to prevent disabling renderer
    } finally {
      g.restore();
    }
  }
  private getPositionMultiLine(x: number, y: number): number | null {
    for (const bound of this._multiLineMonomerBounds) {
    // Check if the point is within the monomer bounds
      if (x >= bound.bounds.x && x <= bound.bounds.x + bound.bounds.width &&
        y >= bound.bounds.y && y <= bound.bounds.y + bound.bounds.height)
        return bound.monomerIdx;
    }
    return null;
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

    // Check if we're in multi-line mode for monomer detection
    if (this._multiLineMonomerBounds.length > 0) {
      for (const bound of this._multiLineMonomerBounds) {
        if (argsX >= bound.bounds.x && argsX <= bound.bounds.x + bound.bounds.width &&
          argsY >= bound.bounds.y && argsY <= bound.bounds.y + bound.bounds.height) {
          left = bound.monomerIdx;
          break;
        }
      }
    } else {
      // Single-line hit detection
      const leftPadding = this.shouldRenderShiftedThreeDots(positionShift) && (this._leftThreeDotsPadding ?? 0) > 0 ? this._leftThreeDotsPadding : 0;
      left = this.getPosition(gridCell.tableRowIndex!, argsX, gridCellBounds.width, leftPadding);
    }

    this.logger.debug(`${logPrefix}, argsX: ${argsX}, argsY: ${argsY}, left: ${left}`);

    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const seqSS = sh.getSplitted(gridCell.tableRowIndex!);

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
}

export function getUpdatedWidth(
  grid: DG.Grid | null | undefined, g: CanvasRenderingContext2D, x: number, w: number, dpr: number
): number {
  return !!grid ? Math.max(Math.min(grid.canvas.width / dpr - x, w)) : Math.max(g.canvas.width / dpr - x, 0);
}
