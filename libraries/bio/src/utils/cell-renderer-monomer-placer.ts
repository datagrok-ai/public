/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {getSeqHelper, ISeqHelper} from './seq-helper';
import {ALPHABET, MonomerToShortFunc, NOTATION, SplitterFunc, TAGS as bioTAGS} from './macromolecule';
import {ISeqSplitted} from './macromolecule/types';
import {CellRendererBackBase} from './cell-renderer-back-base';
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
  font: string,
  fontCharWidth: number,
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

interface IMultilineLayoutElement {
  posIdx?: number;
  x: number;
  width: number;
  om: string; // original monomer or separator text
  isSeparator: boolean;
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

  private _cellBounds: Map<number, IMonomerLayoutData[]> = new Map();

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
        if (df.currentRowIdx === -1) {
          this.tableCol.temp[tempTAGS.referenceSequence] = null;
          this.tableCol.temp[tempTAGS.currentWord] = null;
          this.invalidateGrid();
        }
      }));
      const resetTriggerTags = [
        bioTAGS.positionShift,
        'renderMultiline'
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
    let fontSize = 12; // Default font size
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
    const {lineHeight, monomerSpacing} = this.calculateFontBasedSpacing(g);
    const availableWidth = w - (this.padding * 2);

    // --- 2. Find the widest monomer in the sequence to set a uniform column width ---
    let maxMonomerWidth = 0;
    const monomers: {text: string, posIdx: number}[] = [];
    if (subParts.length > 0) {
      for (let i = positionShift; i < subParts.length; i++) {
        const om = subParts.getOriginal(i);
        const shortMon = this.props.monomerToShort(om, maxLengthOfMonomer);
        monomers.push({text: shortMon, posIdx: i});
        maxMonomerWidth = Math.max(maxMonomerWidth, g.measureText(shortMon).width);
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
          posIdx: monomer.posIdx,
          x: xPos,
          width: uniformColumnWidth, // Use the uniform width for all slots
          om: monomer.text,
          isSeparator: false,
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
    this._cellBounds.clear();
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
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle,
  ) {
    // for cases when we render it on somewhere other than grid, gridRow might be null or incorrect (set to 0).
    //for this case we can just recalculate split sequence without caching
    const isRenderedOnGrid = gridCell.grid?.canvas === g.canvas;

    if (!this.seqHelper) return;
    const tableCol = this.tableCol;
    const positionShift = this.positionShift;

    g.save();
    try {
      const sh = this.seqHelper.getSeqHandler(tableCol);
      let maxLengthOfMonomer: number = this.monomerLengthLimit;
      if (mmcrTAGS.maxMonomerLength in tableCol.tags) {
        const v = parseInt(tableCol.getTag(mmcrTAGS.maxMonomerLength));
        maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
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


      const rowIdx = gridCell.cell.rowIndex;
      const value: any = gridCell.cell.value;
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


      const shouldUseMultiLine = this.shouldUseMultilineRendering(tableCol);

      if (shouldUseMultiLine) {
        const currentCellBounds: IMonomerLayoutData[] = [];
        const layout = this.calculateMultiLineLayoutDynamic(g, w, h, subParts, positionShift, maxLengthOfMonomer);

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

          for (const element of lineLayout.elements) {
            const elementX = x + element.x;
            const monomer = element;
            const monomerIndex = monomer.posIdx!;
            const cm = subParts.getCanonical(monomerIndex);

            let color = undefinedColor;
            const monomerLib = this.getMonomerLib();
            if (monomerLib)
              color = monomerLib.getMonomerTextColor(sh.defaultBiotype, cm);

            let transparencyRate = 0.0;
            if (gridCell.tableRowIndex !== tableCol.dataFrame.currentRowIdx && referenceSequence.length > 0) {
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

            printLeftOrCentered(g, monomer.om, elementX, lineY, element.width, layout.lineHeight, {
              color: color,
              isMultiLineContext: true,
              transparencyRate: transparencyRate,
              selectedPosition: isNaN(selectedPosition) || selectedPosition < 1 ? undefined : selectedPosition,
              wordIdx: monomerIndex
            });
          }
        }
        if (gridCell.tableRowIndex !== null)
          this._cellBounds.set(gridCell.tableRowIndex, currentCellBounds);
      } else {
        // --- Single-line rendering path (is unchanged) ---
        this._leftThreeDotsPadding = this.shouldRenderShiftedThreeDots(positionShift) ? g.measureText(shiftedLeftPaddingText).width : 0;
        let [, maxLengthWordsSum]: [number[], number[]] = this.getCellMonomerLengths(gridCell.tableRowIndex!, w);
        if (!isRenderedOnGrid)
          maxLengthWordsSum = this.getSummedMonomerLengths(this.getCellMonomerLengthsForSeqValue(value, w));
        const minMonomerWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
        const visibleSeqLength = Math.min(subParts.length, Math.ceil(w / (minMonomerWidth)) + positionShift);

        for (let posIdx: number = positionShift; posIdx < visibleSeqLength; ++posIdx) {
          const om: string = posIdx < subParts.length ? subParts.getOriginal(posIdx) : sh.defaultGapOriginal;
          const cm: string = posIdx < subParts.length ? subParts.getCanonical(posIdx) : sh.defaultGapOriginal;

          let color = undefinedColor;
          if (this.getMonomerLib())
            color = this.getMonomerLib()!.getMonomerTextColor(sh.defaultBiotype, cm);

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
        }
        if (this.shouldRenderShiftedThreeDots(positionShift)) {
          const opts = {
            color: undefinedColor, pivot: 0, left: true, transparencyRate: 0, separator: separator, last: false,
            drawStyle: drawStyle, maxWord: maxLengthWordsSum, wordIdx: 0, gridCell: gridCell,
            maxLengthOfMonomer: maxLengthOfMonomer,
            monomerTextSizeMap: this._monomerLengthMap, logger: this.logger,
          };
          printLeftOrCentered(g, shiftedLeftPaddingText, x + this.padding, y, w, h, opts);
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
    const boundsForCell = this._cellBounds.get(gridCell.tableRowIndex!);

    if (boundsForCell) {
      for (const bound of boundsForCell) {
        if (bound.bounds.contains(argsX, argsY)) {
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
