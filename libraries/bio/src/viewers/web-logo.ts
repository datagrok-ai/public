import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import * as rxjs from 'rxjs';

import {Vector} from '@datagrok-libraries/utils/src/type-declarations';
import {vectorLength, vectorDotProduct} from '@datagrok-libraries/utils/src/vector-operations';
import {Aminoacids, AminoacidsPalettes} from '../aminoacids';
import {Nucleotides, NucleotidesPalettes} from '../nucleotides';
import {UnknownSeqPalettes} from '../unknown';
import {SeqPalette} from '../seq-palettes';
import {Subscription} from 'rxjs';
import {NOTATION, UnitsHandler} from '../utils/units-handler';

declare module 'datagrok-api/src/grid' {
  interface Rect {
    contains(x: number, y: number): boolean;
  }
}

enum PositionHeight {
  Entropy = 'Entropy',
  full = '100%',
}

declare global {
  interface HTMLCanvasElement {
    getCursorPosition(event: MouseEvent, r: number): DG.Point;
  }
}

export type MonomerFreqs = { [m: string]: number };
export type SeqColStats = { freq: MonomerFreqs, sameLength: boolean }
export type SplitterFunc = (seq: string) => string[];

/**@param {MouseEvent} event
 * @param {number} r devicePixelRation
 * @return {DG.Point} canvas related cursor position
 */
HTMLCanvasElement.prototype.getCursorPosition = function(event: MouseEvent, r: number): DG.Point {
  const rect = this.getBoundingClientRect();
  return new DG.Point((event.clientX - rect.left) * r, (event.clientY - rect.top) * r);
};

DG.Rect.prototype.contains = function(x: number, y: number): boolean {
  return this.left <= x && x <= this.right && this.top <= y && y <= this.bottom;
};

export class PositionMonomerInfo {
  /** Sequences count with monomer in position */
  count: number;

  /** Remember screen coords rect */
  bounds: DG.Rect;

  constructor(count: number = 0, bounds: DG.Rect = new DG.Rect(0, 0, 0, 0)) {
    this.count = count;
    this.bounds = bounds;
  }
}

export class PositionInfo {
  public readonly name: string;
  freq: { [m: string]: PositionMonomerInfo };
  rowCount: number;
  sumForHeightCalc: number;

  /** freq = {}, rowCount = 0
   * @param {string} name Name of position ('111A', '111.1', etc)
   * @param {number} sumForHeightCalc Sum of all monomer counts for height calculation
   * @param {number} rowCount Count of elements in column
   * @param {string[]} freq frequency of monomers in position
   */
  constructor(name: string, freq: { [m: string]: PositionMonomerInfo } = {}, rowCount: number = 0, sumForHeightCalc: number = 0) {
    this.name = name;
    this.freq = freq;
    this.rowCount = rowCount;
    this.sumForHeightCalc = sumForHeightCalc;
  }
}


export class WebLogo extends DG.JsViewer {
  public static residuesSet = 'nucleotides';
  private static viewerCount: number = -1;

  private readonly viewerId: number = -1;
  private unitsHandler: UnitsHandler | null;
  private initialized: boolean = false;

  // private readonly colorScheme: ColorScheme = ColorSchemes[NucleotidesWebLogo.residuesSet];
  protected cp: SeqPalette | null = null;

  private host?: HTMLDivElement;
  private msgHost?: HTMLElement;
  private canvas?: HTMLCanvasElement;
  private slider?: DG.RangeSlider;
  private readonly textBaseline: CanvasTextBaseline;

  private axisHeight: number = 12;

  private seqCol: DG.Column<string> | null = null;
  private splitter: SplitterFunc | null = null;
  // private maxLength: number = 100;
  private positions: PositionInfo[] = [];

  private rowsMasked: number = 0;
  private rowsNull: number = 0;

  // Viewer's properties (likely they should be public so that they can be set outside)
  private _positionWidth: number;
  public positionWidth: number;
  public minHeight: number;
  public maxHeight: number;
  public considerNullSequences: boolean;
  public sequenceColumnName: string | null;
  public positionMarginState: string;
  public positionMargin: number = 0;
  public startPositionName: string | null;
  public endPositionName: string | null;
  public fixWidth: boolean;
  public verticalAlignment: string | null;
  public horizontalAlignment: string | null;
  public fitArea: boolean;
  public shrinkEmptyTail: boolean;
  public skipEmptyPositions: boolean;
  public positionHeight: string;

  private positionNames: string[] = [];

  private startPosition: number = -1;

  private endPosition: number = -1;

  /** For startPosition equals to endPosition Length is 1 */
  private get Length(): number {
    if (this.skipEmptyPositions) {
      return this.positions.length;
    }
    return this.startPosition <= this.endPosition ? this.endPosition - this.startPosition + 1 : 0;
  }

  /** Calculate new position data basic on {@link positionMarginState} and {@link positionMargin} */
  private get positionWidthWithMargin() {
    if ((this.positionMarginState === 'auto') && (this.unitsHandler?.getAlphabetIsMultichar() === true)) {
      return this._positionWidth + this.positionMargin;
    }
    if (this.positionMarginState === 'enable') {
      return this._positionWidth + this.positionMargin;
    }

    return this._positionWidth;
  }

  private viewSubs: Subscription[] = [];

  constructor() {
    super();

    this.viewerId = WebLogo.viewerCount;
    WebLogo.viewerCount += 1;

    this.textBaseline = 'top';
    this.unitsHandler = null;

    this._positionWidth = this.positionWidth = this.float('positionWidth', 16/*,
      {editor: 'slider', min: 4, max: 64, postfix: 'px'}*/);
    this.minHeight = this.float('minHeight', 50/*,
      {editor: 'slider', min: 25, max: 250, postfix: 'px'}*/);
    this.maxHeight = this.float('maxHeight', 100/*,
      {editor: 'slider', min: 25, max: 500, postfix: 'px'}*/);

    this.considerNullSequences = this.bool('considerNullSequences', false);
    this.sequenceColumnName = this.string('sequenceColumnName', null);

    this.startPositionName = this.string('startPositionName', null);
    this.endPositionName = this.string('endPositionName', null);

    this.fixWidth = this.bool('fixWidth', false);

    this.verticalAlignment = this.string('verticalAlignment', 'middle',
      {choices: ['top', 'middle', 'bottom']});
    this.horizontalAlignment = this.string('horizontalAlignment', 'center',
      {choices: ['left', 'center', 'right']});
    this.fitArea = this.bool('fitArea', true);
    this.shrinkEmptyTail = this.bool('shrinkEmptyTail', true);
    this.skipEmptyPositions = this.bool('skipEmptyPositions', false);
    this.positionMarginState = this.string('positionMarginState', 'auto',
      {choices: ['auto', 'enable', 'off']});
    let defaultValueForPositionMargin = 0;
    if (this.positionMarginState === 'auto') {
      defaultValueForPositionMargin = 4;
    }
    this.positionMargin = this.int('positionMargin', defaultValueForPositionMargin, {min: 0, max: 16});
    this.positionHeight = this.string('positionHeight', PositionHeight.full, {choices: [PositionHeight.full, PositionHeight.Entropy]});
  }

  private async init(): Promise<void> {
    if (this.initialized) {
      console.error('WebLogo second initialization!');
      return;
    }

    this.initialized = true;
    this.helpUrl = '/help/visualize/viewers/web-logo.md';

    this.msgHost = ui.div('No message');
    this.msgHost.style.display = 'none';

    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';

    this.host = ui.div([this.msgHost, this.canvas]);

    this.slider = ui.rangeSlider(0, 100, 0, 100);
    this.slider.root.style.width = '100%';
    this.slider.root.style.height = '12px';
    this.slider.root.focus();

    const getMonomer = (p: DG.Point): [number, string | null, PositionMonomerInfo | null] => {
      const jPos = Math.floor(p.x / this.positionWidthWithMargin);
      const position = this.positions[jPos];

      if (position === void 0)
        return [jPos, null, null];

      const monomer: string | undefined = Object.keys(position.freq)
        .find((m) => position.freq[m].bounds.contains(p.x, p.y));
      if (monomer === undefined)
        return [jPos, null, null];

      return [jPos, monomer, position.freq[monomer]];
    };

    // this.subs.push(rxjs.fromEvent<MouseEvent>(this.canvas, 'mouseover').subscribe((e: MouseEvent) => {
    // }));
    const correctMonomerFilter = (iRow: number, monomer: string, jPos: number) => {
      const seq = this.seqCol!.get(iRow);
      const seqM = seq ? this.splitter!(seq)[this.startPosition + jPos] : null;
      return ((seqM === monomer) || (seqM === '' && monomer === '-')) && this.dataFrame.filter.get(iRow);
    };

    this.viewSubs.push(rxjs.fromEvent<MouseEvent>(this.canvas, 'mousemove').subscribe((e: MouseEvent) => {
      if (!this.canvas)
        return;

      const args = e as MouseEvent;

      const r: number = window.devicePixelRatio;
      const cursorP: DG.Point = this.canvas.getCursorPosition(args, r);
      const [jPos, monomer] = getMonomer(cursorP);
      if (this.dataFrame && this.seqCol && this.splitter && monomer) {
        const rowCount = wu.count().take(this.dataFrame.rowCount).filter(function(iRow) {
          return correctMonomerFilter(iRow, monomer, jPos);
        }).reduce<number>((count, iRow) => count + 1, 0);
        ui.tooltip.show(ui.div([ui.div(`${monomer}`), ui.div(`${rowCount} rows`)]), args.x + 16, args.y + 16);
      } else {
        ui.tooltip.hide();
      }
    }));

    this.viewSubs.push(rxjs.fromEvent<MouseEvent>(this.canvas, 'mousedown').subscribe((e: MouseEvent) => {
      if (!this.canvas || e.button != 0)
        return;

      const args = e as MouseEvent;
      const r: number = window.devicePixelRatio;
      const [jPos, monomer] = getMonomer(this.canvas.getCursorPosition(args, r));

      // prevents deselect all rows if we miss monomer bounds
      if (this.dataFrame && this.seqCol && this.splitter && monomer) {
        this.dataFrame.selection.init(function(iRow) {
          return correctMonomerFilter(iRow, monomer, jPos);
        });
      }
    }));

    this.viewSubs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));

    this.root.append(this.host);
    this.root.appendChild(this.slider.root);

    this.render(true);
  }

  /** Handler of changing size WebLogo */
  private rootOnSizeChanged(): void {
    this.render(true);
  }

  /** Assigns {@link seqCol} and {@link cp} based on {@link sequenceColumnName} and calls {@link render}().
   */
  private updateSeqCol(): void {
    if (this.dataFrame) {
      this.seqCol = this.sequenceColumnName ? this.dataFrame.col(this.sequenceColumnName) : null;
      if (this.seqCol == null) {
        this.seqCol = WebLogo.pickUpSeqCol(this.dataFrame);
        this.sequenceColumnName = this.seqCol ? this.seqCol.name : null;
      }
      if (this.seqCol) {
        const units: string = this.seqCol!.getTag(DG.TAGS.UNITS);
        const separator: string = this.seqCol!.getTag(UnitsHandler.TAGS.separator);
        this.splitter = WebLogo.getSplitter(units, separator);
        this.unitsHandler = new UnitsHandler(this.seqCol);

        this.updatePositions();
        this.cp = WebLogo.pickUpPalette(this.seqCol);
      } else {
        this.splitter = null;
        this.positionNames = [];
        this.startPosition = -1;
        this.endPosition = -1;
        this.cp = null;
      }
    }
    this.render();
  }

  /** Updates {@link positionNames} and calculates {@link startPosition} and {@link endPosition}.
   */
  private updatePositions(): void {
    if (!this.seqCol)
      return;

    let categories: (string | null) [];
    if (this.shrinkEmptyTail) {
      const indices: Int32Array = this.dataFrame.filter.getSelectedIndexes();
      categories = Array.from(new Set(
        Array.from(Array(indices.length).keys()).map((i: number) => this.seqCol!.get(indices[i]))));
    } else {
      categories = this.seqCol.categories;
    }
    const maxLength = categories.length > 0 ? Math.max(...categories.map(
      (s) => s !== null ? this.splitter!(s).length : 0)) : 0;

    // Get position names from data column tag 'positionNames'
    const positionNamesTxt = this.seqCol.getTag('positionNames');
    // Fallback if 'positionNames' tag is not provided
    this.positionNames = positionNamesTxt ? positionNamesTxt.split(', ').map((n) => n.trim()) :
      [...Array(maxLength).keys()].map((jPos) => `${jPos + 1}`);

    this.startPosition = (this.startPositionName && this.positionNames &&
      this.positionNames.includes(this.startPositionName)) ?
      this.positionNames.indexOf(this.startPositionName) : 0;
    this.endPosition = (this.endPositionName && this.positionNames &&
      this.positionNames.includes(this.endPositionName)) ?
      this.positionNames.indexOf(this.endPositionName) : (maxLength - 1);
  }


  /** Handler of property change events. */
  public override onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);

    switch (property.name) {
    case 'considerNullSequences':
      this.render(true);
      break;
    case 'sequenceColumnName':
      this.updateSeqCol();
      break;
    case 'startPositionName':
      this.updateSeqCol();
      break;
    case 'endPositionName':
      this.updateSeqCol();
      break;
    case 'positionWidth':
      this._positionWidth = this.positionWidth;
      this.render(true);
      break;
    case 'minHeight':
      this.render(true);
      break;
    case 'maxHeight':
      this.render(true);
      break;
    case 'fixWidth':
      this.render(true);
      break;
    case 'verticalAlignment':
      this.render(true);
      break;
    case 'horizontalAlignment':
      this.render(true);
      break;
    case 'fitArea':
      this.render(true);
      break;
    case 'shrinkEmptyTail':
      this.updatePositions();
      this.render(true);
      break;
    case 'skipEmptyPositions':
      this.updatePositions();
      this.render(true);
      break;
    case 'positionMargin':
      this.render(true);
      break;
    case 'positionMarginState':
      this.render(true);
      break;
    case 'positionHeight':
      this.render(true);
      break;
    }
  }

  /** Add filter handlers when table is a attached  */
  public override async onTableAttached() {
    const dataFrameTxt: string = this.dataFrame ? 'data' : 'null';
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached( dataFrame = ${dataFrameTxt} ) start`);
    super.onTableAttached();

    this.updateSeqCol();

    if (this.dataFrame !== void 0) {
      this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
      this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => {
        this.updatePositions();
        this.render();
      }));
    }

    await this.init();
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached() end`);
  }

  /** Remove all handlers when table is a detach  */
  public override async detach() {
    const dataFrameTxt = `${this.dataFrame ? 'data' : 'null'}`;
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached( dataFrame = ${dataFrameTxt} ) start`);
    super.detach();

    this.viewSubs.forEach((sub) => sub.unsubscribe());
    this.host!.remove();
    if (this.canvas != null) {
      this.canvas = undefined;
    }
    this.msgHost = undefined;
    this.host = undefined;

    this.initialized = false;
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached() end`);
  }

  /** Helper function for rendering */
  protected _nullSequence(fillerResidue = 'X'): string {
    if (this.considerNullSequences)
      return new Array(this.Length).fill(fillerResidue).join('');

    return '';
  }

  /** Helper function for remove empty positions */
  // TODO: use this function in from core
  protected removeWhere(array: Array<any>, predicate: (T: any) => boolean): Array<any> {
    let length = array.length;
    let updateIterator = 0;
    for (let deleteIterator = 0; deleteIterator < length; deleteIterator++) {
      if (!predicate(array[deleteIterator])) {
        array[updateIterator] = array[deleteIterator];
        updateIterator++;
      }
    }
    array.length = updateIterator;
    return array;
  }

  /** Function for removing empty positions */
  protected _removeEmptyPositions() {
    if (this.skipEmptyPositions) {
      this.removeWhere(this.positions, item => item?.freq['-']?.count === item.rowCount);
    }
  }

  protected _calculate(r: number) {
    if (!this.canvas || !this.host || !this.seqCol || !this.dataFrame)
      return;
    this.unitsHandler = new UnitsHandler(this.seqCol);

    this.calcSize();

    this.positions = new Array(this.startPosition <= this.endPosition ? this.endPosition - this.startPosition + 1 : 0);
    for (let jPos = 0; jPos < this.Length; jPos++) {
      const posName: string = this.positionNames[this.startPosition + jPos];
      this.positions[jPos] = new PositionInfo(posName);
    }

    // 2022-05-05 askalkin instructed to show WebLogo based on filter (not selection)
    const indices = this.dataFrame.filter.getSelectedIndexes();
    // const indices = this.dataFrame.selection.trueCount > 0 ? this.dataFrame.selection.getSelectedIndexes() :
    //   this.dataFrame.filter.getSelectedIndexes();

    this.rowsMasked = indices.length;
    this.rowsNull = 0;

    for (const i of indices) {
      let s: string = <string>(this.seqCol.get(i));

      if (!s) {
        s = this._nullSequence();
        ++this.rowsNull;
      }

      const seqM: string[] = this.splitter!(s);
      for (let jPos = 0; jPos < this.Length; jPos++) {
        const pmInfo = this.positions[jPos].freq;
        const m: string = seqM[this.startPosition + jPos] || '-';
        if (!(m in pmInfo))
          pmInfo[m] = new PositionMonomerInfo();
        pmInfo[m].count++;
      }
    }

    //#region Polish freq counts
    for (let jPos = 0; jPos < this.Length; jPos++) {
      // delete this.positions[jPos].freq['-'];

      this.positions[jPos].rowCount = 0;
      for (const m in this.positions[jPos].freq)
        this.positions[jPos].rowCount += this.positions[jPos].freq[m].count;
      if (this.positionHeight == PositionHeight.Entropy) {
        this.positions[jPos].sumForHeightCalc = 0;
        for (const m in this.positions[jPos].freq) {
          const pn = this.positions[jPos].freq[m].count / this.positions[jPos].rowCount;
          this.positions[jPos].sumForHeightCalc += -pn * Math.log2(pn);
        }
      }
    }
    //#endregion
    this._removeEmptyPositions();

    const absoluteMaxHeight = this.canvas.height - this.axisHeight * r;

    //#region Calculate screen
    for (let jPos = 0; jPos < this.Length; jPos++) {
      const freq: { [c: string]: PositionMonomerInfo } = this.positions[jPos].freq;
      const rowCount = this.positions[jPos].rowCount;
      const alphabetSize = this.getAlphabetSize();
      if ((this.positionHeight == PositionHeight.Entropy) && (alphabetSize == null)) {
        grok.shell.error('WebLogo: alphabet is undefined.');
      }

      const maxHeight = (this.positionHeight == PositionHeight.Entropy) ? (absoluteMaxHeight * (Math.log2(alphabetSize) - (this.positions[jPos].sumForHeightCalc)) / Math.log2(alphabetSize)) : absoluteMaxHeight;

      let y: number = this.axisHeight * r + (absoluteMaxHeight - maxHeight - 1);

      const entries = Object.entries(freq).sort((a, b) => {
        if (a[0] !== '-' && b[0] !== '-')
          return b[1].count - a[1].count;
        else if (a[0] === '-' && b[0] === '-')
          return 0;
        else if (a[0] === '-')
          return -1;
        else /* (b[0] === '-') */
          return +1;
      });
      for (const entry of entries) {
        const pmInfo: PositionMonomerInfo = entry[1];
        // const m: string = entry[0];
        const h: number = maxHeight * pmInfo.count / rowCount;

        pmInfo.bounds = new DG.Rect(jPos * this.positionWidthWithMargin, y, this._positionWidth, h);
        y += h;
      }
    }
    //#endregion

  }

  /** Render WebLogo sensitive to changes in params of rendering
   *@param {boolean} recalc - indicates that need to recalculate data for rendering
   */
  render(recalc = true) {
    if (this.msgHost) {
      if (this.seqCol && !this.cp) {
        this.msgHost!.innerText = `Unknown palette (column semType: '${this.seqCol.semType}').`;
        this.msgHost!.style.display = '';
      } else {
        this.msgHost!.style.display = 'none';
      }
    }

    if (!this.canvas || !this.seqCol || !this.dataFrame || !this.cp || this.startPosition === -1 || this.endPosition === -1 || this.host == null)
      return;

    const g = this.canvas.getContext('2d');
    if (!g) return;

    const r = window.devicePixelRatio;

    if (recalc)
      this._calculate(r);

    g.resetTransform();
    g.fillStyle = 'white';
    g.fillRect(0, 0, this.canvas.width, this.canvas.height);
    const maxCountOfRowsRendered = this.host.clientWidth / this.positionWidthWithMargin + 2;
    const startRendering = Math.floor(this.host.scrollLeft / this.positionWidthWithMargin);
    const checkEndRendering = (jPos: number) => {
      return jPos < this.Length && jPos < maxCountOfRowsRendered + startRendering;
    };
    g.textBaseline = this.textBaseline;

    //#region Plot positionNames
    const positionFontSize = 10 * r;
    g.resetTransform();
    g.fillStyle = 'black';
    g.textAlign = 'center';
    g.font = `${positionFontSize.toFixed(1)}px Roboto, Roboto Local, sans-serif`;
    const posNameMaxWidth = Math.max(...this.positions.map((pos) => g.measureText(pos.name).width));
    const hScale = posNameMaxWidth < (this._positionWidth - 2) ? 1 : (this._positionWidth - 2) / posNameMaxWidth;

    for (let jPos = startRendering; checkEndRendering(jPos); jPos++) {
      const pos: PositionInfo = this.positions[jPos];
      g.resetTransform();
      g.setTransform(
        hScale, 0, 0, 1,
        jPos * this.positionWidthWithMargin + this._positionWidth / 2, 0);
      g.fillText(pos.name, 0, 0);
    }
    //#endregion Plot positionNames
    const fontStyle = '16px Roboto, Roboto Local, sans-serif';
    // Hacks to scale uppercase characters to target rectangle
    const uppercaseLetterAscent = 0.25;
    const uppercaseLetterHeight = 12.2;
    for (let jPos = startRendering; checkEndRendering(jPos); jPos++) {
      for (const [monomer, pmInfo] of Object.entries(this.positions[jPos].freq)) {
        if (monomer !== '-') {
          const monomerTxt = WebLogo.monomerToShort(monomer, 5);
          const b = pmInfo.bounds;

          g.resetTransform();
          g.strokeStyle = 'lightgray';
          g.lineWidth = 1;
          g.rect(b.left, b.top, b.width, b.height);
          g.fillStyle = this.cp.get(monomer) ?? this.cp.get('other');
          g.textAlign = 'left';
          g.font = fontStyle;
          //g.fillRect(b.left, b.top, b.width, b.height);
          const mTm: TextMetrics = g.measureText(monomerTxt);

          g.setTransform(
            b.width / mTm.width, 0, 0, b.height / uppercaseLetterHeight,
            b.left, b.top);
          g.fillText(monomerTxt, 0, -uppercaseLetterAscent);
        }
      }
    }
  }

  /** Calculate canvas size an positionWidth and updates properties */
  private calcSize() {
    if (!this.canvas || !this.host)
      return;

    const r: number = window.devicePixelRatio;

    let width: number = this.Length * this.positionWidthWithMargin / r;
    let height = Math.min(this.maxHeight, Math.max(this.minHeight, this.root.clientHeight));

    if (this.fitArea) {
      const xScale: number = this.root.clientHeight / height;
      const yScale: number = this.root.clientWidth / width;
      const scale = Math.max(1, Math.min(xScale, yScale));
      width = width * scale;
      height = height * scale;
      this._positionWidth = this.positionWidth * scale;
    }

    this.canvas.width = width * r;
    this.canvas.style.width = `${width}px`;

    // const canvasHeight: number = width > this.root.clientWidth ? height - 8 : height;
    this.host.style.setProperty('height', `${height}px`);
    const canvasHeight: number = this.host.clientHeight;
    this.canvas.height = canvasHeight * r;
    this.canvas.style.setProperty('height', `${canvasHeight}px`);

    // Adjust host and root width
    if (this.fixWidth) {
      // full width for canvas host and root
      this.root.style.width = this.host.style.width = `${width}px`;
      this.root.style.height = `${height}px`;
      this.host.style.setProperty('overflow', 'hidden', 'important');
    } else {
      // allow scroll canvas in root
      this.root.style.width = this.host.style.width = '100%';
      this.host.style.overflowX = 'auto!important';
      this.host.style.setProperty('overflow', null);
      this.host.style.setProperty('text-align', this.horizontalAlignment);

      // vertical alignment
      let hostTopMargin = 0;
      switch (this.verticalAlignment) {
      case 'top':
        hostTopMargin = 0;
        break;
      case 'middle':
        hostTopMargin = Math.max(0, (this.root.clientHeight - height) / 2);
        break;
      case 'bottom':
        hostTopMargin = Math.max(0, this.root.clientHeight - height);
        break;
      }
      this.host.style.setProperty('margin-top', `${hostTopMargin}px`, 'important');

      if (this.root.clientHeight < height) {
        this.host.style.setProperty('height', `${this.root.clientHeight}px`);
        this.host.style.setProperty('overflow-y', null);
      } else {
        this.host.style.setProperty('overflow-y', 'hidden', 'important');
      }
    }
  }

  /** Selects a suitable palette based on column data
   * @param {DG.Column} seqCol Column to look for a palette
   * @param {number}  minLength minimum length of sequence to detect palette (empty strings are allowed)
   * @return {SeqPalette} Palette corresponding to the alphabet of the sequences in the column
   */
  public static pickUpPalette(seqCol: DG.Column, minLength: number = 5): SeqPalette {
    let res: SeqPalette | null;
    switch (seqCol.semType) {
    case Aminoacids.SemTypeMultipleAlignment:
      res = AminoacidsPalettes.GrokGroups;
      break;
    case Nucleotides.SemTypeMultipleAlignment:
      res = NucleotidesPalettes.Chromatogram;
      break;
    }
    const stats: SeqColStats = WebLogo.getStats(seqCol, minLength, WebLogo.splitterAsFasta);

    const alphabetCandidates: [Set<string>, SeqPalette][] = [
      [new Set(Object.keys(Nucleotides.Names)), NucleotidesPalettes.Chromatogram],
      [new Set(Object.keys(Aminoacids.Names)), AminoacidsPalettes.GrokGroups],
    ];
    // Calculate likelihoods for alphabet_candidates
    const alphabetCandidatesSim: number[] = alphabetCandidates
      .map((c) => WebLogo.getAlphabetSimilarity(stats.freq, c[0]));
    const maxCos = Math.max(...alphabetCandidatesSim);
    if (maxCos > 0.55)
      res = alphabetCandidates[alphabetCandidatesSim.indexOf(maxCos)][1];
    else
      res = UnknownSeqPalettes.Color;

    return res;
  }

  public getAlphabetSize(): number {
    return this.unitsHandler?.getAlphabetSize() ?? 0;
  }

  /** Stats of sequences with specified splitter func, returns { freq, sameLength }.
   * @param {DG.Column} seqCol
   * @param {number} minLength
   * @param {SplitterFunc} splitter
   * @return { SeqColStats }, sameLength: boolean } stats of column sequences
   */
  static getStats(seqCol: DG.Column, minLength: number, splitter: SplitterFunc): SeqColStats {
    const freq: { [m: string]: number } = {};
    let sameLength = true;
    let firstLength = null;

    for (const seq of seqCol.categories) {
      const mSeq = splitter(seq);

      if (firstLength == null)
        firstLength = mSeq.length;
      else if (mSeq.length !== firstLength)
        sameLength = false;

      if (mSeq.length > minLength) {
        for (const m of mSeq) {
          if (!(m in freq))
            freq[m] = 0;
          freq[m] += 1;
        }
      }
    }
    return {freq: freq, sameLength: sameLength};
  }

  /** Calculate similarity in current sequence and alphabet.
   * @param {MonomerFreqs} freq
   * @param {Set<string>} alphabet
   * @param {string} gapSymbol
   * @return {number} Cosine similarity
   */
  public static getAlphabetSimilarity(freq: MonomerFreqs, alphabet: Set<string>, gapSymbol: string = '-'): number {
    const keys = new Set<string>([...new Set(Object.keys(freq)), ...alphabet]);
    keys.delete(gapSymbol);

    const freqA: number[] = [];
    const alphabetA: number[] = [];
    for (const m of keys) {
      freqA.push(m in freq ? freq[m] : 0);
      alphabetA.push(alphabet.has(m) ? 1 : 0);
    }
    /* There were a few ideas: chi-squared, pearson correlation (variance?), scalar product */
    const freqV: Vector = new Vector(freqA);
    const alphabetV: Vector = new Vector(alphabetA);
    return vectorDotProduct(freqV, alphabetV) / (vectorLength(freqV) * vectorLength(alphabetV));
  }

  // Previous version of pickUpSeqCol
  // /** First try to find column with semType 'alignedSequence'.
  //  * Next look for column with data alphabet corresponding to any of the known palettes.
  //  * @param {DG.DataFrame} dataFrame
  //  * @return {DG.Column} The column we were looking for or null
  //  */
  // private pickUpSeqCol(dataFrame: DG.DataFrame): DG.Column | null {
  //   let res: DG.Column | null = dataFrame.columns.bySemType('alignedSequence');
  //   if (res == null) {
  //     for (const col of dataFrame.columns) {
  //       const cp = WebLogo.pickUpPalette(col as DG.Column, 5);
  //       if (cp !== null && !(cp instanceof UnknownSeqPalette)) {
  //         res = col;
  //         break;
  //       }
  //     }
  //   }
  //   return res;
  // }

  public static pickUpSeqCol(df: DG.DataFrame): DG.Column | null {
    const semTypeColList = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
    let resCol: DG.Column | null = semTypeColList.find((col) => {
      const units = col.getTag(DG.TAGS.UNITS);
      return units ? units.indexOf('MSA') !== -1 : false;
    }) ?? null;
    if (!resCol && semTypeColList.length > 0)
      resCol = semTypeColList[0];
    return resCol;
  }

  private static monomerRe = /\[(\w+)|(\w)|(-)/g;

  /** Only some of the synonyms. These were obtained from the clustered oligopeptide dataset. */
  private static aaSynonyms: { [name: string]: string } = {
    '[MeNle]': 'L', // Nle - norleucine
    '[MeA]': 'A',
    '[MeG]': 'G',
    '[MeF]': 'F',
  };

  /** Split sequence for single character monomers, square brackets multichar monomer names or gap symbol.
   * @param {any} seq object with sequence
   * @return {string[]} array of monomers
   */
  public static splitterAsFasta(seq: any): string[] {
    return wu<RegExpMatchArray>(seq.toString().matchAll(WebLogo.monomerRe))
      .map((ma: RegExpMatchArray) => {
        let mRes: string;
        const m: string = ma[0];
        if (m.length > 1) {
          if (m in WebLogo.aaSynonyms) {
            mRes = WebLogo.aaSynonyms[m];
          } else {
            mRes = '';
            console.debug(`Long monomer '${m}' has not a short synonym.`);
          }
        } else {
          mRes = m;
        }
        return mRes;
      }).toArray();
  }

  private static helmRe = /(PEPTIDE1|DNA1|RNA1)([^}]+)}/g;
  private static helmPp1Re = /\[([^\[\]]+)]/g;

  /** Splits Helm string to monomers, but does not replace monomer names to other notation (e.g. for RNA).
   * Only for linear polymers, does not split RNA for ribose and phosphate monomers.
   * @param {string} seq Source string of HELM notation
   * @return {string[]}
   */
  public static splitterAsHelm(seq: any): string[] {
    WebLogo.helmRe.lastIndex = 0;
    const ea: RegExpExecArray | null = WebLogo.helmRe.exec(seq.toString());
    const inSeq: string | null = ea ? ea[2] : null;

    const mmPostProcess = (mm: string): string => {
      WebLogo.helmPp1Re.lastIndex = 0;
      const pp1M = WebLogo.helmPp1Re.exec(mm);
      if (pp1M && pp1M.length >= 2)
        return pp1M[1];
      else
        return mm;
    };

    const mmList: string[] = inSeq ? inSeq.split('.') : [];
    return mmList.map(mmPostProcess);
  }

  /** Gets method to split sequence by separator
   * @param {string} separator
   * @param limit
   * @return {SplitterFunc}
   */
  public static getSplitterWithSeparator(separator: string, limit: number | undefined = undefined): SplitterFunc {
    return (seq: string) => {
      return seq.split(separator, limit);
    };
  }

  /** Get splitter method to split sequences to monomers
   * @param {string} units
   * @param {string} separator
   * @param limit
   * @return {SplitterFunc}
   */
  public static getSplitter(units: string, separator: string, limit: number | undefined = undefined): SplitterFunc {
    if (units.toLowerCase().startsWith(NOTATION.FASTA))
      return WebLogo.splitterAsFasta;
    else if (units.toLowerCase().startsWith(NOTATION.SEPARATOR))
      return WebLogo.getSplitterWithSeparator(separator, limit);
    else if (units.toLowerCase().startsWith(NOTATION.HELM))
      return WebLogo.splitterAsHelm;
    else
      throw new Error(`Unexpected units ${units} .`);

    // TODO: Splitter for HELM
  }

  /** Generate splitter function for sequence column
   * @param {DG.Column} col
   * @return {SplitterFunc} Splitter function
   */
  public static getSplitterForColumn(col: DG.Column): SplitterFunc {
    if (col.semType !== DG.SEMTYPE.MACROMOLECULE)
      throw new Error(`Get splitter for semType "${DG.SEMTYPE.MACROMOLECULE}" only.`);

    const units = col.getTag(DG.TAGS.UNITS);
    const separator = col.getTag(UnitsHandler.TAGS.separator);
    return WebLogo.getSplitter(units, separator);
  }

  /** Convert long monomer names to short ones */
  public static monomerToShort(amino: string, maxLengthOfMonomer: number): string {
    return amino.length <= maxLengthOfMonomer ? amino : amino.substring(0, maxLengthOfMonomer) + '…';
  }
}
