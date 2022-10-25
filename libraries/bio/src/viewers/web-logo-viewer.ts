import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import * as rxjs from 'rxjs';

import {SeqPalette} from '../seq-palettes';
import {Subscription} from 'rxjs';
import {UnitsHandler} from '../utils/units-handler';
import {SliderOptions} from 'datagrok-api/dg';
import {getSplitter, monomerToShort, pickUpPalette, pickUpSeqCol, SplitterFunc, TAGS} from '../utils/macromolecule';

export enum PositionHeight {
  Entropy = 'Entropy',
  full = '100%',
}

declare global {
  interface HTMLCanvasElement {
    getCursorPosition(event: MouseEvent, r: number): DG.Point;
  }
}

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

export class WebLogoViewer extends DG.JsViewer {
  public static residuesSet = 'nucleotides';
  private static viewerCount: number = -1;

  private readonly viewerId: number = -1;
  private unitsHandler: UnitsHandler | null;
  private initialized: boolean = false;

  // private readonly colorScheme: ColorScheme = ColorSchemes[NucleotidesWebLogo.residuesSet];
  protected cp: SeqPalette | null = null;

  private host?: HTMLDivElement;
  private msgHost?: HTMLElement;
  private canvas: HTMLCanvasElement;
  private slider: DG.RangeSlider;
  private readonly textBaseline: CanvasTextBaseline;

  private axisHeight: number = 12;

  private seqCol: DG.Column<string> | null = null;
  private splitter: SplitterFunc | null = null;
  // private maxLength: number = 100;
  private positions: PositionInfo[] = [];

  private rowsMasked: number = 0;
  private rowsNull: number = 0;
  private visibleSlider: boolean = false;
  private allowResize: boolean = true;
  private turnOfResizeForOneSetValue: boolean = false;

  // Viewer's properties (likely they should be public so that they can be set outside)
  private _positionWidth: number;
  public positionWidth: number;
  public minHeight: number;
  public backgroundColor: number = 0xFFFFFFFF;
  public maxHeight: number;
  public skipEmptySequences: boolean;
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
    return this._positionWidth + this.positionMarginValue;
  }

  private get positionMarginValue() {
    if ((this.positionMarginState === 'auto') && (this.unitsHandler?.getAlphabetIsMultichar() === true)) {
      return this.positionMargin;
    }
    if (this.positionMarginState === 'enable') {
      return this.positionMargin;
    }

    return 0;
  }

  /** Count of position rendered for calculations countOfRenderPositions */
  private get countOfRenderPositions() {
    if (this.host == null) {
      return 0;
    }
    const r = window.devicePixelRatio;
    if (r > 1) {
      return this.canvasWidthWithRatio / this.positionWidthWithMargin;
    } else {
      return this.canvas.width / (this.positionWidthWithMargin * r);
    }
  }

  private get canvasWidthWithRatio() {
    return this.canvas.width * window.devicePixelRatio;
  }


  /** Position of start rendering */
  private get firstVisibleIndex(): number {
    return (this.visibleSlider) ? Math.floor(this.slider.min) : 0;
  }

  private viewSubs: Subscription[] = [];

  constructor() {
    super();

    this.viewerId = WebLogoViewer.viewerCount;
    WebLogoViewer.viewerCount += 1;

    this.textBaseline = 'top';
    this.unitsHandler = null;

    this.backgroundColor = this.int('backgroundColor', 0xFFFFFFFF);
    this._positionWidth = this.positionWidth = this.float('positionWidth', 16/*,
      {editor: 'slider', min: 4, max: 64, postfix: 'px'}*/);
    this.minHeight = this.float('minHeight', 50/*,
      {editor: 'slider', min: 25, max: 250, postfix: 'px'}*/);
    this.maxHeight = this.float('maxHeight', 100/*,
      {editor: 'slider', min: 25, max: 500, postfix: 'px'}*/);

    this.skipEmptySequences = this.bool('skipEmptySequences', true);
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

    const style: SliderOptions = {style: 'barbell'};
    this.slider = ui.rangeSlider(0, 100, 0, 20, false, style);
    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';
  }

  private init(): void {
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

    //this.slider.setShowHandles(false);
    this.slider.root.style.position = 'absolute';
    this.slider.root.style.zIndex = '999';
    this.slider.root.style.display = 'none';
    this.slider.root.style.height = '0.7em';

    this.visibleSlider = false;

    this.slider.onValuesChanged.subscribe(() => {
      if ((this.host == null)) {
        return;
      }
      /* Resize slider if we can resize do that */
      if ((this.allowResize) && (!this.turnOfResizeForOneSetValue) &&
        (this.visibleSlider)) {
        const countOfPositions = Math.ceil(this.slider.max - this.slider.min);
        const calculatedWidth = (this.canvas.width / countOfPositions) - this.positionMarginValue;
        // saving positionWidth value global (even if slider is not visible)
        this.positionWidth = calculatedWidth;
        this._positionWidth = calculatedWidth;
      }
      this.turnOfResizeForOneSetValue = false;
      this.render(true);
    });


    this.host = ui.div([this.msgHost, this.canvas]);

    this.host.style.justifyContent = 'center';
    this.host.style.alignItems = 'center';
    this.host.style.position = 'relative';
    this.host.style.setProperty('overflow', 'hidden', 'important');

    const getMonomer = (p: DG.Point): [number, string | null, PositionMonomerInfo | null] => {
      const calculatedX = p.x + this.firstVisibleIndex * this.positionWidthWithMargin;
      const jPos = Math.floor(p.x / this.positionWidthWithMargin + this.firstVisibleIndex);
      const position = this.positions[jPos];

      if (position === void 0)
        return [jPos, null, null];

      const monomer: string | undefined = Object.keys(position.freq)
        .find((m) => position.freq[m].bounds.contains(calculatedX, p.y));
      if (monomer === undefined)
        return [jPos, null, null];

      return [jPos, monomer, position.freq[monomer]];
    };

    const correctMonomerFilter = (iRow: number, monomer: string, jPos: number) => {
      const seq = this.seqCol!.get(iRow);
      const seqM = seq ? this.splitter!(seq)[this.startPosition + jPos] : null;
      return ((seqM === monomer) || (seqM === '' && monomer === '-')) && this.dataFrame.filter.get(iRow);
    };

    rxjs.fromEvent<MouseEvent>(this.canvas, 'mousemove').subscribe((e: MouseEvent) => {
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
    });

    rxjs.fromEvent<MouseEvent>(this.canvas, 'mousedown').subscribe((e: MouseEvent) => {
      const args = e as MouseEvent;
      const r: number = window.devicePixelRatio;
      const [jPos, monomer] = getMonomer(this.canvas.getCursorPosition(args, r));

      // prevents deselect all rows if we miss monomer bounds
      if (this.dataFrame && this.seqCol && this.splitter && monomer) {
        this.dataFrame.selection.init(function(iRow) {
          return correctMonomerFilter(iRow, monomer, jPos);
        });
      }
    });

    rxjs.fromEvent<WheelEvent>(this.canvas, 'wheel').subscribe((e: WheelEvent) => {
      if (!this.visibleSlider)
        return;
      const countOfScrollPositions = (e.deltaY / 100) * Math.max(Math.floor((this.countOfRenderPositions) / 2), 1);
      this.slider.scrollBy(this.slider.min + countOfScrollPositions);

    });

    this.viewSubs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));

    this.root.append(this.host);
    this.root.append(this.slider.root);

    this._calculate(window.devicePixelRatio);
    this.updateSlider();
    this.render(true);
  }

  /** Handler of changing size WebLogo */
  private rootOnSizeChanged(): void {
    this._calculate(window.devicePixelRatio);
    this.updateSlider();
    this.render(true);
  }

  /** Assigns {@link seqCol} and {@link cp} based on {@link sequenceColumnName} and calls {@link render}().
   */
  private updateSeqCol(): void {
    if (this.dataFrame) {
      this.seqCol = this.sequenceColumnName ? this.dataFrame.col(this.sequenceColumnName) : null;
      if (this.seqCol == null) {
        this.seqCol = pickUpSeqCol(this.dataFrame);
        this.sequenceColumnName = this.seqCol ? this.seqCol.name : null;
      }
      if (this.seqCol) {
        const units: string = this.seqCol!.getTag(DG.TAGS.UNITS);
        const separator: string = this.seqCol!.getTag(TAGS.separator);
        this.splitter = getSplitter(units, separator);
        this.unitsHandler = new UnitsHandler(this.seqCol);

        this.updatePositions();
        this.cp = pickUpPalette(this.seqCol);
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

  private get widthArea() {
    return this.Length * this.positionWidth / window.devicePixelRatio;
  }

  private get heightArea() {
    return Math.min(this.maxHeight, Math.max(this.minHeight, this.root.clientHeight));
  }

  private get xScale() {
    return this.widthArea > 0 ? (this.root.clientWidth - this.Length * this.positionMarginValue) / this.widthArea : 0;
  }

  private get yScale() {
    return this.root.clientHeight / this.heightArea;
  }

  private checkIsHideSlider(): boolean {
    let showSliderWithFitArea = true;
    const minScale = Math.min(this.xScale, this.yScale);

    if (((minScale == this.xScale) || (minScale <= 1)) && (this.fitArea)) {
      showSliderWithFitArea = false;
    }
    return ((this.fixWidth || Math.ceil(this.canvas.width / this.positionWidthWithMargin) >= this.Length) || (showSliderWithFitArea));
  }

  setSliderVisibility(visible: boolean): void {
    if (visible) {
      this.slider.root.style.display = 'inherit';
      this.visibleSlider = true;
    } else {
      this.slider.root.style.display = 'none';
      this.visibleSlider = false;
    }
  }

  /** Updates {@link slider}, needed to set slider options and to update slider position. */
  private updateSlider(): void {
    if (this.checkIsHideSlider()) {
      this.setSliderVisibility(false);
    } else {
      this.setSliderVisibility(true);
    }
    if ((this.slider != null) && (this.canvas != null)) {
      const diffEndScrollAndSliderMin = Math.max(0,
        Math.floor(this.slider.min + this.canvas.width / this.positionWidthWithMargin) - this.Length);
      let newMin = Math.floor(this.slider.min - diffEndScrollAndSliderMin);
      let newMax = Math.floor(this.slider.min - diffEndScrollAndSliderMin) + Math.floor(this.canvas.width / this.positionWidthWithMargin);
      if (this.checkIsHideSlider()) {
        newMin = 0;
        newMax = Math.max(newMin, this.Length - 1);
      }
      this.turnOfResizeForOneSetValue = true;
      this.slider.setValues(0, this.Length,
        newMin, newMax);
    }
  }

  /** Handler of property change events. */
  public override onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);

    switch (property.name) {
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
      this.updateSlider();
      break;
    case 'fixWidth':
      this.updateSlider();
      break;
    case 'fitArea':
      this.updateSlider();
      break;
    case 'shrinkEmptyTail':
      this.updatePositions();
      break;
    case 'skipEmptyPositions':
      this.updatePositions();
      break;
    case 'positionMargin':
      this.updateSlider();
      break;
    }

    this.render(true);
  }

  /** Add filter handlers when table is a attached  */
  public override onTableAttached() {
    super.onTableAttached();

    const dataFrameTxt: string = this.dataFrame ? 'data' : 'null';
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached( dataFrame = ${dataFrameTxt} ) start`);

    this.updateSeqCol();

    if (this.dataFrame !== void 0) {
      this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
      this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => {
        this.updatePositions();
        this.render();
      }));
    }

    this.init();
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached() end`);
  }

  /** Remove all handlers when table is a detach  */
  public override async detach() {
    const dataFrameTxt = `${this.dataFrame ? 'data' : 'null'}`;
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached( dataFrame = ${dataFrameTxt} ) start`);
    super.detach();

    this.viewSubs.forEach((sub) => sub.unsubscribe());
    this.host!.remove();
    this.msgHost = undefined;
    this.host = undefined;

    this.initialized = false;
    console.debug(`bio: WebLogo<${this.viewerId}>.onTableAttached() end`);
  }

  /** Helper function for rendering */
  protected _nullSequence(fillerResidue = 'X'): string {
    if (!this.skipEmptySequences)
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
    if (!this.host || !this.seqCol || !this.dataFrame)
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

    if (!this.seqCol || !this.dataFrame || !this.cp || this.startPosition === -1 || this.endPosition === -1 || this.host == null || this.slider == null)
      return;

    const g = this.canvas.getContext('2d');
    if (!g) return;

    this.slider.root.style.width = `${this.host.clientWidth}px`;

    const r = window.devicePixelRatio;

    if (recalc)
      this._calculate(r);

    g.resetTransform();
    g.fillStyle = DG.Color.toHtml(this.backgroundColor);
    g.fillRect(0, 0, this.canvas.width, this.canvas.height);
    g.textBaseline = this.textBaseline;

    const maxCountOfRowsRendered = this.countOfRenderPositions + 1;
    const firstVisibleIndex = (this.visibleSlider) ? Math.floor(this.slider.min) : 0;
    const lastVisibleIndex = Math.min(this.Length, firstVisibleIndex + maxCountOfRowsRendered);

    //#region Plot positionNames
    const positionFontSize = 10 * r;
    g.resetTransform();
    g.fillStyle = 'black';
    g.textAlign = 'center';
    g.font = `${positionFontSize.toFixed(1)}px Roboto, Roboto Local, sans-serif`;
    const posNameMaxWidth = Math.max(...this.positions.map((pos) => g.measureText(pos.name).width));
    const hScale = posNameMaxWidth < (this._positionWidth - 2) ? 1 : (this._positionWidth - 2) / posNameMaxWidth;

    for (let jPos = this.firstVisibleIndex; jPos < lastVisibleIndex; jPos++) {
      const pos: PositionInfo = this.positions[jPos];
      g.resetTransform();
      g.setTransform(
        hScale, 0, 0, 1,
        jPos * this.positionWidthWithMargin + this._positionWidth / 2 - this.positionWidthWithMargin * firstVisibleIndex, 0);
      g.fillText(pos.name, 0, 0);
    }
    //#endregion Plot positionNames
    const fontStyle = '16px Roboto, Roboto Local, sans-serif';
    // Hacks to scale uppercase characters to target rectangle
    const uppercaseLetterAscent = 0.25;
    const uppercaseLetterHeight = 12.2;
    for (let jPos = this.firstVisibleIndex; jPos < lastVisibleIndex; jPos++) {
      for (const [monomer, pmInfo] of Object.entries(this.positions[jPos].freq)) {
        if (monomer !== '-') {
          const monomerTxt = monomerToShort(monomer, 5);
          const b = pmInfo.bounds;
          const left = b.left - this.positionWidthWithMargin * this.firstVisibleIndex;

          g.resetTransform();
          g.strokeStyle = 'lightgray';
          g.lineWidth = 1;
          g.rect(left, b.top, b.width, b.height);
          g.fillStyle = this.cp.get(monomer) ?? this.cp.get('other');
          g.textAlign = 'left';
          g.font = fontStyle;
          //g.fillRect(b.left, b.top, b.width, b.height);
          const mTm: TextMetrics = g.measureText(monomerTxt);

          g.setTransform(
            b.width / mTm.width, 0, 0, b.height / uppercaseLetterHeight,
            left, b.top);
          g.fillText(monomerTxt, 0, -uppercaseLetterAscent);
        }
      }
    }
  }

  /** Calculate canvas size an positionWidth and updates properties */
  private calcSize() {
    if (!this.host)
      return;

    const r: number = window.devicePixelRatio;

    let width: number = this.widthArea;
    let height = this.heightArea;

    if ((this.fitArea) && (!this.visibleSlider)) {
      const scale = Math.max(1, Math.min(this.xScale, this.yScale));
      width = width * scale;
      height = height * scale;
      this._positionWidth = this.positionWidth * scale;
    }

    width = this.Length * this.positionWidthWithMargin / r;

    this.canvas.width = this.root.clientWidth * r;
    this.canvas.style.width = `${this.root.clientWidth}px`;

    // const canvasHeight: number = width > this.root.clientWidth ? height - 8 : height;
    this.host.style.setProperty('height', `${height}px`);
    const canvasHeight: number = this.host.clientHeight;
    this.canvas.height = canvasHeight * r;

    // Adjust host and root width
    if (this.fixWidth) {
      // full width for canvas host and root
      this.root.style.width = this.host.style.width = `${width}px`;
      this.root.style.height = `${height}px`;
      this.root.style.overflow = 'hidden';
      this.host.style.setProperty('overflow-y', 'hidden', 'important');
    } else {
      // allow scroll canvas in root
      this.root.style.width = this.host.style.width = '100%';
      this.host.style.overflowX = 'auto!important';
      this.host.style.setProperty('text-align', this.horizontalAlignment);

      const sliderHeight = this.visibleSlider ? 10 : 0;

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
        hostTopMargin = Math.max(0, this.root.clientHeight - height - sliderHeight);
        break;
      }
      // horizontal alignment
      let hostLeftMargin = 0;
      switch (this.horizontalAlignment) {
      case 'left':
        hostLeftMargin = 0;
        break;
      case 'center':
        hostLeftMargin = Math.max(0, (this.root.clientWidth - width) / 2);
        break;
      case 'right':
        hostLeftMargin = Math.max(0, this.root.clientWidth - width);
        break;
      }
      this.host.style.setProperty('margin-top', `${hostTopMargin}px`, 'important');
      this.host.style.setProperty('margin-left', `${hostLeftMargin}px`, 'important');
      if (this.slider != null) {
        this.slider.root.style.setProperty('margin-top', `${hostTopMargin + canvasHeight}px`, 'important');
      }

      if (this.root.clientHeight <= height) {
        this.host.style.setProperty('height', `${this.root.clientHeight}px`);
        this.host.style.setProperty('overflow-y', null);
      } else {
        this.host.style.setProperty('overflow-y', 'hidden', 'important');
      }
    }
  }

  public getAlphabetSize(): number {
    return this.unitsHandler?.getAlphabetSize() ?? 0;
  }
}
