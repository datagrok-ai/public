import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import * as rxjs from 'rxjs';

import {SliderOptions} from 'datagrok-api/dg';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {
  getSplitter, monomerToShort, pickUpPalette, pickUpSeqCol, SplitterFunc,
  TAGS as bioTAGS
} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  WebLogoPropsDefault, WebLogoProps, IWebLogoViewer,
  PositionHeight,
  positionSeparator,
} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {TAGS as wlTAGS} from '@datagrok-libraries/bio/src/viewers/web-logo';

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
  /** Position in sequence */
  public readonly pos: number;

  /** Position name from column tag*/
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
  constructor(pos: number, name: string,
    freq: { [m: string]: PositionMonomerInfo } = {}, rowCount: number = 0, sumForHeightCalc: number = 0
  ) {
    this.pos = pos;
    this.name = name;
    this.freq = freq;
    this.rowCount = rowCount;
    this.sumForHeightCalc = sumForHeightCalc;
  }
}

export enum VerticalAlignments {
  TOP = 'top',
  MIDDLE = 'middle',
  BOTTOM = 'bottom',
}

export enum HorizontalAlignments {
  LEFT = 'left',
  CENTER = 'center',
  RIGHT = 'right',
}

export enum PositionMarginStates {
  AUTO = 'auto',
  ON = 'on',
  OFF = 'off',
}

export enum FilterSources {
  Filtered = 'Filtered',
  Selected = 'Selected',
}

export enum PROPS_CATS {
  STYLE = 'Style',
  BEHAVIOR = 'Behavior',
  LAYOUT = 'Layout',
  DATA = 'Data',
}

export enum PROPS {
  // -- Data --
  sequenceColumnName = 'sequenceColumnName',
  startPositionName = 'startPositionName',
  endPositionName = 'endPositionName',
  skipEmptySequences = 'skipEmptySequences',
  skipEmptyPositions = 'skipEmptyPositions',
  shrinkEmptyTail = 'shrinkEmptyTail',

  // -- Style --
  backgroundColor = 'backgroundColor',
  positionHeight = 'positionHeight',
  positionWidth = 'positionWidth',

  // -- Layout --
  verticalAlignment = 'verticalAlignment',
  horizontalAlignment = 'horizontalAlignment',
  fixWidth = 'fixWidth',
  fitArea = 'fitArea',
  minHeight = 'minHeight',
  maxHeight = 'maxHeight',
  positionMarginState = 'positionMarginState',
  positionMargin = 'positionMargin',

  // -- Behavior --
  filterSource = 'filterSource',
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
  // -- Data --
  public sequenceColumnName: string | null;
  public skipEmptySequences: boolean;
  public skipEmptyPositions: boolean;

  // -- Style --
  private _positionWidth: number;
  public positionWidth: number;
  public minHeight: number;
  public backgroundColor: number = 0xFFFFFFFF;
  public maxHeight: number;
  public positionMarginState: string;
  public positionMargin: number = 0;
  public startPositionName: string | null;
  public endPositionName: string | null;
  public fixWidth: boolean;
  public verticalAlignment: string | null;
  public horizontalAlignment: string | null;
  public fitArea: boolean;
  public shrinkEmptyTail: boolean;
  public positionHeight: string;

  // -- Behavior --
  public filterSource: FilterSources;

  private positionNames: string[] = [];
  private startPosition: number = -1;
  private endPosition: number = -1;

  private get filter(): DG.BitSet {
    let res: DG.BitSet;
    switch (this.filterSource) {
    case FilterSources.Filtered:
      res = this.dataFrame.filter;
      break;
    case FilterSources.Selected:
      res = this.dataFrame.selection;
      break;
    }
    return res;
  }

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

  private viewSubs: rxjs.Unsubscribable[] = [];

  constructor() {
    super();

    this.viewerId = WebLogoViewer.viewerCount;
    WebLogoViewer.viewerCount += 1;

    this.textBaseline = 'top';
    this.unitsHandler = null;

    // -- Data --
    this.sequenceColumnName = this.string(PROPS.sequenceColumnName, null,
      {category: PROPS_CATS.DATA});
    this.startPositionName = this.string(PROPS.startPositionName, null,
      {category: PROPS_CATS.DATA});
    this.endPositionName = this.string(PROPS.endPositionName, null,
      {category: PROPS_CATS.DATA});
    this.skipEmptySequences = this.bool(PROPS.skipEmptySequences, true,
      {category: PROPS_CATS.DATA});
    this.skipEmptyPositions = this.bool(PROPS.skipEmptyPositions, false,
      {category: PROPS_CATS.DATA});
    this.shrinkEmptyTail = this.bool(PROPS.shrinkEmptyTail, true,
      {category: PROPS_CATS.DATA});


    // -- Style --
    this.backgroundColor = this.int(PROPS.backgroundColor, 0xFFFFFFFF,
      {category: PROPS_CATS.STYLE});
    this.positionHeight = this.string(PROPS.positionHeight, PositionHeight.full,
      {category: PROPS_CATS.STYLE, choices: Object.values(PositionHeight)});
    this._positionWidth = this.positionWidth = this.float(PROPS.positionWidth, 16,
      {category: PROPS_CATS.STYLE/* editor: 'slider', min: 4, max: 64, postfix: 'px' */});


    // -- Layout --
    this.verticalAlignment = this.string(PROPS.verticalAlignment, VerticalAlignments.MIDDLE,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(VerticalAlignments)});
    this.horizontalAlignment = this.string(PROPS.horizontalAlignment, HorizontalAlignments.CENTER,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(HorizontalAlignments)});
    this.fixWidth = this.bool(PROPS.fixWidth, false,
      {category: PROPS_CATS.LAYOUT});
    this.fitArea = this.bool(PROPS.fitArea, true,
      {category: PROPS_CATS.LAYOUT});
    this.minHeight = this.float(PROPS.minHeight, 50,
      {category: PROPS_CATS.LAYOUT/*, editor: 'slider', min: 25, max: 250, postfix: 'px'*/});
    this.maxHeight = this.float(PROPS.maxHeight, 100,
      {category: PROPS_CATS.LAYOUT/*, editor: 'slider', min: 25, max: 500, postfix: 'px'*/});
    this.positionMarginState = this.string(PROPS.positionMarginState, PositionMarginStates.AUTO,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(PositionMarginStates)});
    let defaultValueForPositionMargin = 0;
    if (this.positionMarginState === 'auto') defaultValueForPositionMargin = 4;
    this.positionMargin = this.int(PROPS.positionMargin, defaultValueForPositionMargin,
      {category: PROPS_CATS.LAYOUT, min: 0, max: 16});

    // -- Behavior --
    this.filterSource = this.string(PROPS.filterSource, FilterSources.Filtered,
      {category: PROPS_CATS.BEHAVIOR, choices: Object.values(FilterSources)}) as FilterSources;

    const style: SliderOptions = {style: 'barbell'};
    this.slider = ui.rangeSlider(0, 100, 0, 20, false, style);
    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';
  }

  private init(): void {
    if (this.initialized) {
      console.error('Bio: WebLogoViewer.init() second initialization!');
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

    this.subs.push(this.slider.onValuesChanged.subscribe(this.sliderOnValuesChanged.bind(this)));

    this.host = ui.div([this.msgHost, this.canvas]);

    this.host.style.justifyContent = 'center';
    this.host.style.alignItems = 'center';
    this.host.style.position = 'relative';
    this.host.style.setProperty('overflow', 'hidden', 'important');

    this.subs.push(
      rxjs.fromEvent<MouseEvent>(this.canvas, 'mousemove').subscribe(this.canvasOnMouseMove.bind(this)));
    this.subs.push(
      rxjs.fromEvent<MouseEvent>(this.canvas, 'mousedown').subscribe(this.canvasOnMouseDown.bind(this)));

    this.subs.push(rxjs.fromEvent<WheelEvent>(this.canvas, 'wheel').subscribe(this.canvasOnWheel.bind(this)));

    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));

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
        const separator: string = this.seqCol!.getTag(bioTAGS.separator);
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
    const positionNamesTxt = this.seqCol.getTag(wlTAGS.positionNames);
    // Fallback if 'positionNames' tag is not provided
    this.positionNames = positionNamesTxt ? positionNamesTxt.split(positionSeparator).map((n) => n.trim()) :
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
    case PROPS.sequenceColumnName:
    case PROPS.startPositionName:
    case PROPS.endPositionName:
    case PROPS.filterSource:
      this.updateSeqCol();
      break;
    case PROPS.positionWidth:
      this._positionWidth = this.positionWidth;
      this.updateSlider();
      break;
    case PROPS.fixWidth:
    case PROPS.fitArea:
    case PROPS.positionMargin:
      this.updateSlider();
      break;
    case PROPS.shrinkEmptyTail:
    case PROPS.skipEmptyPositions:
      this.updatePositions();
      break;
    }

    this.render(true);
  }

  /** Add filter handlers when table is a attached  */
  public override onTableAttached() {
    super.onTableAttached();

    const dataFrameTxt: string = this.dataFrame ? 'data' : 'null';
    console.debug(`Bio: WebLogoViewer<${this.viewerId}>.onTableAttached( dataFrame = ${dataFrameTxt} ) start`);

    this.updateSeqCol();

    if (this.dataFrame !== undefined) {
      this.viewSubs.push(this.dataFrame.filter.onChanged.subscribe(this.dataFrameFilterOnChanged.bind(this)));
      this.viewSubs.push(this.dataFrame.selection.onChanged.subscribe(this.dataFrameSelectionOnChanged.bind(this)));
    }

    this.init();
    console.debug(`Bio: WebLogoViewer<${this.viewerId}>.onTableAttached() end`);
  }

  /** Remove all handlers when table is a detach  */
  public override async detach() {
    const dataFrameTxt = `${this.dataFrame ? 'data' : 'null'}`;
    console.debug(`Bio: WebLogoViewer<${this.viewerId}>.onTableAttached( dataFrame = ${dataFrameTxt} ) start`);
    super.detach();

    this.viewSubs.forEach((sub) => sub.unsubscribe());
    this.host!.remove();
    this.msgHost = undefined;
    this.host = undefined;

    this.initialized = false;
    console.debug(`Bio: WebLogoViewer<${this.viewerId}>.onTableAttached() end`);
  }

  // -- Routines --

  getMonomer(p: DG.Point): [number, string | null, PositionMonomerInfo | null] {
    const calculatedX = p.x + this.firstVisibleIndex * this.positionWidthWithMargin;
    const jPos = Math.floor(p.x / this.positionWidthWithMargin + this.firstVisibleIndex);
    const position = this.positions[jPos];

    if (position == undefined)
      return [jPos, null, null];

    const monomer: string | undefined = Object.keys(position.freq)
      .find((m) => position.freq[m].bounds.contains(calculatedX, p.y));
    if (monomer === undefined)
      return [jPos, null, null];

    return [jPos, monomer, position.freq[monomer]];
  };

  /** Helper function for rendering */
  protected _nullSequence(fillerResidue = 'X'): string {
    if (!this.skipEmptySequences)
      return new Array(this.Length).fill(fillerResidue).join('');

    return '';
  }

  /** Helper function for remove empty positions */
  // TODO: use this function in from core
  protected removeWhere(array: Array<any>, predicate: (T: any) => boolean): Array<any> {
    const length = array.length;
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
      this.removeWhere(this.positions, (item) => item?.freq['-']?.count === item.rowCount);
    }
  }

  protected _calculate(r: number) {
    if (!this.host || !this.seqCol || !this.dataFrame)
      return;
    this.unitsHandler = new UnitsHandler(this.seqCol);

    this.calcSize();

    const posCount: number = this.startPosition <= this.endPosition ? this.endPosition - this.startPosition + 1 : 0;
    this.positions = new Array(posCount);
    for (let jPos = 0; jPos < this.Length; jPos++) {
      const posName: string = this.positionNames[this.startPosition + jPos];
      this.positions[jPos] = new PositionInfo(this.startPosition + jPos, posName);
    }

    // 2022-05-05 askalkin instructed to show WebLogo based on filter (not selection)
    const indices = this.filter.getSelectedIndexes();
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

      const alphabetSizeLog = Math.log2(alphabetSize);
      const maxHeight = (this.positionHeight == PositionHeight.Entropy) ?
        (absoluteMaxHeight * (alphabetSizeLog - (this.positions[jPos].sumForHeightCalc)) / alphabetSizeLog) :
        absoluteMaxHeight;

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

  // -- Handle events --

  private sliderOnValuesChanged(value: any): void {
    if ((this.host == null)) return;

    try {
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
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('Bio: WebLogoViewer.sliderOnValuesChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private dataFrameFilterOnChanged(value: any): void {
    console.debug('Bio: WebLogoViewer.dataFrameFilterChanged()');
    try {
      this.updatePositions();
      this.render();
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('Bio: WebLogoViewer.dataFrameFilterOnChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private dataFrameSelectionOnChanged(value: any): void {
    console.debug('Bio: WebLogoViewer.dataFrameSelectionOnChanged()');
    try {
      this.render();
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('Bio: WebLogoViewer.dataFrameSelectionOnChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private canvasOnMouseMove(e: MouseEvent) {
    try {
      const args = e as MouseEvent;

      const dpr: number = window.devicePixelRatio;
      const cursorP: DG.Point = this.canvas.getCursorPosition(args, dpr);
      const [jPos, monomer] = this.getMonomer(cursorP);
      // if (jPos != undefined && monomer == undefined) {
      //   const preEl = ui.element('pre');
      //   preEl.innerHTML = jPos < this.positions.length ?
      //     JSON.stringify(this.positions[jPos].freq, undefined, 2) : 'NO jPos';
      //   const tooltipEl = ui.div([ui.div(`pos: ${jPos}`), preEl]);
      //   ui.tooltip.show(tooltipEl, args.x + 16, args.y + 16);
      // } else
      if (this.dataFrame && this.seqCol && this.splitter && monomer) {
        const atPI: PositionInfo = this.positions[jPos];
        const monomerAtPosSeqCount = countForMonomerAtPosition(
          this.dataFrame, this.seqCol, this.filter, this.splitter, monomer, atPI);

        const tooltipEl = ui.div([
          // ui.div(`pos ${jPos}`),
          ui.div(`${monomer}`),
          ui.div(`${monomerAtPosSeqCount} rows`)]);
        ui.tooltip.show(tooltipEl, args.x + 16, args.y + 16);
      } else {
        ui.tooltip.hide();
      }
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('Bio: WebLogoViewer.canvasOnMouseMove() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private canvasOnMouseDown(e: MouseEvent): void {
    try {
      const args = e as MouseEvent;
      const r: number = window.devicePixelRatio;
      const [jPos, monomer] = this.getMonomer(this.canvas.getCursorPosition(args, r));

      // prevents deselect all rows if we miss monomer bounds
      if (this.dataFrame && this.seqCol && this.splitter && monomer) {
        const atPI: PositionInfo = this.positions[jPos];

        // this.dataFrame.selection.init((rowI: number) => {
        //   return checkSeqForMonomerAtPos(
        //     this.dataFrame, this.seqCol!, this.filter, rowI, this.splitter!, monomer, atPI);
        // });
        // Calculate a new BitSet object for selection to prevent interfering with existing
        const selBS: DG.BitSet = DG.BitSet.create(this.dataFrame.selection.length, (rowI: number) => {
          return checkSeqForMonomerAtPos(
            this.dataFrame, this.seqCol!, this.filter, rowI, this.splitter!, monomer, atPI);
        });
        this.dataFrame.selection.init((i) => selBS.get(i));
      }
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('Bio: WebLogoViewer.canvasOnMouseDown() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private canvasOnWheel(e: WheelEvent) {
    try {
      if (!this.visibleSlider)
        return;
      const countOfScrollPositions = (e.deltaY / 100) * Math.max(Math.floor((this.countOfRenderPositions) / 2), 1);
      this.slider.scrollBy(this.slider.min + countOfScrollPositions);
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('Bio: WebLogoViewer.canvasOnWheel() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }
}

export function checkSeqForMonomerAtPos(
  df: DG.DataFrame, seqCol: DG.Column, filter: DG.BitSet, rowI: number,
  splitter: SplitterFunc, monomer: string, at: PositionInfo
): boolean {
  // if (!filter.get(rowI)) return false;
  // TODO: Use BitSet.get(idx)
  if (!filter.getSelectedIndexes().includes(rowI)) return false;

  const seq = seqCol!.get(rowI);
  const seqM = seq ? splitter!(seq)[at.pos] : null;
  return ((seqM === monomer) || (seqM === '' && monomer === '-'));
}

export function countForMonomerAtPosition(
  df: DG.DataFrame, seqCol: DG.Column, filter: DG.BitSet,
  splitter: SplitterFunc, monomer: string, at: PositionInfo
): number {
  const posMList: (string | null)[] = wu.count(0).take(df.rowCount)
    .filter((rowI) => filter.get(rowI))
    .map((rowI) => {
      const seq: string | null = seqCol!.get(rowI);
      const seqMList: string[] = seq ? splitter!(seq) : [];
      const seqMPos: number = at.pos;
      const seqM: string | null = seqMPos < seqMList.length ? seqMList[seqMPos] : null;
      return seqM;
    }).toArray();
  // wu.count().take(this.dataFrame.rowCount).filter(function(iRow) {
  //   return correctMonomerFilter(iRow, monomer, jPos);
  // }).reduce<number>((count, iRow) => count + 1, 0);
  const monomerAtPosRowCount = posMList.filter((m) => m == monomer).reduce((count, m) => count + 1, 0);
  return monomerAtPosRowCount;
}
