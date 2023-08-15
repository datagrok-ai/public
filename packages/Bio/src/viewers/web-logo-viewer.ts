import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {fromEvent, Observable, Subject, Unsubscribable} from 'rxjs';

import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {
  monomerToShort,
  pickUpPalette,
  pickUpSeqCol,
  TAGS as bioTAGS,
} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  FilterSources,
  HorizontalAlignments,
  IWebLogoViewer,
  PositionHeight,
  PositionMarginStates,
  positionSeparator,
  TAGS as wlTAGS,
  VerticalAlignments,
  WebLogoProps,
  WebLogoPropsDefault,
} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {intToHtmlA} from '@datagrok-libraries/utils/src/color';

import {_package} from '../package';

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
  public count: number;

  /** Remember screen coords rect */
  public bounds?: DG.Rect;

  constructor(count: number = 0, bounds?: DG.Rect) {
    this.count = count;
    this.bounds = bounds;
  }
}

export class PositionInfo {
  /** Position in sequence */
  public readonly pos: number;

  /** Position name from column tag*/
  public readonly name: string;

  private readonly _label: string | undefined;
  public get label(): string { return !!this._label ? this._label : this.name; }

  private readonly _freqs: { [m: string]: PositionMonomerInfo };

  rowCount: number;
  sumForHeightCalc: number;

  /** freq = {}, rowCount = 0
   * @param {number} pos Position in sequence
   * @param {string} name Name of position ('111A', '111.1', etc)
   * @param {string[]} freqs frequency of monomers in position
   * @param {number} rowCount Count of elements in column
   * @param {number} sumForHeightCalc Sum of all monomer counts for height calculation
   */
  constructor(pos: number, name: string, freqs?: { [m: string]: PositionMonomerInfo },
    options?: { rowCount?: number, sumForHeightCalc?: number, label?: string }
  ) {
    this.pos = pos;
    this.name = name;
    this._freqs = freqs ?? {};

    if (options?.rowCount) this.rowCount = options.rowCount;
    if (options?.sumForHeightCalc) this.sumForHeightCalc = options.sumForHeightCalc;
    if (options?.label) this._label = options.label;
  }

  public getMonomers(): string[] {
    return Object.keys(this._freqs);
  }

  public hasMonomer(m: string): boolean {
    return m in this._freqs;
  }

  /** Creates empty PositionMonomerInfo for {@link m} key if missed in {@link _freqs}. */
  public getFreq(m: string): PositionMonomerInfo {
    let res: PositionMonomerInfo = this._freqs[m];
    if (!res) res = this._freqs[m] = new PositionMonomerInfo();
    return res;
  }

  getMonomerAt(calculatedX: number, y: number): string | undefined {
    const findRes = Object.entries(this._freqs)
      .find(([m, pmInfo]) => {
        return pmInfo.bounds!.contains(calculatedX, y);
      });
    return !!findRes ? findRes[0] : undefined;
  }

  calcHeights(heightMode: PositionHeight): void {
    /*
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
    /**/

    this.rowCount = 0;
    for (const [m, pmInfo] of Object.entries(this._freqs))
      this.rowCount += pmInfo.count;

    this.sumForHeightCalc = 0;
    if (heightMode === PositionHeight.Entropy) {
      for (const [m, pmInfo] of Object.entries(this._freqs)) {
        const pn = pmInfo.count / this.rowCount;
        this.sumForHeightCalc += -pn * Math.log2(pn);
      }
    } else if (heightMode === PositionHeight.full) {
      for (const [m, pmInfo] of Object.entries(this._freqs)) {
        const pn = pmInfo.count / this.rowCount;
        this.sumForHeightCalc += pn;
      }
    }
    const k = 42;
  }

  calcScreen(
    jPos: number, absoluteMaxHeight: number, heightMode: PositionHeight, alphabetSizeLog: number,
    positionWidthWithMargin: number, positionWidth: number, r: number, axisHeight: number
  ): void {
    const dpr = window.devicePixelRatio;

    const maxHeight = (heightMode == PositionHeight.Entropy) ?
      (absoluteMaxHeight * (alphabetSizeLog - (this.sumForHeightCalc)) / alphabetSizeLog) :
      absoluteMaxHeight;
    let y: number = axisHeight * r + (absoluteMaxHeight - maxHeight - 1);

    const entries = Object.entries(this._freqs)
      .sort((a, b) => {
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
      const h: number = maxHeight * pmInfo.count / this.rowCount;

      pmInfo.bounds = new DG.Rect(jPos * dpr * positionWidthWithMargin, y, positionWidth * dpr, h);
      y += h;
    }
  }

  render(g: CanvasRenderingContext2D,
    fontStyle: string, uppercaseLetterAscent: number, uppercaseLetterHeight: number,
    positionWidthWithMargin: number, firstVisibleIndex: number, cp: SeqPalette
  ) {
    for (const [monomer, pmInfo] of Object.entries(this._freqs)) {
      if (monomer !== '-') {
        const monomerTxt = monomerToShort(monomer, 5);
        const b = pmInfo.bounds!;
        const left = b.left - positionWidthWithMargin * firstVisibleIndex;

        g.resetTransform();
        g.strokeStyle = 'lightgray';
        g.lineWidth = 1;
        g.rect(left, b.top, b.width, b.height);
        g.fillStyle = cp.get(monomer) ?? cp.get('other');
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
  showPositionLabels = 'showPositionLabels',
  positionMarginState = 'positionMarginState',
  positionMargin = 'positionMargin',

  // -- Behavior --
  filterSource = 'filterSource',
}

const defaults: WebLogoProps = WebLogoPropsDefault;

enum RecalcLevel {
  None = 0,
  Layout = 1,
  Freqs = 2,
}

const POSITION_LABELS_HEIGHT: number = 12;

export class WebLogoViewer extends DG.JsViewer implements IWebLogoViewer {
  public static residuesSet = 'nucleotides';
  private static viewerCount: number = -1;

  private viewed: boolean = false;

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

  private seqCol: DG.Column<string> | null = null;
  // private maxLength: number = 100;
  private positions: PositionInfo[] = [];

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
  public showPositionLabels: boolean;
  public positionMarginState: PositionMarginStates;
  public positionMargin: number = 0;
  public startPositionName: string | null;
  public endPositionName: string | null;
  public fixWidth: boolean;
  public verticalAlignment: VerticalAlignments;
  public horizontalAlignment: HorizontalAlignments;
  public fitArea: boolean;
  public shrinkEmptyTail: boolean;
  public positionHeight: string;

  // -- Behavior --
  public filterSource: FilterSources;

  private positionNames: string[] = [];
  private positionLabels: string[] | undefined = undefined;
  private startPosition: number = -1;
  private endPosition: number = -1;

  private error: Error | null = null;

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
  public get Length(): number {
    if (this.skipEmptyPositions)
      return this.positions.length;

    return this.startPosition <= this.endPosition ? this.endPosition - this.startPosition + 1 : 0;
  }

  /** Calculate new position data basic on {@link positionMarginState} and {@link positionMargin} */
  private get positionWidthWithMargin() {
    return this._positionWidth + this.positionMarginValue;
  }

  private get positionMarginValue() {
    if (this.positionMarginState === PositionMarginStates.AUTO &&
      this.unitsHandler?.getAlphabetIsMultichar() === true
    ) return this.positionMargin;

    if (this.positionMarginState === PositionMarginStates.ON)
      return this.positionMargin;

    return 0;
  }

  /** Count of position rendered for calculations countOfRenderPositions */
  private get countOfRenderPositions() {
    if (this.host == null)
      return 0;

    const r = window.devicePixelRatio;
    if (r > 1)
      return this.canvasWidthWithRatio / this.positionWidthWithMargin;
    else
      return this.canvas.width / (this.positionWidthWithMargin * r);
  }

  private get canvasWidthWithRatio() {
    return this.canvas.width * window.devicePixelRatio;
  }


  /** Position of start rendering */
  private get firstVisibleIndex(): number {
    return (this.visibleSlider) ? Math.floor(this.slider.min) : 0;
  }

  constructor() {
    super();

    this.viewerId = WebLogoViewer.viewerCount;
    WebLogoViewer.viewerCount += 1;

    this.textBaseline = 'top';
    this.unitsHandler = null;

    // -- Data --
    this.sequenceColumnName = this.string(PROPS.sequenceColumnName, defaults.sequenceColumnName,
      {category: PROPS_CATS.DATA});
    this.startPositionName = this.string(PROPS.startPositionName, defaults.startPositionName,
      {category: PROPS_CATS.DATA});
    this.endPositionName = this.string(PROPS.endPositionName, defaults.endPositionName,
      {category: PROPS_CATS.DATA});
    this.skipEmptySequences = this.bool(PROPS.skipEmptySequences, defaults.skipEmptySequences,
      {category: PROPS_CATS.DATA});
    this.skipEmptyPositions = this.bool(PROPS.skipEmptyPositions, defaults.skipEmptyPositions,
      {category: PROPS_CATS.DATA});
    this.shrinkEmptyTail = this.bool(PROPS.shrinkEmptyTail, defaults.shrinkEmptyTail,
      {category: PROPS_CATS.DATA});

    // -- Style --
    this.backgroundColor = this.int(PROPS.backgroundColor, defaults.backgroundColor,
      {category: PROPS_CATS.STYLE});
    this.positionHeight = this.string(PROPS.positionHeight, defaults.positionHeight,
      {category: PROPS_CATS.STYLE, choices: Object.values(PositionHeight)});
    this._positionWidth = this.positionWidth = this.float(PROPS.positionWidth, defaults.positionWidth,
      {category: PROPS_CATS.STYLE/* editor: 'slider', min: 4, max: 64, postfix: 'px' */});

    // -- Layout --
    this.verticalAlignment = this.string(PROPS.verticalAlignment, defaults.verticalAlignment,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(VerticalAlignments)}) as VerticalAlignments;
    this.horizontalAlignment = this.string(PROPS.horizontalAlignment, defaults.horizontalAlignment,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(HorizontalAlignments)}) as HorizontalAlignments;
    this.fixWidth = this.bool(PROPS.fixWidth, defaults.fixWidth,
      {category: PROPS_CATS.LAYOUT});
    this.fitArea = this.bool(PROPS.fitArea, defaults.fitArea,
      {category: PROPS_CATS.LAYOUT});
    this.minHeight = this.float(PROPS.minHeight, defaults.minHeight,
      {category: PROPS_CATS.LAYOUT/*, editor: 'slider', min: 25, max: 250, postfix: 'px'*/});
    this.maxHeight = this.float(PROPS.maxHeight, defaults.maxHeight,
      {category: PROPS_CATS.LAYOUT/*, editor: 'slider', min: 25, max: 500, postfix: 'px'*/});
    this.showPositionLabels = this.bool(PROPS.showPositionLabels, defaults.showPositionLabels,
      {category: PROPS_CATS.LAYOUT});
    this.positionMarginState = this.string(PROPS.positionMarginState, defaults.positionMarginState,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(PositionMarginStates)}) as PositionMarginStates;
    let defaultValueForPositionMargin = 0;
    if (this.positionMarginState === 'auto') defaultValueForPositionMargin = 4;
    this.positionMargin = this.int(PROPS.positionMargin, defaultValueForPositionMargin,
      {category: PROPS_CATS.LAYOUT, min: 0, max: 16});

    // -- Behavior --
    this.filterSource = this.string(PROPS.filterSource, defaults.filterSource,
      {category: PROPS_CATS.BEHAVIOR, choices: Object.values(FilterSources)}) as FilterSources;

    const style: DG.SliderOptions = {style: 'barbell'};
    this.slider = ui.rangeSlider(0, 100, 0, 20, false, style);
    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';
  }

  // -- Data --

  setData(): void {
    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }

    this.updateSeqCol();

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private viewSubs: Unsubscribable[] = [];

  private destroyView() {
    for (const sub of this.viewSubs) sub.unsubscribe();
    this.viewSubs = [];

    const dataFrameTxt = `${this.dataFrame ? 'data' : 'null'}`;
    _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.destroyView( dataFrame = ${dataFrameTxt} ) start`);
    super.detach();

    this.viewSubs.forEach((sub) => sub.unsubscribe());
    this.host!.remove();
    this.msgHost = undefined;
    this.host = undefined;

    _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.destroyView() end`);
  }

  private buildView() {
    const dataFrameTxt: string = this.dataFrame ? 'data' : 'null';
    _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.buildView( dataFrame = ${dataFrameTxt} ) start`);

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

    this.host = ui.div([this.msgHost, this.canvas]);

    this.host.style.justifyContent = 'center';
    this.host.style.alignItems = 'center';
    this.host.style.position = 'relative';
    this.host.style.setProperty('overflow', 'hidden', 'important');

    this.root.append(this.host);
    this.root.append(this.slider.root);

    this.updateSlider();
    this.render(RecalcLevel.Freqs, 'init');

    if (!!this.error) {
      this.msgHost!.innerText = this.error.message;
      ui.tooltip.bind(this.msgHost!, this.error.stack);
      this.msgHost!.style.setProperty('display', null);
    }

    if (this.dataFrame) {
      this.viewSubs.push(this.dataFrame.filter.onChanged.subscribe(this.dataFrameFilterOnChanged.bind(this)));
      this.viewSubs.push(this.dataFrame.selection.onChanged.subscribe(this.dataFrameSelectionOnChanged.bind(this)));
    }
    this.viewSubs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
    this.viewSubs.push(this.slider.onValuesChanged.subscribe(this.sliderOnValuesChanged.bind(this)));
    this.viewSubs.push(
      fromEvent<MouseEvent>(this.canvas, 'mousemove').subscribe(this.canvasOnMouseMove.bind(this)));
    this.viewSubs.push(
      fromEvent<MouseEvent>(this.canvas, 'mousedown').subscribe(this.canvasOnMouseDown.bind(this)));
    this.viewSubs.push(
      fromEvent<WheelEvent>(this.canvas, 'wheel').subscribe(this.canvasOnWheel.bind(this)));

    _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.buildView() end`);
  }

  /** Handler of changing size WebLogo */
  private rootOnSizeChanged(): void {
    _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.rootOnSizeChanged(), start `);

    this.updateSlider();
    this.render(RecalcLevel.Layout, 'rootOnSizeChanged');
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
        try {
          const units: string = this.seqCol!.getTag(DG.TAGS.UNITS);
          const separator: string = this.seqCol!.getTag(bioTAGS.separator);
          this.unitsHandler = UnitsHandler.getOrCreate(this.seqCol);

          this.updatePositions();
          this.cp = pickUpPalette(this.seqCol);
          this.error = null;
        } catch (err: any) {
          this.seqCol = null;
          this.error = err instanceof Error ? err : new Error(err.toString());
          throw err;
        }
        if (!this.seqCol) {
          this.unitsHandler = null;
          this.positionNames = [];
          this.positionLabels = [];
          this.startPosition = -1;
          this.endPosition = -1;
          this.cp = null;
        }
      }
    }
  }

  /** Updates {@link positionNames} and calculates {@link startPosition} and {@link endPosition}.
   */
  private updatePositions(): void {
    if (!this.seqCol)
      return;

    const dfFilter = this.filterSource === FilterSources.Filtered ? this.dataFrame.filter :
      this.dataFrame.selection;
    const maxLength = wu.enumerate(this.unitsHandler!.splitted).map(([mList, rowI]) => {
      return dfFilter.get(rowI) && !!mList ? mList.length : 0;
    }).reduce((max, l) => Math.max(max, l), 0);

    /** positionNames and positionLabel can be set up through the column's tags only */
    const positionNamesTxt = this.seqCol.getTag(wlTAGS.positionNames);
    const positionLabelsTxt = this.seqCol.getTag(wlTAGS.positionLabels);
    this.positionNames = !!positionNamesTxt ? positionNamesTxt.split(positionSeparator).map((v) => v.trim()) :
      [...Array(maxLength).keys()].map((jPos) => `${jPos + 1}`)/* fallback if tag is not provided */;
    this.positionLabels = !!positionLabelsTxt ? positionLabelsTxt.split(positionSeparator).map((v) => v.trim()) :
      undefined;

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

    if (((minScale == this.xScale) || (minScale <= 1)) && (this.fitArea))
      showSliderWithFitArea = false;

    return ((this.fixWidth || Math.ceil(this.canvas.width / this.positionWidthWithMargin) >= this.Length) ||
      (showSliderWithFitArea));
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
    if (this.checkIsHideSlider())
      this.setSliderVisibility(false);
    else
      this.setSliderVisibility(true);

    if ((this.slider != null) && (this.canvas != null)) {
      const diffEndScrollAndSliderMin = Math.max(0,
        Math.floor(this.slider.min + this.canvas.width / this.positionWidthWithMargin) - this.Length);
      let newMin = Math.floor(this.slider.min - diffEndScrollAndSliderMin);
      let newMax = Math.floor(this.slider.min - diffEndScrollAndSliderMin) +
        Math.floor(this.canvas.width / this.positionWidthWithMargin);
      if (this.checkIsHideSlider()) {
        newMin = 0;
        newMax = Math.max(newMin, this.Length - 1);
      }
      this.turnOfResizeForOneSetValue = true;
      this.slider.setValues(0, this.Length,
        newMin, newMax);
    }
  }

  /** Handler of property change events.
   * @param {DG.Property} property - property which was changed.
   */
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
      case PROPS.showPositionLabels:
      case PROPS.shrinkEmptyTail:
      case PROPS.skipEmptyPositions:
        this.updatePositions();
        break;
    }

    switch (property.name) {
      case PROPS.fixWidth:
      case PROPS.fitArea:
      case PROPS.showPositionLabels:
      case PROPS.positionWidth:
      case PROPS.positionMargin:
        this.render(RecalcLevel.Layout, 'onPropertyChanged');
        break;

      default:
        this.render(RecalcLevel.Freqs, 'onPropertyChanged');
        break;
    }
  }

  /** Add filter handlers when table is a attached  */
  public override onTableAttached() {
    _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.onTableAttached(), `);

    // -- Props editors --

    super.onTableAttached();
    this.setData();
  }

  /** Remove all handlers when table is a detach  */
  public override async detach() {
    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }
  }

  private _onSizeChanged: Subject<void> = new Subject<void>();

  public get onSizeChanged(): Observable<void> { return this._onSizeChanged; }

  // -- Routines --

  getMonomer(p: DG.Point, dpr: number): [number, string | null, PositionMonomerInfo | null] {
    const calculatedX = p.x + this.firstVisibleIndex * this.positionWidthWithMargin * dpr;
    const jPos = Math.floor(p.x / (this.positionWidthWithMargin * dpr) + this.firstVisibleIndex);
    const position: PositionInfo = this.positions[jPos];

    if (position === undefined)
      return [jPos, null, null];

    const monomer: string | undefined = position.getMonomerAt(calculatedX, p.y);
    if (monomer === undefined)
      return [jPos, null, null];

    return [jPos, monomer, position.getFreq(monomer)];
  };

  /** Helper function for rendering
   * @param {boolean} fillerResidue
   * @return {string} - string with null sequence
   */
  protected _nullSequence(fillerResidue = 'X'): string {
    if (!this.skipEmptySequences)
      return new Array(this.Length).fill(fillerResidue).join('');

    return '';
  }

  /** Function for removing empty positions */
  protected _removeEmptyPositions() {
    if (this.skipEmptyPositions) {
      this.positions = wu(this.positions)
        .filter((pi) => !pi.hasMonomer('-') || pi.getFreq('-').count !== pi.rowCount).toArray();
    }
  }

  private renderRequested: boolean = false;
  /** default value of RecalcLevel.Freqs is for recalc from the scratch at the beginning */
  private recalcLevelRequested: RecalcLevel = RecalcLevel.Freqs;

  /** Renders requested repeatedly will be performed once on window.requestAnimationFrame() */
  render(recalcLevel: RecalcLevel, reason: string): void {
    _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>` +
      `.render( recalcLevel=${recalcLevel}, reason='${reason}' )`);

    /** Calculate freqs of monomers */
    const calculateFreqsInt = (): void => {
      _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.render.calculateFreqsInt(), start `);
      if (!this.host || !this.seqCol || !this.dataFrame) return;

      const length: number = this.startPosition <= this.endPosition ? this.endPosition - this.startPosition + 1 : 0;
      this.unitsHandler = UnitsHandler.getOrCreate(this.seqCol);
      const posCount: number = this.startPosition <= this.endPosition ? this.endPosition - this.startPosition + 1 : 0;
      this.positions = new Array(posCount);
      for (let jPos = 0; jPos < length; jPos++) {
        const posName: string = this.positionNames[this.startPosition + jPos];
        const posLabel: string | undefined = this.positionLabels ?
          this.positionLabels[this.startPosition + jPos] : undefined;
        this.positions[jPos] = new PositionInfo(this.startPosition + jPos, posName, {}, {label: posLabel});
      }

      // 2022-05-05 askalkin instructed to show WebLogo based on filter (not selection)
      const dfFilter =
        this.filterSource === FilterSources.Filtered ? this.dataFrame.filter : this.dataFrame.selection;
      const dfRowCount = this.dataFrame.rowCount;
      const splitted = this.unitsHandler.splitted;
      for (let rowI = 0; rowI < dfRowCount; ++rowI) {
        if (dfFilter.get(rowI)) {
          const seqMList: string[] = splitted[rowI];
          for (let jPos = 0; jPos < length; ++jPos) {
            const m: string = seqMList[this.startPosition + jPos] || '-';
            const pmInfo = this.positions[jPos].getFreq(m);
            pmInfo.count++;
          }
        }
      }

      //#region Polish freq counts
      for (let jPos = 0; jPos < length; jPos++) {
        // delete this.positions[jPos].freq['-'];
        this.positions[jPos].calcHeights(this.positionHeight as PositionHeight);
      }
      //#endregion
      this._removeEmptyPositions();
    };

    /** Calculate layout of monomers on screen (canvas) based on freqs, required to handle mouse events */
    const calculateLayoutInt = (dpr: number, positionLabelsHeight: number): void => {
      _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.render.calculateLayoutInt(), start `);

      const length = this.positions.length;
      this.calcSize();
      const absoluteMaxHeight = this.canvas.height - positionLabelsHeight * dpr;
      const alphabetSize = this.getAlphabetSize();
      if ((this.positionHeight == PositionHeight.Entropy) && (alphabetSize == null))
        grok.shell.error('WebLogo: alphabet is undefined.');
      const alphabetSizeLog = Math.log2(alphabetSize);

      for (let jPos = 0; jPos < length; jPos++) {
        this.positions[jPos].calcScreen(jPos, absoluteMaxHeight, this.positionHeight as PositionHeight,
          alphabetSizeLog, this.positionWidthWithMargin, this._positionWidth, dpr, positionLabelsHeight);
      }
    };

    /** Render WebLogo sensitive to changes in params of rendering
     *@param {boolean} recalcFreqs - indicates that need to recalculate data for rendering
     */
    const renderInt = (recalcLevel: RecalcLevel) => {
      _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.render.renderInt( recalcLevel=${recalcLevel} ), ` +
        `start `);
      if (this.msgHost) {
        if (this.seqCol && !this.cp) {
          this.msgHost!.innerText = `Unknown palette (column semType: '${this.seqCol.semType}').`;
          this.msgHost!.style.display = '';
        } else {
          this.msgHost!.style.display = 'none';
        }
      }

      if (!this.seqCol || !this.dataFrame || !this.cp || this.startPosition === -1 ||
        this.endPosition === -1 || this.host == null || this.slider == null)
        return;

      const g = this.canvas.getContext('2d');
      if (!g) return;

      this.slider.root.style.width = `${this.host.clientWidth}px`;

      const dpr: number = window.devicePixelRatio;
      /** 0 is for no position labels */
      const positionLabelsHeight = this.showPositionLabels ? POSITION_LABELS_HEIGHT : 0;
      if (recalcLevel >= RecalcLevel.Freqs) calculateFreqsInt();
      if (recalcLevel >= RecalcLevel.Layout) calculateLayoutInt(window.devicePixelRatio, positionLabelsHeight);

      const length: number = this.Length;
      g.resetTransform();
      g.fillStyle = intToHtmlA(this.backgroundColor);
      g.fillRect(0, 0, this.canvas.width, this.canvas.height);
      g.textBaseline = this.textBaseline;

      const maxCountOfRowsRendered = this.countOfRenderPositions + 1;
      const firstVisibleIndex = (this.visibleSlider) ? Math.floor(this.slider.min) : 0;
      const lastVisibleIndex = Math.min(length, firstVisibleIndex + maxCountOfRowsRendered);

      //#region Plot positionNames
      const positionFontSize = 10 * dpr;
      g.resetTransform();
      g.fillStyle = 'black';
      g.textAlign = 'center';
      g.font = `${positionFontSize.toFixed(1)}px Roboto, Roboto Local, sans-serif`;
      const posNameMaxWidth = Math.max(...this.positions.map((pos) => g.measureText(pos.name).width));
      const hScale = posNameMaxWidth < (this._positionWidth * dpr - 2) ? 1 :
        (this._positionWidth * dpr - 2) / posNameMaxWidth;

      if (positionLabelsHeight > 0) {
        renderPositionLabels(g, dpr, hScale, this.positionWidthWithMargin, this._positionWidth,
          this.positions, this.firstVisibleIndex, lastVisibleIndex);
      }
      //#endregion Plot positionNames
      const fontStyle = '16px Roboto, Roboto Local, sans-serif';
      // Hacks to scale uppercase characters to target rectangle
      const uppercaseLetterAscent = 0.25;
      const uppercaseLetterHeight = 12.2;
      for (let jPos = this.firstVisibleIndex; jPos < lastVisibleIndex; jPos++) {
        this.positions[jPos].render(g, fontStyle, uppercaseLetterAscent, uppercaseLetterHeight,
          this.positionWidthWithMargin, firstVisibleIndex, this.cp);
      }

      _package.logger.debug(`Bio: WebLogoViewer<${this.viewerId}>.render.renderInt( recalcLevel=${recalcLevel} ), end`);
    };

    this.recalcLevelRequested = Math.max(this.recalcLevelRequested, recalcLevel);
    if (!this.renderRequested) {
      this.renderRequested = true;
      // requestAnimationFrame callback will be executed after this.render()
      switch (this.recalcLevelRequested) {
        case RecalcLevel.Freqs:
          /* Avoiding [Violation] 'requestAnimationFrame' handler took too much */
          window.setTimeout(() => {
            renderInt(this.recalcLevelRequested);
            this.recalcLevelRequested = RecalcLevel.None;
            this.renderRequested = false;
          }, 0 /* next event cycle */);
          break;

        case RecalcLevel.Layout:
        case RecalcLevel.None:
          window.requestAnimationFrame((time: number) => {
            renderInt(this.recalcLevelRequested);
            this.recalcLevelRequested = RecalcLevel.None;
            this.renderRequested = false;
          });
          break;
      }
    }
  }

  private _lastWidth: number;
  private _lastHeight: number;

  /** Calculate canvas size an positionWidth and updates properties */
  private calcSize() {
    if (!this.host) return;

    const dpr: number = window.devicePixelRatio;

    let width: number = this.widthArea;
    let height = this.heightArea;

    if ((this.fitArea) && (!this.visibleSlider)) {
      const scale = Math.max(1, Math.min(this.xScale, this.yScale));
      width = width * scale;
      height = height * scale;
      this._positionWidth = this.positionWidth * scale;
    }

    width = this.Length * this.positionWidthWithMargin;

    this.canvas.width = this.root.clientWidth * dpr;
    this.canvas.style.width = `${this.root.clientWidth}px`;

    // const canvasHeight: number = width > this.root.clientWidth ? height - 8 : height;
    this.host.style.setProperty('height', `${height}px`);
    const canvasHeight: number = this.host.clientHeight;
    this.canvas.height = canvasHeight * dpr;

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
        case HorizontalAlignments.LEFT:
          hostLeftMargin = 0;
          break;
        case HorizontalAlignments.CENTER:
          hostLeftMargin = Math.max(0, (this.root.clientWidth - width) / 2);
          break;
        case HorizontalAlignments.RIGHT:
          hostLeftMargin = Math.max(0, this.root.clientWidth - width);
          break;
      }
      this.host.style.setProperty('margin-top', `${hostTopMargin}px`, 'important');
      this.host.style.setProperty('margin-left', `${hostLeftMargin}px`, 'important');
      if (this.slider != null)
        this.slider.root.style.setProperty('margin-top', `${hostTopMargin + canvasHeight}px`, 'important');


      if (this.root.clientHeight <= height) {
        this.host.style.setProperty('height', `${this.root.clientHeight}px`);
        this.host.style.setProperty('overflow-y', null);
      } else {
        this.host.style.setProperty('overflow-y', 'hidden', 'important');
      }
    }

    if (this._lastWidth !== this.root.clientWidth && this._lastHeight !== this.root.clientHeight) {
      this._lastWidth = this.root.clientWidth;
      this._lastHeight = this.root.clientHeight;
      this._onSizeChanged.next();
    }
  }

  public getAlphabetSize(): number {
    return this.unitsHandler?.getAlphabetSize() ?? 0;
  }

  // -- Handle events --

  private sliderOnValuesChanged(_value: any): void {
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
      this.render(RecalcLevel.Layout, 'sliderOnValuesChanged');
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      _package.logger.error('Bio: WebLogoViewer<${this.viewerId}>.sliderOnValuesChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private dataFrameFilterOnChanged(_value: any): void {
    _package.logger.debug('Bio: WebLogoViewer<${this.viewerId}>.dataFrameFilterChanged()');
    try {
      this.updatePositions();
      if (this.filterSource === FilterSources.Filtered)
        this.render(RecalcLevel.Freqs, 'dataFrameFilterOnChanged');
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      _package.logger.error('Bio: WebLogoViewer<${this.viewerId}>.dataFrameFilterOnChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private dataFrameSelectionOnChanged(_value: any): void {
    _package.logger.debug('Bio: WebLogoViewer<${this.viewerId}>.dataFrameSelectionOnChanged()');
    try {
      if (this.filterSource === FilterSources.Selected)
        this.render(RecalcLevel.Freqs, 'dataFrameSelectionOnChanged');
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      _package.logger.error('Bio: WebLogoViewer<${this.viewerId}>.dataFrameSelectionOnChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private canvasOnMouseMove(e: MouseEvent) {
    const dpr = window.devicePixelRatio;
    try {
      const args = e as MouseEvent;

      const cursorP: DG.Point = this.canvas.getCursorPosition(args, dpr);
      const [jPos, monomer] = this.getMonomer(cursorP, dpr);
      // if (jPos != undefined && monomer == undefined) {
      //   const preEl = ui.element('pre');
      //   preEl.innerHTML = jPos < this.positions.length ?
      //     JSON.stringify(this.positions[jPos].freq, undefined, 2) : 'NO jPos';
      //   const tooltipEl = ui.div([ui.div(`pos: ${jPos}`), preEl]);
      //   ui.tooltip.show(tooltipEl, args.x + 16, args.y + 16);
      // } else
      if (this.dataFrame && this.seqCol && monomer) {
        const atPI: PositionInfo = this.positions[jPos];
        const monomerAtPosSeqCount = countForMonomerAtPosition(
          this.dataFrame, this.unitsHandler!, this.filter, monomer, atPI);

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
      _package.logger.error('Bio: WebLogoViewer<${this.viewerId}>.canvasOnMouseMove() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private canvasOnMouseDown(e: MouseEvent): void {
    try {
      const args = e as MouseEvent;
      const dpr: number = window.devicePixelRatio;
      const [jPos, monomer] = this.getMonomer(this.canvas.getCursorPosition(args, dpr), dpr);

      // prevents deselect all rows if we miss monomer bounds
      if (this.dataFrame && this.seqCol && this.unitsHandler && monomer) {
        const atPI: PositionInfo = this.positions[jPos];

        // Calculate a new BitSet object for selection to prevent interfering with existing
        const selBS: DG.BitSet = DG.BitSet.create(this.dataFrame.selection.length, (rowI: number) => {
          return checkSeqForMonomerAtPos(this.dataFrame, this.unitsHandler!, this.filter, rowI, monomer, atPI);
        });
        this.dataFrame.selection.init((i) => selBS.get(i));
      }
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      _package.logger.error('Bio: WebLogoViewer<${this.viewerId}>.canvasOnMouseDown() error:\n' + errMsg);
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
      _package.logger.error('Bio: WebLogoViewer<${this.viewerId}>.canvasOnWheel() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }
}

function renderPositionLabels(g: CanvasRenderingContext2D,
  dpr: number, hScale: number, positionWidthWithMargin: number, positionWidth: number,
  positions: PositionInfo[], firstVisibleIndex: number, lastVisibleIndex: number): void {
  for (let jPos = firstVisibleIndex; jPos < lastVisibleIndex; jPos++) {
    const pos: PositionInfo = positions[jPos];
    g.resetTransform();
    g.setTransform(
      hScale, 0, 0, 1,
      jPos * positionWidthWithMargin * dpr + positionWidth * dpr / 2 -
      positionWidthWithMargin * firstVisibleIndex, 0);
    g.fillText(pos.label, 0, 1);
  }
}

export function checkSeqForMonomerAtPos(
  df: DG.DataFrame, unitsHandler: UnitsHandler, filter: DG.BitSet, rowI: number, monomer: string, at: PositionInfo,
): boolean {
  const seqMList: string[] = unitsHandler.splitted[rowI];
  const seqM = at.pos < seqMList.length ? seqMList[at.pos] : null;
  return ((seqM === monomer) || (seqM === '' && monomer === '-'));
}

export function countForMonomerAtPosition(
  df: DG.DataFrame, uh: UnitsHandler, filter: DG.BitSet, monomer: string, at: PositionInfo
): number {
  let count = 0;
  let rowI = -1;
  while ((rowI = filter.findNext(rowI, true)) != -1) {
    const seqMList: string[] = uh.splitted[rowI];
    const seqMPos: number = at.pos;
    const seqM: string | null = seqMPos < seqMList.length ? seqMList[seqMPos] : null;
    if (seqM === monomer) count++;
  }
  return count;
}
