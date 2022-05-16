import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

import {
  Aminoacids, AminoacidsPalettes, Nucleotides, NucleotidesPalettes,
  PositionInfo, PositionMonomerInfo
} from '@datagrok-libraries/bio';


import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';


// Using color schemes from chem-palette

declare module 'datagrok-api/src/grid' {
  interface Rect {
    contains(x: number, y: number): boolean;
  }
}

declare global {
  interface HTMLCanvasElement {
    getCursorPosition(event: MouseEvent): DG.Point;
  }
}

HTMLCanvasElement.prototype.getCursorPosition = function(event: MouseEvent): DG.Point {
  const rect = this.getBoundingClientRect();
  return new DG.Point(event.clientX - rect.left, event.clientY - rect.top);
};

DG.Rect.prototype.contains = function(x: number, y: number): boolean {
  return this.left <= x && x <= this.right && this.top <= y && y <= this.bottom;
};

export class WebLogoViewer extends DG.JsViewer {
  public static residuesSet = 'nucleotides';

  // private readonly colorScheme: ColorScheme = ColorSchemes[NucleotidesWebLogo.residuesSet];
  protected cp: StringDictionary | null = null;

  // private readonly host: HTMLDivElement;
  private readonly canvas: HTMLCanvasElement;
  private readonly slider: DG.RangeSlider;
  private readonly textBaseline: CanvasTextBaseline;

  private minHeight: number;
  private maxHeight: number;
  private axisHeight: number = 12;

  private positionWidth: number = 16;

  private seqCol: DG.Column | null = null;
  // private maxLength: number = 100;
  private positions: PositionInfo[] = [];

  private rowsMasked: number = 0;
  private rowsNull: number = 0;

  // Viewer's properties (likely they should be public so that they can be set outside)
  public considerNullSequences: boolean;
  public sequenceColumnName: string;
  public startPositionName: string;
  public endPositionName: string;

  private positionNames: string[] = [];

  private startPosition: number = -1;

  private endPosition: number = -1;

  private get Length(): number { return this.endPosition - this.startPosition + 1; }

  constructor() {
    super();
    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';

    this.slider = ui.rangeSlider(0, 20, 2, 5);
    this.slider.root.style.width = '100%';
    this.slider.root.style.height = '12px';

    // this.host = ui.divV([/*this.slider,*/this.canvas]);

    // this.canvas.style.height = '100%';
    this.textBaseline = 'top';

    this.minHeight = this.float('minHeight', 50);
    this.maxHeight = this.float('maxHeight', 100);

    this.considerNullSequences = this.bool('considerNullSequences', false);
    this.sequenceColumnName = this.string('sequenceColumnName', null);
    this.startPositionName = this.string('startPositionName', null);
    this.endPositionName = this.string('endPositionName', null);

    const getMonomer = (p: DG.Point): [number, string | null, PositionMonomerInfo | null] => {
      const jPos = Math.floor(p.x / this.positionWidth);
      const position = this.positions[jPos];

      if (position === void 0)
        return [jPos, null, null];

      const monomer: string | undefined = Object.keys(position.freq)
        .find((m) => position.freq[m].bounds.contains(p.x, p.y));
      if (monomer === undefined)
        return [jPos, null, null];

      return [jPos, monomer, position.freq[monomer]];
    };

    this.canvas.onmouseover = (e: MouseEvent) => {

    };

    rxjs.fromEvent(this.canvas, 'mousemove').subscribe((e: Event) => {
      const args = e as MouseEvent;
      const [jPos, monomer] = getMonomer(this.canvas.getCursorPosition(args));

      if (this.dataFrame && this.seqCol && monomer) {
        ui.tooltip.showRowGroup(this.dataFrame, (iRow) => {
          const seq = this.seqCol!.get(iRow);
          const mSeq = seq ? seq[jPos] : null;
          return mSeq === monomer;
        }, args.x + 16, args.y + 16);
      } else {
        ui.tooltip.hide();
      }
    });

    rxjs.fromEvent(this.canvas, 'mousedown').subscribe((e: Event) => {
      const args = e as MouseEvent;
      const [jPos, monomer] = getMonomer(this.canvas.getCursorPosition(args));

      // prevents deselect all rows if we miss monomer bounds
      if (this.dataFrame && this.seqCol && monomer) {
        this.dataFrame.selection.init((iRow) => {
          const seq = this.seqCol!.get(iRow);
          const mSeq = seq ? seq[jPos] : null;
          return mSeq === monomer;
        });
      }
    });

    this.root.append(this.canvas);
    // this.root.appendChild(this.slider.root);
    // this.root.append(this.host);

    // ui.onSizeChanged(this.canvas).subscribe(this.canvasOnSizeChanged.bind(this));
    ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this));
  }

  // canvasOnSizeChanged(args: any) {
  //   this.canvas.width = this.canvas.clientWidth;
  //   this.canvas.height = this.canvas.clientHeight;
  //   console.debug(`WebLogo.onSizeChanged() width=${this.canvas.width}, height=${this.canvas.height}`);
  //   this.render(false);
  // }

  rootOnSizeChanged(args: any) {
    // this.canvas.width = this.root.clientWidth;
    // this.canvas.style.width = `${this.root.clientWidth}px`;
    let height = Math.min(this.maxHeight, Math.max(this.minHeight, this.root.clientHeight)) - 2;
    if (this.canvas.width > this.root.clientWidth) /* horizontal scroller is enabled */
      height -= 6; /* free some space for horizontal scroller */
    this.canvas.height = height;
    this.canvas.style.height = `${height}px`;

    this.root.style.width = `${this.canvas.width}px`;

    console.debug(`WebLogo.onRootSizeChanged() ` +
      `root.width=${this.root.clientWidth}, root.height=${this.root.clientHeight}, ` +
      `canvas.width=${this.canvas.width}, canvas.height=${this.canvas.height} .`);

    this.render(true);
  }

  /** Assigns {@link seqCol} and {@link cp} based on {@link sequenceColumnName} and calls {@link render}().
   */
  updateSeqCol(): void {
    if (this.dataFrame) {
      this.seqCol = this.dataFrame.col(this.sequenceColumnName);
      if (this.seqCol) {
        let maxLength = 0;
        for (const category of this.seqCol.categories)
          maxLength = Math.max(maxLength, category.length);

        // Get position names from data column tag 'positionNames'
        const positionNamesTxt = this.seqCol.getTag('positionNames');
        // Fallback if 'positionNames' tag is not provided
        this.positionNames = positionNamesTxt ? positionNamesTxt.split(', ').map((n) => n.trim()) :
          [...Array(maxLength).keys()].map((jPos) => `${jPos + 1}`);

        this.startPosition = this.startPositionName && this.positionNames ?
          this.positionNames.indexOf(this.startPositionName) : 0;
        this.endPosition = this.endPositionName && this.positionNames ?
          this.positionNames.indexOf(this.endPositionName) : maxLength;

        //#region -- palette --
        switch (this.seqCol.semType) {
        case Aminoacids.SemTypeMultipleAlignment:
          this.cp = AminoacidsPalettes.GrokGroups;
          break;
        case Nucleotides.SemTypeMultipleAlignment:
          this.cp = NucleotidesPalettes.Chromatogram;
          break;
        }
        //#endregion -- palette --
      } else {
        this.cp = null;
        this.positionNames = [];
        this.startPosition = -1;
        this.endPosition = -1;
      }
    }
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);

    switch (property.name) {
    case 'considerNullSequences':
      this.render();
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
    }
  }

  onTableAttached() {
    this.updateSeqCol();

    if (this.dataFrame !== void 0) {
      // There are two approaches:
      // first  - look in the dataFrame for the first matching column by semType of the
      //          corresponding viewer (but we want only one class of a more universal viewer),
      // second - draw column data if the passed column is of suitable semType
      // We decided that we will not search, but we will display asked data if we can
      // const semType = (<typeof NucleotidesWebLogo>(this.constructor)).residuesSet;
      // this.seqCol = (this.dataFrame.columns as DG.ColumnList).bySemType(semType);

      this.dataFrame.selection.onChanged.subscribe((_) => this.render());
      this.dataFrame.filter.onChanged.subscribe((_) => this.render());
    }
  }

  protected _nullSequence(fillerResidue = 'X'): string {
    if (this.considerNullSequences)
      return new Array(this.Length).fill(fillerResidue).join('');

    return '';
  }

  protected _calculate() {
    if (!this.seqCol || !this.dataFrame || this.startPosition === -1 || this.endPosition === -1)
      return;

    const width = this.Length * this.positionWidth;
    this.canvas.width = width;
    this.canvas.style.width = `${width}px`;

    this.positions = new Array(this.Length);
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

      for (let jPos = 0; jPos < this.Length; jPos++) {
        const pmInfo = this.positions[jPos].freq;
        const m: string = s[this.startPosition + jPos];
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
    }
    //#endregion

    const maxHeight = this.canvas.height - this.axisHeight;
    // console.debug(`WebLogo._calculate() maxHeight=${maxHeight}.`);

    //#region Calculate screen
    for (let jPos = 0; jPos < this.Length; jPos++) {
      const freq: { [c: string]: PositionMonomerInfo } = this.positions[jPos].freq;
      const rowCount = this.positions[jPos].rowCount;

      let y: number = this.axisHeight;

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

        pmInfo.bounds = new DG.Rect(jPos * this.positionWidth, y, this.positionWidth, h);
        y += h;
      }
    }
    //#endregion
  }

  // reflect changes made to filter/selection
  render(recalc = true) {
    const g = this.canvas.getContext('2d');
    if (!this.seqCol || !this.dataFrame || !g || !this.cp || this.startPosition === -1 || this.endPosition === -1)
      return;

    if (recalc)
      this._calculate();

    // let rowCount = this.rowsMasked;
    // if (!this.considerNullSequences)
    //   rowCount -= this.rowsNull;

    g.resetTransform();
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    g.textBaseline = this.textBaseline;

    for (let jPos = 0; jPos < this.Length; jPos++) {
      const pos: PositionInfo = this.positions[jPos];
      g.resetTransform();
      g.fillStyle = 'black';
      g.textAlign = 'center';
      g.font = '10px Roboto, Roboto Local, sans-serif';
      const posNameTm: TextMetrics = g.measureText(pos.name);
      const jPosWidth = posNameTm.width < (this.positionWidth - 2) ? posNameTm.width : (this.positionWidth - 2);
      g.setTransform(
        jPosWidth / posNameTm.width, 0, 0, 1,
        jPos * this.positionWidth + this.positionWidth / 2, 0);
      g.fillText(pos.name, 0, 0);

      for (const [monomer, pmInfo] of Object.entries(this.positions[jPos].freq)) {
        if (monomer !== '-') {
          const b = pmInfo.bounds;

          const fontStyle = '16px Roboto, Roboto Local, sans-serif';
          // Hacks to scale uppercase characters to target rectangle
          const uppercaseLetterAscent = 0.25;
          const uppercaseLetterHeight = 12.2;

          g.resetTransform();
          g.strokeStyle = 'lightgray';
          g.lineWidth = 1;
          g.rect(b.left, b.top, b.width, b.height);
          g.fillStyle = this.cp[monomer] ?? this.cp['other'];
          g.textAlign = 'left';
          g.font = fontStyle;
          //g.fillRect(b.left, b.top, b.width, b.height);
          const mTm: TextMetrics = g.measureText(monomer);

          // if (mM.actualBoundingBoxAscent != 0)
          //   console.debug(`m: ${m}, mM.actualBoundingBoxAscent: ${mM.actualBoundingBoxAscent}`);

          g.setTransform(
            b.width / mTm.width, 0, 0, b.height / uppercaseLetterHeight,
            b.left, b.top);
          g.fillText(monomer, 0, -uppercaseLetterAscent);
        }
      }
    }
  }
}
