import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {ChemPalette} from '../utils/chem-pallete';
import {Point, Rect} from 'datagrok-api/dg';
import {fromEvent} from 'rxjs';

// Using color schemes from chem-palette

declare module 'datagrok-api/src/grid' {
  interface Rect {
    contains(x: number, y: number): boolean;
  }
}

declare global {
  interface HTMLCanvasElement {
    getCursorPosition(event: MouseEvent): Point;
  }
}

HTMLCanvasElement.prototype.getCursorPosition = function(event: MouseEvent): Point {
  const rect = this.getBoundingClientRect();
  return new Point(event.clientX - rect.left, event.clientY - rect.top);
};

Rect.prototype.contains = function(x: number, y: number): boolean {
  return this.left <= x && x <= this.right && this.top <= y && y <= this.bottom;
};

class PositionMonomerInfo {
  /** Sequences count with monomer in position
   */
  count: number;

  /** Remember screen coords rect
   */
  bounds: DG.Rect;
}

class PositionInfo {
  freq: { [m: string]: PositionMonomerInfo };
  rowCount: number;
}

export class NucleotidesWebLogo extends DG.JsViewer {
  public static residuesSet = 'nucleotides';

  // private readonly colorScheme: ColorScheme = ColorSchemes[NucleotidesWebLogo.residuesSet];
  protected cp: StringDictionary;
  private readonly canvas: HTMLCanvasElement;
  private readonly textBaseline: CanvasTextBaseline;

  private maxHeight: number;
  private axisHeight: number = 14;

  private positionWidth: number = 16;

  private seqCol: DG.Column = null;
  private maxLength: number;
  private positions: PositionInfo[];

  private rowsMasked: number;
  private rowsNull: number;

  // Viewer's properties (likely they should be public so that they can be set outside)
  public considerNullSequences: boolean;
  public sequenceColumnName: string;

  constructor() {
    super();
    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';
    this.canvas.style.height = '60px';
    this.textBaseline = 'top';
    this.maxHeight = 60;

    const semType = (<typeof NucleotidesWebLogo>(this.constructor)).residuesSet;
    // this.colorScheme = ColorSchemes[semType];

    this.considerNullSequences = this.bool('considerNullSequences', false);
    this.sequenceColumnName = this.string('sequenceColumnName', null);

    const getMonomer = (p: Point): [number, string, PositionMonomerInfo] => {
      const jPos = Math.floor(p.x / this.positionWidth);
      const position = this.positions[jPos];

      console.debug(`getMonomer( ${p} )`);

      if (position === void 0)
        return [jPos, null, null];

      const monomer: string = Object.keys(position.freq).find((m) => position.freq[m].bounds.contains(p.x, p.y));
      return [jPos, monomer, position[monomer]];
    };

    this.canvas.onmouseover = (e: MouseEvent) => {

    };

    fromEvent(this.canvas, 'mousedown').subscribe((e: MouseEvent) => {
      const [jPos, monomer] = getMonomer(this.canvas.getCursorPosition(e));

      // prevents deselect all rows if we miss monomer bounds
      if (monomer) {
        this.dataFrame.selection.init((iRow) => {
          const seq = this.seqCol.get(iRow);
          const mSeq = seq ? seq[jPos] : null;
          return mSeq === monomer;
        });
      }
    });

    // this.canvas.onmousedown = ;

    this.root.appendChild(this.canvas);
    ui.onSizeChanged(this.canvas).subscribe(this.onSizeChanged.bind(this));
  }

  onSizeChanged(args: any) {
    this.canvas.width = this.canvas.clientWidth;
    this.canvas.height = this.canvas.clientHeight;
    console.debug(`WebLogo.onSizeChanged() width=${this.canvas.width}, height=${this.canvas.height}`);
    this.render(false);
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);

    switch (property.name) {
    case 'considerNullSequences':
      this.render();
      break;
    case 'sequenceColumnName':
      this.seqCol = this.dataFrame.col(this.sequenceColumnName);
      this.render();
      break;
    }
  }

  onTableAttached() {
    // Trying to find a column of appropriate semType
    const semType = (<typeof NucleotidesWebLogo>(this.constructor)).residuesSet;
    this.seqCol = (this.dataFrame.columns as DG.ColumnList).bySemType(semType);

    this.dataFrame.selection.onChanged.subscribe((_) => this.render());
    this.dataFrame.filter.onChanged.subscribe((_) => this.render());
    this.render();
  }

  protected _nullSequence(fillerResidue = 'X'): string {
    if (this.considerNullSequences)
      return new Array(this.maxLength).fill(fillerResidue).join('');

    return '';
  }

  protected _calculate() {
    if (!this.seqCol)
      return;

    this.maxLength = 0;
    for (const category of this.seqCol.categories)
      this.maxLength = Math.max(this.maxLength, category.length);

    this.positions = new Array(this.maxLength);
    for (let i = 0; i < this.maxLength; i++)
      this.positions[i] = {freq: {}, rowCount: 0};

    const indices = this.dataFrame.selection.trueCount > 0 ?
      this.dataFrame.selection.getSelectedIndexes() :
      this.dataFrame.filter.getSelectedIndexes();

    this.rowsMasked = indices.length;
    this.rowsNull = 0;

    for (const i of indices) {
      let s: string = <string>(this.seqCol.get(i));

      if (!s) {
        s = this._nullSequence();
        ++this.rowsNull;
      }

      for (let jPos = 0; jPos < s.length; jPos++) {
        const pmInfo = this.positions[jPos].freq;
        const m: string = s[jPos];
        if (!(m in pmInfo))
          pmInfo[m] = {count: 0, bounds: null};
        pmInfo[m].count++;
      }
    }

    //#region Polish freq counts
    for (let jPos = 0; jPos < this.positions.length; jPos++) {
      delete this.positions[jPos].freq['-'];

      this.positions[jPos].rowCount = 0;
      for (const m in this.positions[jPos].freq)
        this.positions[jPos].rowCount += this.positions[jPos].freq[m].count;
    }
    //#endregion

    const maxHeight = this.canvas.clientHeight - this.axisHeight;

    //#region Calculate screen
    for (let jPos = 0; jPos < this.positions.length; jPos++) {
      const freq: { [c: string]: PositionMonomerInfo } = this.positions[jPos].freq;
      const rowCount = this.positions[jPos].rowCount;

      let y: number = this.axisHeight;

      const entries = Object.entries(freq).sort((a, b) => b[1].count - a[1].count);
      for (const entry of entries) {
        const pmInfo: PositionMonomerInfo = entry[1];
        // const m: string = entry[0];
        const h: number = maxHeight * pmInfo.count / rowCount;

        pmInfo.bounds = new Rect(jPos * this.positionWidth, y, this.positionWidth, h);
        y += h;
      }
    }
    //#endregion
  }

  // reflect changes made to filter/selection
  render(recalc = true) {
    if (!this.seqCol)
      return;

    if (recalc)
      this._calculate();

    const g = this.canvas.getContext('2d');

    // let rowCount = this.rowsMasked;
    // if (!this.considerNullSequences)
    //   rowCount -= this.rowsNull;

    g.resetTransform();
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    g.textBaseline = this.textBaseline;

    for (let jPos = 0; jPos < this.maxLength; jPos++) {
      g.resetTransform();
      g.fillStyle = 'black';
      g.textAlign = 'center';
      g.font = '10px Roboto, Roboto Local, sans-serif';
      g.fillText(jPos.toString(), jPos * this.positionWidth + this.positionWidth / 2, 0);

      for (const [m, pmInfo] of Object.entries(this.positions[jPos].freq)) {
        const b = pmInfo.bounds;

        const fontStyle = '16px Roboto, Roboto Local, sans-serif';
        // Hacks to scale uppercase characters to target rectangle
        const uppercaseLetterAscent = 0.25;
        const uppercaseLetterHeight = 12.2;

        g.resetTransform();
        g.strokeStyle = 'lightgray';
        g.lineWidth = 1;
        g.rect(b.left, b.top, b.width, b.height);
        g.fillStyle = this.cp[m] ?? this.cp['other'];
        g.textAlign = 'left';
        g.font = fontStyle;
        //g.fillRect(b.left, b.top, b.width, b.height);
        const mM: TextMetrics = g.measureText(m);

        if (mM.actualBoundingBoxAscent != 0)
          console.debug(`m: ${m}, mM.actualBoundingBoxAscent: ${mM.actualBoundingBoxAscent}`);

        g.setTransform(
          b.width / mM.width, 0, 0, b.height / uppercaseLetterHeight,
          b.left, b.top);
        g.fillText(m, 0, -uppercaseLetterAscent);
      }
    }
  }
}

export class AminoacidsWebLogo extends NucleotidesWebLogo {
  public static override residuesSet = 'aminoacids';

  constructor() {
    super();

    this.cp = ChemPalette.getDatagrok();
  }
}
