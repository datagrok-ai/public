import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const ColorSchemes: {[residueSet: string]: ColorScheme} = {
  nucleotides: {'A': 'red', 'G': 'orange', 'T': 'green', 'C': 'blue'},
  aminoacids: {
    'A': 'red',
    'G': 'orange',
    'T': 'green',
    'C': 'blue',
    'X': 'gray',
    'other': 'black',
  },
};

export class NucleotidesWebLogo extends DG.JsViewer {
  static residuesSet = 'nucleotides';

  colorScheme: ColorScheme = ColorSchemes[NucleotidesWebLogo.residuesSet];
  canvas: HTMLCanvasElement;
  textBaseline: CanvasTextBaseline;
  maxHeight: number;
  seqCol: DG.Column;
  maxLength: number;
  freq: any[];
  rowsMasked: number;
  rowsNull: number;
  considerNullSequences: boolean;

  constructor() {
    super();
    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';
    this.canvas.style.height = '50px';
    this.textBaseline = 'top';
    this.maxHeight = 50;

    this.considerNullSequences = this.bool('considerNullSequences', false);

    this.root.appendChild(this.canvas);
    ui.onSizeChanged(this.canvas).subscribe(this.onSizeChanged.bind(this));
  }

  onSizeChanged(args: any) {
    this.canvas.width = this.canvas.clientWidth;
    this.canvas.height = this.canvas.clientHeight;
    this.render(false);
  }

  onPropertyChanged(property: DG.Property): void {
    if (property.name == 'considerNullSequences')
      this.render();
  }

  onTableAttached() {
    const semType = (<typeof NucleotidesWebLogo> this.constructor).residuesSet;
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

    this.freq = new Array(this.maxLength);
    for (let i = 0; i < this.maxLength; i++)
      this.freq[i] = {};

    const indices = this.dataFrame.selection.trueCount > 0 ?
      this.dataFrame.selection.getSelectedIndexes() :
      this.dataFrame.filter.getSelectedIndexes();

    this.rowsMasked = indices.length;
    this.rowsNull = 0;

    for (const i of indices) {
      let s = this.seqCol.get(i);

      if (!s) {
        s = this._nullSequence();
        ++this.rowsNull;
      }

      for (let j = 0; j < s.length; j++) {
        const f = this.freq[j];
        const c = s[j];
        f[c] = f[c] == null ? 1 : f[c] + 1;
      }
    }
  }

  // reflect changes made to filter/selection
  render(recalc = true) {
    if (!this.seqCol)
      return;

    if (recalc)
      this._calculate();

    const g = this.canvas.getContext('2d');
    let rowCount = this.rowsMasked;

    if (!this.considerNullSequences)
      rowCount -= this.rowsNull;

    g.resetTransform();
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    g.textBaseline = this.textBaseline;

    for (let i = 0; i < this.maxLength; i++) {
      const freq = this.freq[i];
      let y = 0;
      for (const c of Object.keys(freq)) {
        const count = freq[c];
        const maxHeight = this.canvas.clientHeight;
        const h = maxHeight * count / rowCount;
        g.setTransform(1, 0, 0, h / 16, 0, y);
        g.fillStyle = this.colorScheme[c] ?? this.colorScheme.other;
        g.fillText(c, i * 16, 0);
        y += h;
      }
    }
  }
}

type ColorScheme = {[code: string]: string};

export class AminoacidsWebLogo extends NucleotidesWebLogo {
  static residuesSet = 'aminoacids';
}
