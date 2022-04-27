import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

type ColorScheme = { [code: string]: string };

const ColorSchemes: { [residueSet: string]: ColorScheme } = {
  // nucleotides: {'A': 'red', 'G': 'orange', 'T': 'green', 'C': 'blue'},

  nucleotides: {
    // Chromatograms nucleotide colours
    'A': 'green',
    'C': 'blue',
    'G': 'black', // orange ?
    'T': 'red', 'U': 'red',
    'others': '',
  },

  aminoacids: {
    // RasMol Amino Colours
    // http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
    'D': '#E60A0A', // asp, aspartic acid, asp
    'E': '#E60A0A', // glu, glutamic acid
    'C': '#E6E600', // cys, cysteine
    'M': '#E6E600', // met, methionine
    'K': '#145AFF', // lys, lysine
    'R': '#145AFF', // arg, arginine
    'S': '#FA9600', // ser, serine
    'T': '#FA9600', // thr, threonine
    'F': '#3232AA', // phe, phenylalanine
    'Y': '#3232AA', // tyr, tyrosine
    'N': '#00DCDC', // asn, asparagine
    'Q': '#00DCDC', // gln, glutamine
    'G': '#EBEBEB', // gly, glycine
    'L': '#0F820F', // leu, leucine
    'V': '#0F820F', // val, valine
    'I': '#0F820F', // ile, isoleucine
    'A': '#C8C8C8', // ala, alanine
    'W': '#B45AB4', // trp, tryptophan
    'H': '#8282D2', // his, histidine
    'P': '#DC9682', // pro, proline
    'others': '#BEA06E',
  },
};

export class NucleotidesWebLogo extends DG.JsViewer {
  public static residuesSet = 'nucleotides';

  private readonly colorScheme: ColorScheme = ColorSchemes[NucleotidesWebLogo.residuesSet];
  private readonly canvas: HTMLCanvasElement;
  private readonly textBaseline: CanvasTextBaseline;

  private maxHeight: number;
  private seqCol: DG.Column = null;
  private maxLength: number;
  private freqs: { freq: { [c: string]: number, }, rowCount: number; }[];

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
    this.colorScheme = ColorSchemes[semType];

    this.considerNullSequences = this.bool('considerNullSequences', false);
    this.sequenceColumnName = this.string('sequenceColumnName', null);

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

    this.freqs = new Array(this.maxLength);
    for (let i = 0; i < this.maxLength; i++)
      this.freqs[i] = {freq: {}, rowCount: 0};

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

      for (let j = 0; j < s.length; j++) {
        const f: { [c: string]: number } = this.freqs[j].freq;
        const c: string = s[j];
        f[c] = (c in f) ? f[c] + 1 : 1;
      }
    }

    //#region Polish freq counts
    for (let j = 0; j < this.freqs.length; j++) {
      delete this.freqs[j].freq['-'];

      this.freqs[j].rowCount = 0;
      for (const c in this.freqs[j].freq)
        this.freqs[j].rowCount += this.freqs[j].freq[c];
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
    const maxHeight = this.canvas.clientHeight;

    // let rowCount = this.rowsMasked;
    // if (!this.considerNullSequences)
    //   rowCount -= this.rowsNull;

    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    g.textBaseline = this.textBaseline;

    for (let j = 0; j < this.maxLength; j++) {
      const freq: { [c: string]: number } = this.freqs[j].freq;
      const rowCount = this.freqs[j].rowCount;

      // positions
      // g.setTransform(1, 0, 0, 1, 0, 0);
      g.resetTransform();
      g.fillStyle = 'black';
      g.fillText(j.toString(), j * 16, 0);
      let y: number = 16;

      const entries = Object.entries(freq).sort((a, b) => b[1] - a[1]);
      for (const entry of entries) {
        const count: number = entry[1];
        const c: string = entry[0];
        const h: number = maxHeight * count / rowCount;
        g.setTransform(1, 0, 0, h / 16, 0, y);
        g.fillStyle = this.colorScheme[c] ?? this.colorScheme['others'];
        g.fillText(c, j * 16, 0);
        y += h;
      }
    }
  }
}

export class AminoacidsWebLogo extends NucleotidesWebLogo {
  public static override residuesSet = 'aminoacids';
}
