class WebLogoViewer extends DG.JsViewer {

  constructor() {
    super();
    this.canvas = ui.canvas();
    this.canvas.style.width = '100%';
    this.canvas.style.height = '50px';

    this.root.appendChild(this.canvas);
    ui.onSizeChanged(this.canvas).subscribe((_) => {
      this.canvas.width = this.canvas.clientWidth;
      this.canvas.height = this.canvas.clientHeight;
      this.render(false);
    });
  }

  onTableAttached() {
    this.seqCol = this.dataFrame.columns.bySemType('nucleotides');
    this.dataFrame.selection.onChanged.subscribe((_) => this.render());
    this.dataFrame.filter.onChanged.subscribe((_) => this.render());
    this.render();
  }

  _calculate() {
    if (!this.seqCol)
      return;

    this.maxLength = 0;
    for (let category of this.seqCol.categories)
      this.maxLength = Math.max(this.maxLength, category.length);

    this.freq = new Array(this.maxLength);
    for (let i = 0; i < this.maxLength; i++)
      this.freq[i] = {};

    this.indexes = this.dataFrame.selection.trueCount > 0
      ? this.dataFrame.selection.getSelectedIndexes()
      : this.dataFrame.filter.getSelectedIndexes()

    for (let i of this.indexes) {
      const s = this.seqCol.get(i);
      for (let j = 0; j < s.length; j++) {
        const f = this.freq[j]
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

    const colors = {A: 'red', 'G': 'orange', T: 'green', C: 'blue'};

    const g = this.canvas.getContext('2d');
    const maxHeight = 50;
    const rowCount = this.indexes.length;
    g.resetTransform();
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    g.textBaseline = 'top';

    for (let i = 0; i < this.maxLength; i++) {
      let y = 0;
      for (const c in this.freq[i]) {
        const count = this.freq[i][c];
        const h = maxHeight * count / rowCount;
        g.setTransform(1, 0, 0, h / 16, 0 , y);
        g.fillStyle = colors[c];
        g.fillText(c, i * 16, 0);
        y += h;
      }
    }
  }
}