class MultiFormViewer extends DG.JsViewer {
  constructor() {
    super();

    // properties
    this.fields = this.addProperty('fields', DG.TYPE.COLUMN_LIST);
    this.colorCode = this.bool('colorCode', true);
    this.showCurrentRow = this.bool('showCurrentRow', true);
    this.showMouseOverRow = this.bool('showMouseOverRow', true);

    //fields
    this.indexes = [];

    // init
    let dis = this;
    this.root.classList.add('d4-multi-form');
    this.columnHeadersDiv = ui.div([], 'd4-multi-form-header');
    this.virtualView = ui.virtualView(0, i => dis.renderForm(i), false, 1);
    this.root.appendChild(ui.divH([this.columnHeadersDiv, this.virtualView.root]));
  }

  onTableAttached() {
    if (this.fields == null)
      this.fields = this.dataFrame.columns.names();

    const sub = (stream, action) => {
      this.subs.push(DG.debounce(stream, 50).subscribe((_) => action()));
    }

    sub(this.dataFrame.selection.onChanged, () => this.render());
    sub(this.dataFrame.filter.onChanged, () => this.render());

    sub(this.dataFrame.onCurrentRowChanged.pipe(rxjs.operators.filter((_) => this.showCurrentRow)),
      () => this.virtualView.refreshItem(this.currentRowPos));

    sub(this.dataFrame.onMouseOverRowChanged.pipe(rxjs.operators.filter((_) => this.showMouseOverRow)),
      () => this.virtualView.refreshItem(this.mouseOverPos));

    this.render();
  }

  onPropertyChanged(prop) {
    this.render();
  }

  renderHeader() {
    ui.empty(this.columnHeadersDiv);
    const form = this.renderForm(0);
    form.classList.add('temp');
    document.body.appendChild(form);

    for (let name of this.fields) {
      const columnLabel = ui.bind(this.dataFrame.col(name), ui.divText(name, 'd4-multi-form-column-name'));
      let formField = form.querySelector('[column="' + name + '"]');
      columnLabel.style.top = `${formField.offsetTop + 10}px`;
      this.columnHeadersDiv.appendChild(columnLabel);
    }

    form.remove();
  }

  get currentRowPos() { return this.showCurrentRow ? 0 : null; }
  get mouseOverPos() { return this.showMouseOverRow ? (this.showCurrentRow ? 1 : 0) : null; }

  renderForm(row) {
    if (this.showCurrentRow && row === this.currentRowPos)
      row = this.dataFrame.currentRowIdx;
    else if (this.showMouseOverRow && row === this.mouseOverPos)
      row = this.dataFrame.mouseOverRowIdx;
    else
      row = this.indexes[row - (this.showCurrentRow ? 1 : 0) - (this.showMouseOverRow ? 1 : 0)];

    const form = ui.divV(
      this.fields.map((name) => {
        if (row === -1)
          return ui.div();

        const input = DG.InputBase.forColumn(this.dataFrame.col(name));
        input.input.setAttribute('column', name);
        input.value = this.dataFrame.get(name, row);

        if (this.colorCode) {
          const grid = this.view.grid;
          const color = grid.cell(name, row).color;
          input.input.style.color = DG.Color.toHtml(DG.Color.getContrastColor(color));
          input.input.style.backgroundColor = DG.Color.toHtml(color);
        }

        return input.input;
      }
    ), 'd4-multi-form-form');

    form.onclick = () => this.dataFrame.currentRowIdx = row;
    form.onmouseenter = () => this.dataFrame.mouseOverRowIdx = row;
    form.onmouseleave = () => this.dataFrame.mouseOverRowIdx = -1;
    return form;
  }

  render() {
    let dis = this;
    this.indexes
      = this.dataFrame.selection.trueCount > 0 ? this.dataFrame.selection.getSelectedIndexes()
      : [];

    ui.empty(this.columnHeadersDiv);

    this.renderHeader();

    this.virtualView.setData(
      this.indexes.length + (this.showCurrentRow ? 1 : 0) + (this.showMouseOverRow ? 1 : 0),
      i => dis.renderForm(i));
  }
}
