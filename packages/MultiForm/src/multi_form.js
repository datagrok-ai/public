class MultiFormViewer extends DG.JsViewer {
  constructor() {
    super();

    // properties
    this.fields = this.addProperty('fields', DG.TYPE.COLUMN_LIST);
    this.colorCode = this.bool('colorCode', true);

    // init
    let dis = this;
    this.root.classList.add('d4-multi-form');
    this.columnHeadersDiv = ui.div([], 'd4-multi-form-header');
    this.virtualView = ui.virtualView(0, function (i) {
      return dis.renderForm(i);
    }, false, 1);
    this.root.appendChild(ui.divH([this.columnHeadersDiv, this.virtualView.root]));

    //ui.onSizeChanged(this.root).subscribe((_) => grok.shell.info('resized'));
  }

  onTableAttached() {
    if (this.fields == null)
      this.fields = this.dataFrame.columns.names();

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 20).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));

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
      columnLabel.style.top = `${formField.offsetTop}px`;
      this.columnHeadersDiv.appendChild(columnLabel);
    }

    form.remove();
  }

  renderForm(row) {
    const form = ui.divV(
      this.fields.map((name) => {
        const input = DG.InputBase.forColumn(this.dataFrame.col(name));
        input.input.setAttribute('column', name);
        input.value = this.dataFrame.get(name, row);

        if (this.colorCode) {
          const grid = this.view?.grid;
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
    let indexes
      = this.dataFrame.selection.anyTrue ? this.dataFrame.selection.getSelectedIndexes()
      : this.dataFrame.currentRowIdx > -1 ? [this.dataFrame.currentRowIdx]
      : this.dataFrame.rowCount > 0 ? [0]
      : [];

    ui.empty(this.columnHeadersDiv);

    this.renderHeader();

    this.virtualView.setData(indexes.length, function (i) {
      return dis.renderForm(i);
    });
  }
}
