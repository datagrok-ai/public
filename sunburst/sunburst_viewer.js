class Sunburst2Viewer extends DG.JsViewer {
    constructor() {
        super();

        // properties
        this.level1ColumnName = this.string('level1ColumnName');
        this.level2ColumnName = this.string('level2ColumnName');
        this.level3ColumnName = this.string('level3ColumnName');
        this.valueColumnName = this.string('valueColumnName');

        this.treeData = {};
        this.containerId = "sunburst-id-123456";
        this.containerDiv = null;

        this.clearTreeData();
    }

    init() {
        this.containerDiv = ui.div([], 'd4-viewer-host');
        this.containerDiv.setAttribute("id", this.containerId);
        this.root.appendChild(this.containerDiv);
    }

    onTableAttached() {
        this.init();

        // this.level1ColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);
        // this.level2ColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);
        // this.level3ColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);
        // this.valueColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));

        this.render(true);
    }

    onPropertyChanged(prop) {
        this.render();
    }

    onSizeChanged(w, h) {
        // if (w !== 0 && h !== 0)
        //     this.map.invalidateSize();
    }

    detach() {
        this.map.remove();
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    clearTreeData() {
        // for (var x in this.treeData) if (this.treeData.hasOwnProperty(x)) delete this.treeData[x];
        this.treeData.id = d3.hierarchy({ id: this.containerId, parent: '', children: [], data: {}});
    }

    getKey(parent, id) {
        return JSON.stringify({ parent, id })
    }

    buildTreeData() {
        this.clearTreeData();

        const indexes = this.dataFrame.filter.getSelectedIndexes();
        console.error(indexes);

        const columns = [];
        this.level1ColumnName && columns.push(this.dataFrame.getCol(this.level1ColumnName));
        this.level2ColumnName && columns.push(this.dataFrame.getCol(this.level2ColumnName));
        this.level3ColumnName && columns.push(this.dataFrame.getCol(this.level3ColumnName));
        const valueColumn = this.valueColumnName && this.dataFrame.getCol(this.valueColumnName);
        console.error(columns, valueColumn);

        if (!columns.length) {
            return;
        }

        const table = new Map();
        table.set(this.getKey('', this.containerId), 0);
        const [lastColumn] = columns.slice(-1);
        for (const index of indexes) {
            let prevPrevColumnValue = null;
            let prevColumnValue = this.containerId;
            for (const column of columns) {
                let columnValue = column.get(index);
                let isLast = column === lastColumn;
                if (columnValue === '') {
                    if (prevPrevColumnValue == null) {
                        break;
                    }
                    isLast = true;
                    columnValue = prevColumnValue;
                    prevColumnValue = prevPrevColumnValue;
                }
                const key = this.getKey(prevColumnValue, columnValue);
                const value = table.get(key) || 0;
                const valueIncrement = isLast ? (valueColumn ? valueColumn.get(index) : 1) : 0;
                table.set(key, value + valueIncrement);
                if (isLast) {
                    break;
                }
                prevPrevColumnValue = prevColumnValue;
                prevColumnValue = columnValue;
            }
        }
        console.error(table);

        const list = [];
        table.forEach((value, k) => list.push({...JSON.parse(k), value}));
        console.error(list);

        const treeData = d3.stratify()
            .id(d => d.id)
            .parentId(d => d.parent)
            (list)
            .sum(d => d.value);
        this.treeData = treeData;
    }

    render() {
        console.error('rendering!', this);
        this.buildTreeData();
        this.containerDiv.innerHTML = '';
        d3sunburst(this.containerId, this.treeData);
    }
}
