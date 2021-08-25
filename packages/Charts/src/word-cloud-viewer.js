import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as echarts from 'echarts';
import 'echarts-wordcloud';
// See also https://datagrok.ai/help/develop/how-to/develop-custom-viewer
// This viewer does the following:
// * listens to changes of filter and selection in the attached table,
// * updates the number of filtered/selected rows accordingly.
export class WordCloudViewer extends DG.JsViewer {
    constructor() {
        super();
        // Register properties and define fields initialized to properties' default values
        // Properties that represent columns should end with the 'ColumnName' postfix
        this.strColumnName = this.string('columnColumnName');
        this.strColumns = [];
        this.initialized = false;
    }
    init() {
        this.initialized = true;
    }
    testColumns() {
        return (this.strColumns.length >= 1);
    }
    onTableAttached() {
        this.init();
        let columns = this.dataFrame.columns.toList();
        this.strColumns = columns.filter((col) => col.type === DG.TYPE.STRING);
        if (this.testColumns()) {
            this.strColumnName = this.strColumns.reduce((prev, curr) => prev.categories.length < curr.categories.length ? prev : curr).name;
        }
        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));
        this.render();
    }
    onPropertyChanged(property) {
        super.onPropertyChanged(property);
        if (this.initialized && this.testColumns()) {
            console.log(this.strColumnName);
            if (property.name === 'columnColumnName') {
                this.strColumnName = property.get();
            }
            this.render();
        }
    }
    detach() {
        this.subs.forEach(sub => sub.unsubscribe());
    }
    render() {
        if (!this.testColumns()) {
            this.root.innerText = "Not enough data to produce the result.";
            return;
        }
        $(this.root).empty();
        let margin = { top: 10, right: 10, bottom: 10, left: 10 };
        let width = this.root.parentElement.clientWidth - margin.left - margin.right;
        let height = this.root.parentElement.clientHeight - margin.top - margin.bottom;
        let strColumn = this.dataFrame.getCol(this.strColumnName);
        let words = strColumn.categories;
        let data = [];
        words.forEach(w => data.push({ name: w, value: strColumn.toList().filter(row => row === w).length }));
        let fontColor = () => {
            let nColors = DG.Color.categoricalPalette.length;
            let color = DG.Color.getCategoricalColor(Math.floor(Math.random() * nColors));
            return DG.Color.toRgb(color);
        };
        let table = this.dataFrame;
        if (this.chart !== undefined) {
            this.chart.dispose();
        }
        this.chart = echarts.init(this.root);
        this.chart.setOption({
            width: width + margin.left + margin.right,
            height: height + margin.top + margin.bottom,
            series: [{
                    type: 'wordCloud',
                    shape: 'circle',
                    left: 'center',
                    top: 'center',
                    width: `${width}`,
                    height: `${height}`,
                    right: null,
                    bottom: null,
                    sizeRange: [14, 100],
                    gridSize: 15,
                    rotationRange: [0, 0],
                    rotationStep: 0,
                    drawOutOfBound: true,
                    textStyle: {
                        fontFamily: 'sans-serif',
                        fontWeight: 'bold',
                        // Color can be a callback function or a color string
                        color: fontColor
                    },
                    emphasis: {
                        focus: 'self',
                        textStyle: {
                            shadowBlur: 5,
                            shadowColor: '#333'
                        }
                    },
                    data: data
                }]
        });
        this.chart
            .on('mouseover', d => ui.tooltip.showRowGroup(table, i => {
            console.log(d);
            return d.name === strColumn.get(i);
        }, 10, 10))
            .on('mouseout', () => ui.tooltip.hide())
            .on('mousedown', d => {
            table.selection.handleClick(i => {
                return d.name === strColumn.get(i);
            }, d);
        });
    }
}
