import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as d3 from 'd3';
import * as cloud from 'd3-cloud';


export class WordCloudViewer extends DG.JsViewer {

    constructor() {
        super();

        // Properties
        this.strColumnName = this.string('strColumnName');

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
        this.strColumns = columns.filter(col => col.type === 'string');

        if (this.testColumns()) {
            // Find a string column with the smallest number of unique values
            this.strColumnName = this.strColumns.reduce((prev, curr) =>
                                 prev.categories.length < curr.categories.length ? prev : curr).name;
        }

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.onSizeChanged, 50).subscribe((_) => this.render()));

        this.render();
    }

    onPropertyChanged(property) {
        super.onPropertyChanged(property);
        if (this.initialized && this.testColumns()) {
            this.render();
        }
    }

    detach() {
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    render() {

        if (!this.testColumns()) {
            this.root.innerText = "Not enough data to produce the result.";
            return;
        }

        this.root.classList.add('viewer-window');
        let margin = {top: 10, right: 10, bottom: 10, left: 10};
        let width = this.root.clientWidth - margin.left - margin.right;
        let height = this.root.clientHeight - margin.top - margin.bottom;

        let svg = d3.select(this.root).append("svg")
                        .attr("width", width + margin.left + margin.right)
                        .attr("height", height + margin.top + margin.bottom);

        svg
            .append("g")
            .attr("transform", `translate(${margin.left}, ${margin.top})`);

        let strColumn = this.dataFrame.getCol(this.strColumnName);
        let words = strColumn.categories;
        this.freqMap = {};
        strColumn.toList().forEach(w => this.freqMap[w] = (this.freqMap[w] || 0) + 1);
        let sortedWords = Object.keys(this.freqMap).sort((a, b) => {
            this.freqMap[b] - this.freqMap[a];
        });
        this.max = this.freqMap[sortedWords[0]];  

        let layout = cloud()
            .size([width, height])
            .words(words.map(d => {
                return { text: d,
                         size: (this.freqMap[d] * (100 - 10))/this.max + 10,
                         color: DG.Color.toRgb(DG.Color.getCategoricalColor(Math.floor(
                                Math.random() * DG.Color.categoricalPalette.length)))
                };
            }))
            .padding(10)
            .rotate(() => (~~(Math.random() * 6) - 3) * 30)
            .fontSize(d => d.size)
            .on("end", draw);
        layout.start();

        function draw(words) {
            svg
                .append("g")
                .attr("transform", `translate(${layout.size()[0] / 2}, ${layout.size()[1] / 2})`)
                .selectAll("text")
                    .data(words)
                .enter().append("text")
                    .style("font-size", d => d.size + "px")
                    .style("fill", d => d.color)
                    .attr("text-anchor", "middle")
                    .attr("transform", d => `translate(${[d.x, d.y]}) rotate(${d.rotate})`)
                    .text(d => d.text);
        }

        this.root.appendChild(svg.node());
    }
}
