class SunburstViewer extends grok.JsViewer {
    init() {
        let data = {
            name: "root",
            children: [
                {
                    name: "leafA", value: 3
                },
                {
                    name: "nodeB",
                    children: [{name: "leafBA", value: 5}, {name: "leafBB", value: 1}]
                }
            ]
        };

        const myChart = Sunburst();
        myChart.data(data)(this.root);
    }

    onFrameAttached(dataFrameHandle) {
        // https://stackoverflow.com/questions/43858714/typeerror-object-entries-is-not-a-function
        if (!Object.entries) {
            Object.entries = function (obj) {
                var ownProps = Object.keys(obj),
                    i = ownProps.length,
                    resArray = new Array(i); // preallocate the Array
                while (i--)
                    resArray[i] = [ownProps[i], obj[ownProps[i]]];
                return resArray;
            };
        }

        this.table = new grok.DataFrame(dataFrameHandle);
        this.init();

        this.subs.push(this.table.selection.onChanged.subscribe((_) => this.render()));
        this.subs.push(this.table.filter.onChanged.subscribe((_) => this.render()));
        this.render();
    }

    render() {
    }
}
