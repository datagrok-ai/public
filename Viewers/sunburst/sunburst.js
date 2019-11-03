class SunburstViewer extends JsViewer {
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

        this.table = new DataFrame(dataFrameHandle);
        this.init();

        this.table.selection.onChanged(() => this.render());
        this.table.filter.onChanged(() => this.render());
        this.render();
    }

    render() {
    }
}
