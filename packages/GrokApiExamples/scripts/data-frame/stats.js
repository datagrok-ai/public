// Column statistics

var demog = gr.testData('demog', 5000);
var view = gr.newView('stats');

for (var col of demog.columns.toList()) {
    view.append(ui.tableFromMap({
        name: col.name,
        totalCount: col.stats.totalCount,
        missingValueCount: col.stats.missingValueCount,
        valueCount: col.stats.valueCount,

        min: col.stats.min,
        max: col.stats.max,
        avg: col.stats.avg,
        stdev: col.stats.stdev,
        variance: col.stats.variance,
        skew: col.stats.skew,
        kurt: col.stats.kurt,
        med: col.stats.med,
        q1: col.stats.q1,
        q2: col.stats.q2,
        q3: col.stats.q3
    }))
}