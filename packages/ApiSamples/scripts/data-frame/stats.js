// Column statistics
let demog = grok.data.demo.demog();
let view = grok.shell.newView('stats');

for (let col of demog.columns.toList()) {
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
  }));
}
