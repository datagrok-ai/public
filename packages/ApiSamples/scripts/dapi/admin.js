let result = await grok.dapi.admin.getServiceInfos();
let v = grok.shell.newView('list');

v.root.appendChild(
  ui.table(result, (item, idx) => [`${item.key}:`, item.status])
);

// Server metrics for a window: request latency per route, the function-call queue,
// and database statistics. `date` is any Datagrok datetime pattern.
let metrics = await grok.dapi.admin.getMetrics({date: 'last 7 days', limit: 5});
let top = metrics.http.routes[0];
grok.shell.info(`p95 ${metrics.http.now.p95} ms over ${metrics.http.now.count} requests; ` +
  (top ? `slowest route: ${top.method} ${top.route} (p95 ${top.p95} ms)` : 'no routes'));

// How the platform's data grows: per table (and category, storage nature, disk), the size now,
// its change over 1, 7 and 30 days and a daily trend. The server samples it on a schedule.
let growth = await grok.dapi.admin.getGrowthMetrics({area: 'database', days: 7});
let tables = growth.rows.filter((r) => r.table !== '(total)' && r.dimension === null && r.bytes?.d7 != null);
tables.sort((a, b) => b.bytes.d7 - a.bytes.d7);
grok.shell.info(tables.length === 0 ? 'No week of samples yet' :
  `Fastest-growing table this week: ${tables[0].table} (+${tables[0].bytes.d7} bytes, ${tables[0].bytes.now} now)`);
