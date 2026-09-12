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
