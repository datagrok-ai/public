let result = await grok.dapi.admin.getServiceInfos();
let v = grok.shell.newView('list');

v.root.appendChild(
  ui.table(result, (item, idx) => [`${item.key}:`, item.status])
);

// Server metrics for a window: request latency per route, the errors recorded, the function-call
// queue, and database statistics. `date` is any Datagrok datetime pattern.
let metrics = await grok.dapi.admin.getMetrics({date: 'last 7 days', limit: 5});
let top = metrics.http.routes[0];
grok.shell.info(`p95 ${metrics.http.now.p95} ms over ${metrics.http.now.count} requests; ` +
  (top ? `slowest route: ${top.method} ${top.route} (p95 ${top.p95} ms)` : 'no routes'));
// The top error group: what most users hit and what keeps recurring comes first.
let error = metrics.errors.top[0];
grok.shell.info(`${metrics.errors.now.count} errors hit ${metrics.errors.now.users} users; ` +
  (error ? `worst: "${error.message}" (${error.count} times, ${error.users} users, ${error.hours} h)` : 'none'));
