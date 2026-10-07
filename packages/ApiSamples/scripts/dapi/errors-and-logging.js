// The platform's errors as data, the timeline of a session and the logging policy.
// Errors and the timeline need the ViewTelemetry permission, the policy EditPluginsSettings.

// The signatures that hit the most users in the last day
let top = await grok.dapi.log.getErrors({since: '1d', by: 'signature', limit: 10});
grok.shell.info(top.map((r) => `${r.signature.substring(0, 6)}: ${r.count} times, ${r.users} users`).join('\n') || 'No errors');

// What the current session did, oldest first
let session = await grok.dapi.users.currentSession();
let timeline = await grok.dapi.log.getTimeline({session: session.id, limit: 100});
if (timeline.length > 0)
  grok.shell.addTableView(DG.DataFrame.fromObjects(timeline));

// The debug flags a capture rule can turn on
let policy = await grok.dapi.log.getLoggingPolicy();
grok.shell.info(`Debug flags: ${policy.debugFlags.join(', ')}`);
