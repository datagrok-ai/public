// What happened in one request, oldest first: its events and the request itself.
// An action id instead ({action: id}) also brings in the action's requests (`<action>.<n>`).
// A server-side event carries the id of the request it was logged in. Needs the ViewTelemetry permission.

let events = await grok.dapi.log.list({pageSize: 100});
let event = events.find((e) => e.requestId != null);
if (event == null)
  grok.shell.info('No recent events logged in a request');
else {
  let timeline = await grok.dapi.log.getTimeline({request: event.requestId});
  grok.shell.addTableView(DG.DataFrame.fromObjects(timeline) ?? DG.DataFrame.create());
}
