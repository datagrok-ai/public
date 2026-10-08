// What the current session did, oldest first: its events, requests and function calls.
// Needs the ViewTelemetry permission.

let session = await grok.dapi.users.currentSession();
let timeline = await grok.dapi.log.getTimeline({session: session.id, limit: 100});
grok.shell.addTableView(DG.DataFrame.fromObjects(timeline) ?? DG.DataFrame.create());
