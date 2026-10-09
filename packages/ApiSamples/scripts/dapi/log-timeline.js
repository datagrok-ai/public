// What happened in the current session in the last 15 minutes, oldest first: its events and HTTP requests.
// `to` defaults to `from` + 10 minutes, at most 2 hours later. Another session needs the ViewTelemetry permission.

let session = (await grok.dapi.users.currentSession()).id;
let timeline = await grok.dapi.log.getTimeline({session, from: new Date(Date.now() - 15 * 60000), to: new Date()});
grok.shell.addTableView(DG.DataFrame.fromObjects(timeline) ?? DG.DataFrame.create());
