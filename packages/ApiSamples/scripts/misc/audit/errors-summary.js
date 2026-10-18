// The errors of this week, grouped by stack signature, most users first. Needs the ViewTelemetry permission.
let errors = (await grok.dapi.admin.getMetrics({date: 'this week', errorsBy: 'signature', limit: 15})).errors;
grok.shell.addTableView(DG.DataFrame.fromObjects(errors.top));
