// A problem the platform detected and its current alert, by `kind:key` or id; null when there is none.
// Needs the ManageAlerts permission. Error incidents are keyed by the first 12 characters of the stack hash,
// or, for errors without a stack trace, of their event type id.

let problem = await grok.dapi.admin.getProblem('error-incident:0123456789ab');
grok.shell.info(problem === null ? 'No such problem' : `${problem.status}, alert ${problem.alertStatus ?? 'none'}: ${problem.summary}`);
