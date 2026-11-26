//api: DG.Model.inputs, DG.Model.outputs
// Predictive models and the columns they take and return

const models = await grok.dapi.models.list({pageSize: 10});
for (const model of models) {
  const columns = (list) => list.map((c) => `${c.name} (${c.type})`).join(', ');
  grok.shell.info(`${model.friendlyName}: ${columns(model.inputs)} -> ${columns(model.outputs)}`);
}
