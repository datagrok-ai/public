// `grok.dapi.entities.save` is polymorphic: it accepts any Entity subclass
// and dispatches client-side to the matching typed endpoint. Pass a Project
// and it lands the same way `grok.dapi.projects.save` would; pass a
// DataConnection and it goes through `grok.dapi.connections.save`.

const project = DG.Project.create();
project.name = 'sample_polymorphic_save_project';
const savedProject = await grok.dapi.entities.save(project);
grok.shell.info(`Saved Project as ${savedProject.id}`);

const dc = DG.DataConnection.create('sample_polymorphic_save_conn', {
  dataSource: 'PostgresDart', server: 'localhost:5432', db: 'datagrok_dev',
  login: 'datagrok_dev', password: '123',
});
const savedConn = await grok.dapi.entities.save(dc);
grok.shell.info(`Saved DataConnection as ${savedConn.id}`);

await grok.dapi.projects.delete(savedProject);
await grok.dapi.connections.delete(savedConn);
