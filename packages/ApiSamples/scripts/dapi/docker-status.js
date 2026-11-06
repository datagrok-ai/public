// Whether the server can run containers: Grok Spawner is switched on and its last poll found no connection failure
let status = await grok.dapi.docker.getStatus();
if (status.isActive == null)
  grok.shell.warning('The Grok Spawner status could not be read');
else if (!status.enabled)
  grok.shell.info('Grok Spawner is switched off');
else if (!status.isActive)
  grok.shell.info('Grok Spawner is unreachable');
else
  grok.shell.info('Containers can run');
