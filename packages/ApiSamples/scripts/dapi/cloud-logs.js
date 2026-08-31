// Reading the cloud log tiers: CloudWatch (hot, 30 days) and the write-once S3
// archive. Admins only; see core/docs/LOG_ARCHITECTURE.md.

// Log groups the instance can reach. Pass {connection: <S3 connection id>} to use
// stored credentials instead of the instance's own AWS role.
let groups = await grok.dapi.log.getCloudLogGroups({prefix: '/datagrok/'});
grok.shell.info(`${groups.length} log group(s): ${groups.join(', ')}`);

// Events from one group. `filter` is a CloudWatch filter pattern, not a regex.
let end = dayjs();
let events = await grok.dapi.log.getCloudLogEvents(groups[0], end.subtract(1, 'hour'), end,
  {filter: 'error', limit: 500});
grok.shell.addTableView(events);

// The archive is two-step: list keys under a prefix, then decode the one you want.
let connection = await grok.dapi.connections.filter('LogArchive').first();
if (connection != null) {
  let objects = await grok.dapi.log.getArchiveObjects(connection.id, {prefix: 'cloudwatch/', limit: 50});
  grok.shell.addTableView(objects);

  if (objects.rowCount > 0)
    grok.shell.addTableView(await grok.dapi.log.getArchiveEvents(connection.id, objects.get('key', 0)));
}
