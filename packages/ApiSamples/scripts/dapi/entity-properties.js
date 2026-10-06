// Discover what the platform knows about one of its own types, then filter by it.
// The names grok.meta returns are the names the platform's grids, filter panels and
// server-side facets use — so discovery and enforcement never drift apart.

// The curated catalog, in the order the platform renders it.
const props = await grok.meta.propertiesOf('User');
grok.shell.info(props.map((p) => `${p.friendlyName} (${p.name}: ${p.type})`).join('\n'));

// The database-backed subset — the properties a filter string may name.
const filterable = await grok.meta.propertiesOf('User', {filterable: true});
const me = await grok.dapi.users.current();
const by = filterable[0].name;
const sameUser = await grok.dapi.users.filter(`${by} = "${me[by]}"`).list();
grok.shell.info(`filterable: ${filterable.map((p) => p.name).join(', ')} — ` +
  `${by} resolves ${sameUser.length} user`);

// References say what they point at, so a related-entity filter can be built blind.
const sessionProps = await grok.meta.propertiesOf('UserSession');
const ref = sessionProps.find((p) => p.refType === 'User');
grok.shell.info(`UserSession.${ref.name} ${ref.relationKind} ${ref.refType}`);

// Where those rows are served as a read-only domain table — facets, hops and
// filters over the SAME names, no codegen. Null for a type the platform does not serve.
const core = await grok.meta.coreLocationOf('UserSession');
if (core === null)
  grok.shell.info('UserSession rows are not served as a domain table here');
else {
  const sessions = grok.dapi.domains.table(`${core.schema}.${core.table}`);
  const mine = await sessions.count({property: `${ref.name}.${by}`, operator: '=', value: me[by]});
  grok.shell.info(`${core.schema}.${core.table}: ${mine} of ${await sessions.count()} sessions are mine`);
}
