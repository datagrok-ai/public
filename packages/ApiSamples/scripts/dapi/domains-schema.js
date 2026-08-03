// Schema lifecycle handle for a domain schema: manifest export, whole-schema audit,
// and schema-entity grants. (createSchema / apply / delete manage user-owned schemas.)

const schema = grok.dapi.domains.schema('apitests');

const manifest = await schema.manifest();
grok.shell.info(`tables: ${Object.keys(manifest.tables).join(', ')}`);

const audit = await schema.audit({limit: 5});
grok.shell.info(`latest events: ${audit.map((a) => a.op).join(', ')}`);

// Schema grants gate schema-level OPERATIONS (apply requires Edit, delete requires
// Delete, sharing requires Share) — they do NOT grant access to row data; grant per
// table for that: grok.dapi.domains.table('apitests.item').grant(groupId, 'View').
const grants = await schema.grants().catch(() => null); // requires Share on the schema
grok.shell.info(grants == null ? 'no Share on this schema' :
  (grants.map((g) => `${g.group.friendlyName}: ${g.permission}`).join(', ') || 'no direct grants'));

// Column security lives on the table client: shareColumn(column, group) restricts a
// column to its own per-column schema and grants View; restoreColumnVisibility(column)
// undoes it. Not run here — restricting columns affects every consumer of the schema.
