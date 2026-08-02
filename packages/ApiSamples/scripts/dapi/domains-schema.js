// Schema lifecycle handle for a domain schema: manifest export and whole-schema audit.

const schema = grok.dapi.domains.schema('apitests');

const manifest = await schema.manifest();
grok.shell.info(`tables: ${Object.keys(manifest.tables).join(', ')}`);

const audit = await schema.audit({limit: 5});
grok.shell.info(`latest events: ${audit.map((a) => a.op).join(', ')}`);
