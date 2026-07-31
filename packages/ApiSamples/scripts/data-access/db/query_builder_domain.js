// Visual (structured) query over a DOMAIN schema — entity-mapped tables managed
// by the platform (grok.dapi.domains), not an external database.
//
// Every registered domain schema gets one credential-less, system DataConnection
// (dataSource 'Domain'). It never appears in connection lists; it is addressed by
// its deterministic id, the platform hash of 'System:Domain.<schema>'. Datlas
// compiles the structured query and runs it through the domain repository, so the
// per-user row predicate and column security apply exactly as they do to the
// /domains REST API — and raw SQL on such a connection is refused.
//
// Requires a deployed domain schema; this sample uses 'apitests'
// (packages/ApiTests/databases/apitests/schema.json).

const schema = 'apitests';

// The JS equivalent of the platform's entity-id hash (uuid v5, SHA-1).
async function hashId(name) {
  const nsHex = '0bb1a16e-df07-45f5-99e0-8e67000ecc22'.replace(/-/g, '');
  const encoded = new TextEncoder().encode(name);
  const buf = new Uint8Array(16 + encoded.length);
  for (let i = 0; i < 16; i++)
    buf[i] = parseInt(nsHex.substring(i * 2, i * 2 + 2), 16);
  buf.set(encoded, 16);
  const digest = new Uint8Array(await crypto.subtle.digest('SHA-1', buf)).subarray(0, 16);
  digest[6] = (digest[6] & 0x0f) | 0x50;
  digest[8] = (digest[8] & 0x3f) | 0x80;
  const s = Array.from(digest).map((b) => b.toString(16).padStart(2, '0')).join('');
  return `${s.substring(0, 8)}-${s.substring(8, 12)}-${s.substring(12, 16)}-${s.substring(16, 20)}-${s.substring(20, 32)}`;
}

try {
  let conn = await grok.dapi.connections.find(await hashId(`System:Domain.${schema}`));
  if (conn == null)
    throw new Error(`Domain schema '${schema}' is not deployed on this server`);

  let df = await DG.TableQueryBuilder.from('item', conn)
    .select(['sku', 'name', 'quantity'])
    .where('quantity', '> 0', DG.TYPE.INT)
    .sortBy('sku', true)
    .limit(10)
    .build()
    .executeTable();
  grok.shell.addTableView(df);
} catch (e) {
  grok.shell.error(e);
}
