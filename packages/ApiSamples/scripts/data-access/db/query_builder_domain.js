// Visual (structured) query over a DOMAIN schema — entity-mapped tables managed by the
// platform. Its connection is credential-less, absent from every listing, and addressed
// by the hash of 'System:Domain.<schema>'; per-user row and column security applies.

const schema = 'apitests'; // packages/ApiTests/databases/apitests/schema.json

// The JS equivalent of the platform's entity-id hash (uuid v5, SHA-1).
async function hashId(name) {
  const ns = '0bb1a16edf0745f599e08e67000ecc22';
  const encoded = new TextEncoder().encode(name);
  const buf = new Uint8Array(16 + encoded.length);
  for (let i = 0; i < 16; i++)
    buf[i] = parseInt(ns.substring(i * 2, i * 2 + 2), 16);
  buf.set(encoded, 16);
  const d = new Uint8Array(await crypto.subtle.digest('SHA-1', buf)).subarray(0, 16);
  d[6] = (d[6] & 0x0f) | 0x50;
  d[8] = (d[8] & 0x3f) | 0x80;
  const s = Array.from(d).map((b) => b.toString(16).padStart(2, '0')).join('');
  return `${s.slice(0, 8)}-${s.slice(8, 12)}-${s.slice(12, 16)}-${s.slice(16, 20)}-${s.slice(20)}`;
}

try {
  let conn = await grok.dapi.connections.find(await hashId(`System:Domain.${schema}`));
  if (conn == null)
    throw new Error(`Domain schema '${schema}' is not deployed on this server`);

  let df = await conn.buildQuery('item')
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
