import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';

// GROK-20572 (domain-visual-queries) — TableQueryBuilder over a *domain* connection.
// Each registered domain schema gets one credential-less, isSystem DataConnection
// (dataSource 'Domain'); a structured TableQuery on it is compiled and executed
// inside Datlas by DomainDataProvider through DomainRepository, so the per-user row
// predicate and column security apply to visual queries exactly as they do to the
// /domains REST API. Raw SQL on such a connection is refused server-side.
//
// The fixture is this package's own 'apitests' domain schema
// (databases/apitests/schema.json), deployed on publish. Everything self-skips when
// the schema (or the virtual connection — i.e. an older Datlas) is not there.
category('Dapi: domain visual queries', () => {
  const schema = 'apitests';
  let conn: _DG.DataConnection | null = null;
  let skipReason = `the '${schema}' domain schema or its virtual connection is not deployed ` +
    '(older Datlas without domain-visual-queries)';
  const scratchSkus: string[] = [];

  const items = () => grok.dapi.domains.table(`${schema}.item`);
  const sku = () => `VQ-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  /** RFC 4122 v5 (SHA-1) UUID over the platform's entity-id namespace — the JS
   * equivalent of `NamedModel.getHashId`. There is no JS API for it yet, so the
   * recipe is spelled out here; see the report for the follow-up. */
  async function hashId(name: string): Promise<string> {
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
    return `${s.substring(0, 8)}-${s.substring(8, 12)}-${s.substring(12, 16)}-` +
      `${s.substring(16, 20)}-${s.substring(20, 32)}`;
  }

  before(async () => {
    if (crypto.subtle == null) {
      skipReason = 'crypto.subtle is unavailable (a non-secure origin), so the connection id ' +
        'cannot be computed';
      return;
    }
    // Only the presence check is guarded: an older Datlas has neither the domain
    // schema nor the virtual connection, and that is the one case worth skipping.
    try {
      const schemas = await grok.dapi.domains.schemas.list();
      if (!schemas.some((s) => s.name === schema)) {
        skipReason = `the '${schema}' domain schema is not deployed`;
        return;
      }
      // Domain connections are isSystem and therefore absent from every connection
      // listing; they are addressed by their deterministic id, hash of
      // 'System:Domain.<schema>' (DataConnection.domainConnectionId).
      conn = await grok.dapi.connections.find(await hashId(`System:Domain.${schema}`)) ?? null;
    } catch (e: any) {
      conn = null;
      skipReason = `${skipReason}: ${e.message ?? e}`;
      return;
    }
    if (conn == null)
      return;
    // Past the guard: a failure here is a real regression, not a reason to skip.
    for (const s of [sku(), sku()])
      scratchSkus.push(s);
    await items().insert([
      {sku: scratchSkus[0], name: 'VQ Alpha', quantity: 7},
      {sku: scratchSkus[1], name: 'VQ Beta', quantity: 3},
    ]);
  });

  after(async () => {
    for (const s of scratchSkus) {
      try {
        const rows = await items().query({filter: `sku = "${s}"`});
        for (const r of rows)
          await items().delete(r.id);
      } catch (e: any) {
        console.log(`cleanup of '${s}' failed: ${e.message ?? e}`);
      }
    }
  });

  test('buildQuery: select, where, sortBy, limit', async () => {
    if (conn == null) { console.log(`skipped: ${skipReason}`); return; }
    const df = await DG.TableQueryBuilder.from('item', conn)
      .select(['sku', 'name', 'quantity'])
      .where('sku', `in (${scratchSkus.join(', ')})`, DG.TYPE.STRING)
      .sortBy('quantity', false)
      .limit(10)
      .build()
      .executeTable();
    // Exactly the requested fields — no system columns, no hidden 'note'.
    expect(df.columns.names().sort().join(','), 'name,quantity,sku');
    expect(df.rowCount, 2);
    expect(df.col('quantity')!.get(0), 7);
    expect(df.col('sku')!.get(0), scratchSkus[0]);
  });

  test('buildQuery: aggregation over a domain table', async () => {
    if (conn == null) { console.log(`skipped: ${skipReason}`); return; }
    const df = await DG.TableQueryBuilder.from('item', conn)
      .where('sku', `in (${scratchSkus.join(', ')})`, DG.TYPE.STRING)
      .groupBy(['name'])
      .sum('quantity', 'total')
      .sortBy('name', true)
      .build()
      .executeTable();
    expect(df.rowCount, 2);
    expect(df.col('name')!.get(0), 'VQ Alpha');
    expect(df.col('total')!.get(0), 7);
    expect(df.col('total')!.get(1), 3);
  });

  test('raw SQL on a domain connection is refused', async () => {
    if (conn == null) { console.log(`skipped: ${skipReason}`); return; }
    let message = '';
    try {
      await conn.query('domain adhoc', 'select 1').apply();
    } catch (e: any) {
      message = e.message ?? `${e}`;
    }
    expect(message.includes('Raw SQL is not supported on domain connections'), true,
      `expected the provider refusal, got: '${message}'`);
  });
}, {owner: 'askalkin@datagrok.ai'});
