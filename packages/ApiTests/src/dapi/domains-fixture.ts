import type * as _grok from 'datagrok-api/grok';
declare let grok: typeof _grok;

// Shared fixture for the domain tests: the 'apitests' schema that this package
// declares in databases/apitests/schema.json (deployed on publish).

/** Typed client for the `apitests.item` table. */
export const items = () => grok.dapi.domains.table('apitests.item');

/** Unique business key for test rows: `<prefix>-<millis>-<random>`. */
export const uniqueKey = (prefix: string) => `${prefix}-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

/** Deletes all items whose sku starts with [prefix]. */
export async function purge(prefix: string): Promise<void> {
  const rows = await items().query({filter: `sku starts "${prefix}"`, limit: 1000});
  for (const r of rows)
    await items().delete(r.id);
}
