import * as grok from 'datagrok-api/grok';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

/** Northwind on domain databases: entity-mapped tables shipped with this package
 * (databases/northwind/schema.json + seed data), browsed through the standard
 * domain UI (/domains routing) and queried through grok.dapi.domains. */
export async function northwindDemo(): Promise<void> {
  const script = new DemoScript('Domain Databases',
    'Northwind as a domain database: tables shipped by a plugin become platform entities ' +
    'with security, audit history, and standard browsing UI — no backend code',
    false, {autoStartFirstStep: true, path: 'Data Access/Domain Databases'});

  await script
    .step('Enable domain databases', async () => {
      if ((grok.shell.settings as any).enableDomainDatabases !== true) {
        (grok.shell.settings as any).enableDomainDatabases = true;
        grok.shell.info('Enabled the "Domain databases" beta feature (Settings > Beta). ' +
          'After a page refresh, it also appears under Browse > Platform > Domains.');
      }
    }, {
      description: 'Domain databases are PostgreSQL schemas declared by a plugin manifest ' +
        '(databases/northwind/schema.json in the Tutorials package). The platform creates the tables, ' +
        'registers every table as a securable entity, and provides CRUD with row/column security and ' +
        'an audit trail. This step enables the beta feature flag.',
      delay: 2000,
    })
    .step('Browse domains', async () => {
      grok.shell.route('/domains');
    }, {
      description: 'The Domains gallery lists every registered domain schema. ' +
        'The "northwind" schema ships with the Tutorials package: eleven tables with ' +
        'foreign keys, business keys, and referential actions — deployed automatically on package install.',
      delay: 2000,
    })
    .step('Northwind tables', async () => {
      grok.shell.route('/domains/northwind');
    }, {
      description: 'Tables of the northwind schema. Products reference suppliers and categories; ' +
        'orders reference customers, employees, and shippers. Open the schema diagram from the ribbon ' +
        'to see the relationships, including the security delegation from order details to orders.',
      delay: 2000,
    })
    .step('Browse orders', async () => {
      grok.shell.route('/domains/northwind/orders');
    }, {
      description: 'The Domain View shows the table rows with search, facet filters, and ' +
        'card / brief / grid render modes. Orders use row-level security: an order can be shared ' +
        'with a user or group individually, while order details inherit access from their order.',
      delay: 2000,
    })
    .step('Open an order', async () => {
      grok.shell.route('/domains/northwind/orders/10248');
    }, {
      description: 'Each row is addressable by its business key — this is order #10248. ' +
        'The entity view shows the order fields with referenced entities resolved to display names, ' +
        'a tab with its order details (line items), and a History tab where every change made ' +
        'through the platform is recorded as an audit trail.',
      delay: 2000,
    })
    .step('Query from code', async () => {
      const df = await grok.dapi.domains.table('northwind.orders')
        .queryDf({filter: 'ship_country = "Germany"', sort: '!order_date', limit: 100});
      df.name = 'German orders';
      grok.shell.addTableView(df);
    }, {
      description: 'The same data is one call away in JavaScript: ' +
        "grok.dapi.domains.table('northwind.orders').queryDf({filter: 'ship_country = \"Germany\"'}) — " +
        'with row and column security applied server-side. Typed per-table clients can be generated ' +
        'from the manifest with "grok api".',
      delay: 2000,
    })
    .start();
}
