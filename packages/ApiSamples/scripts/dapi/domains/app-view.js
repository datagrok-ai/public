// The domain app framework (@datagrok-libraries/domain-ui): a browse/CRUD app
// over a domain table is ONE expression — list page, row page, deep links,
// permission gating and the unsaved-changes prompt included.
//
// This is the CLASS layer, which stays public as the extension point. The way in
// for new code is the handle — see facade-app.js for the same app in one
// line.
//
// In a package you write `import {DomainAppView} from '@datagrok-libraries/domain-ui'`.
// A sample script has no bundler, so it borrows the copy ApiTests bundles and
// re-exports; calling any of its functions loads the package.
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {DomainAppView, DomainEntityAppView, domainHandler} = apitests.domainUi;

const items = grok.dapi.domains.table('apitests.item');

// 1. The list page: cards / brief / grid, search over the identity columns, New,
// per-row commands, and the query round-tripped through the URL.
const list = new DomainAppView(items, {name: 'Items'});
grok.shell.addView(list);

// 2. The page of ONE row: its property form, a tab per child table (editable
// grids with the parent FK prefilled), history, and the actions this user may
// perform. `rowFrom` wraps a row the query already returned — no round trip.
const first = (await items.query({limit: 1}))[0];
if (first == null)
  grok.shell.info('apitests.item has no rows yet — the list page offers "New Item...".');
else
  grok.shell.addView(DomainEntityAppView.forRow(items, domainHandler('apitests.item').rowFrom(first)));
