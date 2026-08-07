// UI-FACADE §1.4 — a list page scoped by a query, in ONE line; and the reference
// picker, in one more.
//
// Free: everything the list page of §1.1 has (cards / brief / editable grid,
// search over the identity columns, New, per-item commands, the row cap and its
// banner, the unsaved-changes gate), scoped by the query and deep-linkable as
// opened.
//
// In a package: `import {domains} from '@datagrok-libraries/domain-ui';`
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

// The example, verbatim (modulo the fixture table and its columns):
grok.shell.addView((await domains.table('apitests.item')).listView({query: {filter: 'quantity > 0'}}));

// The reference picker, anywhere — a one-shot dialog, so it stays
// string-addressable: the target table's full Domain View in a dialog.
const item = await domains.pick('apitests.item');   // DomainRow | null
grok.shell.info(item == null ? 'Nothing picked' : `Picked ${item.displayName}`);
