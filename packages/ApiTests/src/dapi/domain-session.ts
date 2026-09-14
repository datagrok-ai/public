import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// The unit of work: DG.DomainSession saves several DomainFrameEditors as ONE
// /transaction; a row that does not exist yet carries a draft id ('~new:…') the
// other editor's rows may reference. Fixture: apitests.item (parent) and
// apitests.item_event (child, item_id ref). Every test cleans its prefix up.
category('Dapi: domain session', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const events = () => grok.dapi.domains.table('apitests.item_event');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const like = (property: string, prefix: string): any =>
    ({property, operator: 'like', value: `${prefix}%`});

  type Editor = _DG.DomainFrameEditor;

  const itemsEditor = (prefix: string): Promise<Editor> =>
    DG.DomainFrameEditor.create(items() as any, {query: {filter: like('sku', prefix), sort: 'sku'}, quiet: true});
  const eventsEditor = (prefix: string): Promise<Editor> =>
    DG.DomainFrameEditor.create(events() as any, {query: {filter: like('kind', prefix), sort: 'kind'}, quiet: true});

  async function cleanup(prefix: string): Promise<void> {
    try {
      await events().deleteWhere(like('kind', prefix));
      await items().deleteWhere(like('sku', prefix));
    } catch (e) {
      console.error(`session fixture ${prefix} not cleaned up: ${e}`);
    }
  }

  test('draft ids: addRow stamps ~new:, buildOps emits ref and $-rewrites a sibling draft reference', async () => {
    const prefix = `ds-draft-${stamp()}`;
    const parent = await itemsEditor(prefix);
    const child = await eventsEditor(prefix);
    try {
      const row = parent.addRow({sku: `${prefix}-0`, name: 'Draft item'});
      const draft = parent.dataFrame.get('id', row);
      expect(DG.DomainFrameEditor.isDraftId(draft), true, `addRow did not stamp a draft id: ${draft}`);
      child.addRow({item_id: draft, kind: `${prefix}-ev`, amount: 1});
      const [insert] = parent.buildOps();
      expect(insert.op.ref, draft, 'the insert op does not name its draft id as ref');
      expect('id' in (insert.op.values as any), false, 'the draft id leaked into the insert values');
      const [event] = child.buildOps();
      expect((event.op.values as any).item_id, `$${draft}`, 'the sibling draft reference was not rewritten to $ref');
    } finally {
      parent.detach();
      child.detach();
    }
  });

  test('two editors, one transaction: a draft item referenced by a draft event', async () => {
    const prefix = `ds-tx-${stamp()}`;
    const parent = await itemsEditor(prefix);
    const child = await eventsEditor(prefix);
    try {
      const row = parent.addRow({sku: `${prefix}-0`, name: 'Draft item'});
      const draft = parent.dataFrame.get('id', row);
      child.addRow({item_id: draft, kind: `${prefix}-ev`, amount: 1});
      const session = new DG.DomainSession([parent, child]);
      expect(session.changeCount, 2, 'the session does not aggregate changeCount');
      let saved: _DG.DomainSaveResult | null = null;
      session.onSaved.subscribe((r: _DG.DomainSaveResult) => saved = r);
      expect(await session.save(), true, 'the session save did not land');
      expect(parent.isDirty || child.isDirty, false, 'an editor stayed dirty after the session save');
      const itemId = parent.dataFrame.get('id', 0);
      expect(DG.DomainFrameEditor.isDraftId(itemId), false, 'the draft id was not replaced by the server id');
      expect(saved!.inserted, 2, `the session result does not sum the slices: ${JSON.stringify(saved)}`);
      expect(saved!.assigned[draft], itemId, 'assigned does not map the draft id to the server id');
      expect(child.dataFrame.get('item_id', 0), itemId, 'the event does not point at the inserted item');
      const server = await events().get(child.dataFrame.get('id', 0));
      expect(server.item_id, itemId, 'the server resolved the $ref to another id');
      const audit = await items().audit(itemId);
      const eventAudit = await events().audit(server.id);
      expect(audit[0]?.tx_id != null && audit[0].tx_id === eventAudit[0]?.tx_id, true,
        `the two inserts do not share one tx_id: ${audit[0]?.tx_id} vs ${eventAudit[0]?.tx_id}`);
    } finally {
      parent.detach();
      child.detach();
      await cleanup(prefix);
    }
  });

  test('a failing post-save re-read still resolves the draft reference in the child frame', async () => {
    const prefix = `ds-ref-${stamp()}`;
    const parent = await itemsEditor(prefix);
    const child = await eventsEditor(prefix);
    try {
      const row = parent.addRow({sku: `${prefix}-0`, name: 'Draft item'});
      const draft = parent.dataFrame.get('id', row);
      child.addRow({item_id: draft, kind: `${prefix}-ev`, amount: 1});
      // The post-save enrichment of the child: it must not be what resolves the ref.
      (child.client as any).query = async () => { throw new Error('post-save re-read refused'); };
      const session = new DG.DomainSession([parent, child]);
      expect(await session.save(), true, 'the session save did not land');
      const itemId = parent.dataFrame.get('id', 0);
      expect(DG.DomainFrameEditor.isDraftId(child.dataFrame.get('item_id', 0)), false,
        'the draft reference survived the save in the child frame');
      expect(child.dataFrame.get('item_id', 0), itemId, 'the child does not point at the inserted item');
      expect(parent.isDirty || child.isDirty, false, 'a failed re-read left an editor dirty');
      const server = await events().get(child.dataFrame.get('id', 0));
      expect(server.item_id, itemId, 'the transaction did not commit the resolved reference');
    } finally {
      delete (child.client as any).query;
      parent.detach();
      child.detach();
      await cleanup(prefix);
    }
  });

  test('a 409 in the child retries the whole batch after reload; dismiss leaves both editors pending', async () => {
    const prefix = `ds-409-${stamp()}`;
    const [item] = await items().insert({sku: `${prefix}-0`, name: 'Item'});
    const [event] = await events().insert({item_id: item.id, kind: `${prefix}-ev`, amount: 1});
    const original = (DG.DomainObjectHandler as any).showConflictDialog;
    let decision: 'reload' | null = 'reload';
    let asked = 0;
    (DG.DomainObjectHandler as any).showConflictDialog = async () => {
      asked++;
      return decision;
    };
    try {
      let parent = await itemsEditor(prefix);
      let child = await eventsEditor(prefix);
      parent.setValue(0, 'name', 'Renamed');
      child.setValue(0, 'amount', 2);
      await events().update(event.id, {amount: 9});
      let session = new DG.DomainSession([parent, child]);
      expect(await session.save(), true, 'the reload path did not finish the save');
      expect(asked, 1, 'the conflict dialog was not consulted once');
      expect(parent.isDirty || child.isDirty, false, 'the retried batch left edits pending');
      expect((await items().get(item.id)).name, 'Renamed', 'the parent edit did not land on the retry');
      expect((await events().get(event.id)).amount, 9, 'reload did not keep the server value');
      parent.detach();
      child.detach();

      decision = null;
      asked = 0;
      parent = await itemsEditor(prefix);
      child = await eventsEditor(prefix);
      parent.setValue(0, 'name', 'Never saved');
      child.setValue(0, 'amount', 3);
      await events().update(event.id, {amount: 10});
      session = new DG.DomainSession([parent, child]);
      expect(await session.save(), false, 'a dismissed conflict reported success');
      expect(child.errorOf(0, 'amount')?.kind, 'conflict', 'the dismissed conflict left no marker');
      expect(parent.isDirty && child.isDirty, true, 'a dismissed conflict dropped an edit');
      expect((await items().get(item.id)).name, 'Renamed', 'the parent edit landed despite the dismissal');
      parent.detach();
      child.detach();
    } finally {
      (DG.DomainObjectHandler as any).showConflictDialog = original;
      await cleanup(prefix);
    }
  });

  test('a validation error in the second editor lands on its cell and nothing landed for the first', async () => {
    const prefix = `ds-invalid-${stamp()}`;
    const [item] = await items().insert({sku: `${prefix}-0`, name: 'Item'});
    await events().insert({item_id: item.id, kind: `${prefix}-ev`, amount: 1});
    const child = await eventsEditor(prefix);
    const parent = await itemsEditor(prefix);
    try {
      child.setValue(0, 'amount', 5);
      // A duplicate sku is unique server-side only: the client validator lets it through.
      const row = parent.addRow({sku: `${prefix}-0`, name: 'Duplicate'});
      expect(parent.errorCount, 0, 'the client rejected the duplicate before the server saw it');
      const session = new DG.DomainSession([child, parent]);
      expect(await session.save(), false, 'a duplicate sku was accepted');
      expect(parent.errorOf(row, 'sku')?.kind, 'error', 'the server error did not reach the second editor');
      expect(child.isDirty, true, 'the first editor lost its edit');
      expect((await items().get(item.id)).name, 'Item', 'the transaction did not roll back');
    } finally {
      parent.detach();
      child.detach();
      await cleanup(prefix);
    }
  });

  test('confirmDiscardChanges(session) saves once (one transaction)', async () => {
    const prefix = `ds-gate-${stamp()}`;
    const [item] = await items().insert({sku: `${prefix}-0`, name: 'Item'});
    const parent = await itemsEditor(prefix);
    const child = await eventsEditor(prefix);
    try {
      parent.setValue(0, 'name', 'Gated');
      child.addRow({item_id: item.id, kind: `${prefix}-ev`, amount: 1});
      const session = new DG.DomainSession([parent, child]);
      let saves = 0;
      session.onSaved.subscribe(() => saves++);
      const gate = DG.confirmDiscardChanges(session, {action: 'leave'});
      await DG.delay(300);
      const save = Array.from(document.querySelectorAll('.d4-dialog button'))
        .find((b) => b.textContent?.trim() === 'SAVE') as HTMLElement | undefined;
      expect(save != null, true, 'the unsaved-changes dialog did not open');
      save!.click();
      expect(await gate, true, 'the gate did not let the caller proceed after the save');
      expect(saves, 1, 'the session did not save exactly once');
      expect(session.isDirty, false, 'the gate left the session dirty');
      // The trail is oldest first, and the item was inserted before the gate:
      // its LAST entry is the update the gate save wrote.
      const last = (trail: any[]): any => trail[trail.length - 1];
      const audit = last(await items().audit(item.id));
      const eventAudit = last(await events().audit(child.dataFrame.get('id', 0)));
      expect(audit?.tx_id != null && audit.tx_id === eventAudit?.tx_id, true,
        `the gate save was not one transaction: ${audit?.tx_id} vs ${eventAudit?.tx_id}`);
    } finally {
      parent.detach();
      child.detach();
      await cleanup(prefix);
    }
  });

  test('editor.save() is a session of one — same counters as before', async () => {
    const prefix = `ds-one-${stamp()}`;
    await items().insert({sku: `${prefix}-0`, name: 'Item'});
    const editor = await itemsEditor(prefix);
    try {
      let result: _DG.DomainSaveResult | null = null;
      editor.onSaved.subscribe((r) => result = r);
      editor.setValue(0, 'name', 'Alone');
      expect(await editor.save(), true, 'the lone save did not land');
      expect(JSON.stringify(result), JSON.stringify({inserted: 0, updated: 1, deleted: 0, assigned: {}}),
        'the lone save reports other counters');
    } finally {
      editor.detach();
      await cleanup(prefix);
    }
  });

  test('a literal leading $ round-trips', async () => {
    const prefix = `ds-dollar-${stamp()}`;
    const [item] = await items().insert({sku: `${prefix}-0`, name: 'Item'});
    const editor = await itemsEditor(prefix);
    try {
      editor.setValue(0, 'name', '$100');
      expect((editor.buildOps()[0].op.values as any).name, '$$100', 'a leading $ was not escaped');
      expect(await editor.save(), true, 'the save did not land');
      expect((await items().get(item.id)).name, '$100', 'the literal $ did not round-trip');
    } finally {
      editor.detach();
      await cleanup(prefix);
    }
  });

  test('a draft deleted before save produces no op; unmarkDeleted restores new', async () => {
    const prefix = `ds-drop-${stamp()}`;
    const editor = await itemsEditor(prefix);
    try {
      const row = editor.addRow({sku: `${prefix}-0`, name: 'Draft'});
      editor.markDeleted(row);
      expect(editor.buildOps().length, 0, 'a deleted draft produced an op');
      editor.unmarkDeleted(row);
      expect(editor.stateOf(row), 'new', 'unmarkDeleted did not restore the draft to new');
      editor.markDeleted(row);
      expect(editor.prepareSave()?.length, 0, 'prepareSave did not resolve the deleted draft locally');
      expect(editor.dataFrame.rowCount, 0, 'the deleted draft survived prepareSave');
      expect(editor.isDirty, false, 'the resolved draft held the editor dirty');
    } finally {
      editor.detach();
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
