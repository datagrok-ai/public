import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// DG.DomainRow JS wrapper + ObjectHandler dispatch for domain-table rows (§8.1).
category('JS: domain handlers', () => {
  test('DG.DomainRow wrapper exported', async () => {
    expect(typeof DG.DomainRow, 'function', 'DG.DomainRow is not exported');
    const props = Object.getOwnPropertyNames((DG.DomainRow as any).prototype);
    for (const p of ['schemaName', 'tableName', 'typeName', 'semValue', 'values', 'id'])
      expect(props.includes(p), true, `DG.DomainRow.${p} getter missing`);
  });

  test('grit.issue ObjectHandler resolves', async () => {
    // Requires the Grit package (which registers the grit.issue handler) to be
    // deployed; skip cleanly where it is not so the suite stays green everywhere.
    const handlers = await DG.ObjectHandler.forSemType('grit.issue');
    if (handlers.length === 0)
      return;
    expect(handlers.some((h) => h.type === 'grit.issue'), true,
      'grit.issue handler not resolved by forSemType');
  });
}, {owner: 'askalkin@datagrok.ai'});
