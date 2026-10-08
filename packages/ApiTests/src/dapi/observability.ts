import type * as _grok from 'datagrok-api/grok';
declare let grok: typeof _grok;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi: observability', () => {
  test('errors: occurrences', async () => {
    const rows = await grok.dapi.log.getErrors({since: '1h', limit: 5});
    expect(Array.isArray(rows), true);
    expect(rows.length <= 5, true);
    for (const r of rows)
      expect(typeof r.time, 'string');
  }, {owner: 'aparamonov@datagrok.ai'});

  test('errors: figures by signature', async () => {
    const rows = await grok.dapi.log.getErrors({since: '1d', by: 'signature', limit: 3});
    expect(rows.length <= 3, true);
    for (const r of rows) {
      expect(typeof r.count, 'number');
      expect(typeof r.signature, 'string');
    }
  }, {owner: 'aparamonov@datagrok.ai'});

  test('logging: policy', async () => {
    const policy = await grok.dapi.log.getLoggingPolicy();
    expect(policy.settings != null, true);
    expect(Array.isArray(policy.debugFlags), true);
  }, {owner: 'aparamonov@datagrok.ai'});

  test('logging: timeline of the current session', async () => {
    const session = await grok.dapi.users.currentSession();
    const rows = await grok.dapi.log.getTimeline({session: session.id});
    expect(Array.isArray(rows), true);
    for (const r of rows)
      expect(typeof r.time, 'string');
  }, {owner: 'aparamonov@datagrok.ai'});

  test('logging: add and stop a capture rule', async () => {
    const rule = await grok.dapi.log.addCaptureRule({subject: {type: 'user', value: grok.shell.user.id},
      capture: {requests: true}, forMinutes: 5, reason: 'api test'});
    let stopped = false;
    try {
      expect(rule.status, 'active');
      expect(rule.capture.requests, true);
      expect((await grok.dapi.log.stopCaptureRule(rule.id, 'api test')).status, 'stopped');
      stopped = true;
    } finally {
      if (!stopped)
        await grok.dapi.log.stopCaptureRule(rule.id, 'api test').catch(() => {});
    }
  }, {owner: 'aparamonov@datagrok.ai'});
});
