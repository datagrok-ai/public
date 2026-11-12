import {line} from '../logger';

describe('line', () => {
  test('GROK_JSON_LOGS=true prints one JSON line in the server shape', () => {
    const out = line('WARN', 'Send failed for wss://bot:pw-123456@pipe:3000 with Bearer abc.def-0123456789', 't-1',
      {GROK_JSON_LOGS: 'true', DATAGROK_CELERY_NAME: 'chem-celery'});
    expect(out).not.toContain('\n');
    const json = JSON.parse(out);
    expect(Object.keys(json).sort()).toEqual(['level', 'message', 'params', 'service', 'time', 'v']);
    expect(json).toMatchObject({v: 1, level: 'warning', service: 'chem-celery', params: {taskId: 't-1'}});
    expect(new Date(json.time).toISOString()).toBe(json.time);
    expect(out).not.toContain('pw-123456');
    expect(out).not.toContain('abc.def-0123456789');
    expect(out).toContain('pipe:3000');
  });

  test('text by default', () => {
    expect(line('INFO', 'Started', undefined, {})).toMatch(/^\S+ \[INFO\] Started$/);
  });
});
