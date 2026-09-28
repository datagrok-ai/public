import {describe, it, expect} from 'vitest';
import {parseArgs} from '../utils/cli-args';

describe('parseArgs', () => {
  it('keeps grok s positionals verbatim', () => {
    expect(parseArgs(['s', 'alerts', 'ack', '08327370', '--host', 'obs-b'])._).toEqual(['s', 'alerts', 'ack', '08327370']);
    expect(parseArgs(['s', 'alerts', 'get', '12345678'])._).toEqual(['s', 'alerts', 'get', '12345678']);
    expect(parseArgs(['server', 'errors', 'show', '012345'])._).toEqual(['server', 'errors', 'show', '012345']);
    expect(parseArgs(['s', 'capture', 'show', '17'])._[3]).toBe('17');
  });

  it('keeps id options verbatim', () => {
    const argv = parseArgs(['s', 'timeline', '--action', '0123456789', '--rule', '017', '--until-version', '1.10']);
    expect([argv.action, argv.rule, argv['until-version']]).toEqual(['0123456789', '017', '1.10']);
  });

  it('keeps a dash-leading value with its option, for grok s only', () => {
    const argv = parseArgs(['s', 'logger', 'set', '--save-levels', '-audit', '--since', '-7d']);
    expect([argv['save-levels'], argv.since]).toEqual(['-audit', '-7d']);
    expect(parseArgs(['publish', '--from', '-x']).from).toBe(true);
  });

  it('leaves other commands as they were', () => {
    expect(parseArgs(['test', '42'])._).toEqual(['test', 42]);
  });
});
