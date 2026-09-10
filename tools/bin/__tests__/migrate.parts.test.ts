import {describe, it, expect} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {NodeDapi} from '../utils/node-dapi';
import {missingPackages, missingUsers, namespacesOf, readState, writeState} from '../utils/migrate/parts';

function makeDapi(responder: (path: string) => any) {
  const client: any = {
    baseUrl: 'http://h/api',
    async request(_method: string, p: string) { return responder(p); },
    get(p: string) { return this.request('GET', p); },
    post(p: string, body?: any) { return this.request('POST', p, body); },
  };
  return new NodeDapi(client);
}

describe('namespacesOf', () => {
  it('lists team spaces and personal roots, and leaves the platform\'s own projects out', async () => {
    const dapi = makeDapi(() => [
      {name: 'Reports', namespace: 'Chem:'},
      {name: 'Dashboard', namespace: 'Chem:Sub:'},
      {name: 'Skalkin', isRoot: true},
      {name: 'Chem', isRoot: true, isEntity: true, isPackage: true},
      {name: 'SavedFilter', isEntity: true},
    ]);
    expect(await namespacesOf(dapi)).toEqual(['Chem', 'Skalkin']);
  });

  it('does not repeat a namespace many projects share', async () => {
    const dapi = makeDapi(() => [
      {name: 'A', namespace: 'Team:'}, {name: 'B', namespace: 'Team:'}, {name: 'C', namespace: 'Team:Deep:'},
    ]);
    expect(await namespacesOf(dapi)).toEqual(['Team']);
  });
});

describe('prerequisites', () => {
  it('names the users the target does not have', async () => {
    const from = makeDapi(() => [{login: 'alice'}, {login: 'bob'}, {login: 'carol'}]);
    const to = makeDapi(() => [{login: 'bob'}]);
    expect(await missingUsers(from, to)).toEqual(['alice', 'carol']);
  });

  it('names the packages the target does not have', async () => {
    const from = makeDapi(() => [{name: 'Chem'}, {name: 'Bio'}]);
    const to = makeDapi(() => [{name: 'Chem'}]);
    expect(await missingPackages(from, to)).toEqual(['Bio']);
  });

  it('reports nothing when the target already has everything', async () => {
    const both = () => [{login: 'alice', name: 'Chem'}];
    expect(await missingUsers(makeDapi(both), makeDapi(both))).toEqual([]);
    expect(await missingPackages(makeDapi(both), makeDapi(both))).toEqual([]);
  });
});

describe('run state', () => {
  it('round-trips what finished so a re-run can skip it', () => {
    const file = path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'parts-')), 'nested', 'state.json');
    writeState(file, {Chem: {name: 'Chem', entities: 12, failed: 0, seconds: 3}});
    expect(readState(file).Chem).toMatchObject({entities: 12, failed: 0});
  });

  it('treats a missing or unreadable state file as a fresh run', () => {
    expect(readState(path.join(os.tmpdir(), 'grok-parts-does-not-exist.json'))).toEqual({});
  });
});
