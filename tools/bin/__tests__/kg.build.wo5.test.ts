/// `grok kg build` WO-5 (build-plan.md): the process layer over the mini monorepo under
/// fixtures/kg/build — tickets from the backlog snapshot, the release record with its picked
/// commits, and the people the two of them name.
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {spawnSync} from 'child_process';

import {currentDir} from '../utils/kg/build/write';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const RECORD = 'core/docs/release/1.0.1.yaml';

interface Built {
  repo: string;
  rows: (file: string) => any[];
  report: (name: string) => any;
  manifest: any;
}

async function build(extra: Record<string, unknown> = {}, prepare?: (repo: string) => void): Promise<Built> {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-wo5-'));
  fs.cpSync(fixture, repo, {recursive: true});
  prepare?.(repo);
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), backlog: path.join(repo, 'backlog'), only: 'homes,process', db: false, output: 'json', ...extra});
    const out = currentDir(path.join(repo, '.kg'))!;
    const rows = (file: string) => {
      const p = path.join(out, `data/${file}.jsonl`);
      return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
    };
    const report = (name: string) => {
      const p = path.join(out, 'reports', `${name}.json`);
      return fs.existsSync(p) ? JSON.parse(fs.readFileSync(p, 'utf8')) : undefined;
    };
    return {repo, rows, report, manifest: JSON.parse(String(log.mock.calls[0][0]))};
  } finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

/** Turns the fixture into two git repositories and rewrites the release record's picks with the shas they got. */
function commits(repo: string): void {
  const git = (cwd: string, ...args: string[]) => spawnSync('git', ['-c', 'user.name=t', '-c', 'user.email=t@x', ...args], {cwd, encoding: 'utf8'});
  const shas: string[] = [];
  for (const dir of [repo, path.join(repo, 'public')]) {
    git(dir, 'init', '-q');
    git(dir, 'add', '-A');
    git(dir, 'commit', '-q', '-m', 'fixture');
    shas.push(git(dir, 'rev-parse', '--short=10', 'HEAD').stdout.trim());
  }
  const record = path.join(repo, ...RECORD.split('/'));
  fs.writeFileSync(record, fs.readFileSync(record, 'utf8').replace(/aaaaaaaaaa/g, shas[0]).replace(/bbbbbbbbbb/g, shas[1]));
}

const graph = build();
const byId = (rows: any[], id: string) => rows.find((r) => r.id === id);
const edges = (rows: any[], from?: string, to?: string) => rows.filter((e) => (from === undefined || e.from === from) && (to === undefined || e.to === to));

describe('tickets from the backlog snapshot (build-plan.md WO-5)', () => {
  it('copies the tracker columns, maps status to state and the markdown status to raw_status', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/ticket'), 'GROK-100')).toMatchObject({
      type: 'ticket', name: 'Legend clipped in the scatter plot', tracker: 'jira', key: 'GROK-100', kind: 'bug', state: 'done', raw_status: 'Done',
      priority: 'high', components: ['Core'], labels: ['ui'], resolution: 'completed', closed: '2026-01-09T11:00:00Z', area: 'viewers',
      url: 'https://jira.example/GROK-100', created: '2026-01-02T10:00:00Z', updated: '2026-01-09T11:00:00Z',
      assignee: 'P:jane', reporter: 'P:rob', provenance: 'external', source_layer: 'process', visibility: 'dev',
    });
    expect(byId(rows('nodes/ticket'), 'GROK-101')).not.toHaveProperty('components');
  });

  it('writes a GitHub issue as gh:public#n with the # key and its own raw status', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/ticket'), 'gh:public#7')).toMatchObject({tracker: 'github', key: '#7', state: 'wontfix', raw_status: 'closed/not_planned',
      labels: ['Ins', 'enhancement'], resolution: 'not_planned'});
    expect(rows('nodes/ticket').map((t) => t.id)).toEqual(['GROK-100', 'GROK-101', 'GROK-102', 'gh:public#7']);
  });

  it('takes a numeric fix version to a release stub and leaves a bucket such as v1 alone', async () => {
    const {rows} = await graph;
    expect(edges(rows('edges/targets-release'), 'GROK-100')).toEqual([expect.objectContaining({
      to: 'Rel:1.0.1', kind: 'fix-version', derived_by: 'external', confidence: 1, evidence: ['backlog/index.jsonl'],
    })]);
    expect(rows('nodes/release').map((r) => r.id)).toEqual(['Rel:1.0.0', 'Rel:1.0.1']);
  });

  it('resolves a customer label through the backlog taxonomy into a Cust: stub and a requested-by edge', async () => {
    const {rows} = await graph;
    expect(rows('nodes/customer')).toEqual([
      expect.objectContaining({id: 'Cust:insitro', name: 'Insitro', jira_client: 'Insitro', relationship: 'customer', provenance: 'external', visibility: 'internal'}),
      expect.objectContaining({id: 'Cust:nx', name: 'NX', jira_client: 'NX'}),
    ]);
    expect(edges(rows('edges/requested-by'), 'gh:public#7')).toEqual([expect.objectContaining({to: 'Cust:insitro', derived_by: 'external', evidence: ['backlog/index.jsonl']})]);
  });

  it('reads the ticket body the snapshot has no column for: ~id as affects, ticket keys as counted mentions, never itself', async () => {
    const {rows} = await graph;
    expect(edges(rows('edges/affects'), 'GROK-100')).toEqual([expect.objectContaining({
      to: 'visualize/viewers/scatter-plot', derived_by: 'annotation', confidence: 1, evidence: ['backlog/jira/GROK-100.md'],
    })]);
    expect(edges(rows('edges/affects'), 'gh:public#7').map((e) => e.to)).toEqual(['domains/bio']);
    expect(edges(rows('edges/mentions'), 'GROK-100')).toEqual([expect.objectContaining({to: 'GROK-101', count: 2, evidence: ['backlog/jira/GROK-100.md']})]);
  });

  it('watermarks the backlog source with the newest update it saw, and reports it missing when the folder is not there', async () => {
    const {manifest} = await graph;
    expect(manifest.sources.backlog).toBe('ok@2026-01-12T07:00:00Z');
    const absent = await build({backlog: path.join(os.tmpdir(), 'grok-kg-no-backlog')});
    expect(absent.manifest.sources.backlog).toBe('missing');
    expect(absent.rows('nodes/ticket').map((t) => t.id)).toEqual(['GROK-101']);
  });
});

describe('people resolution (build-plan.md WO-5)', () => {
  it('resolves a handle a home declares, a roster identity a home matches, and the roster itself as a stub', async () => {
    const {rows} = await graph;
    expect(edges(rows('edges/assignee'), 'GROK-102')).toEqual([expect.objectContaining({type: 'ref', name: 'assignee', to: 'P:jane', derived_by: 'external'})]);
    expect(edges(rows('edges/assignee'), 'GROK-100').map((e) => e.to)).toEqual(['P:jane']);
    expect(byId(rows('nodes/person'), 'P:rob')).toMatchObject({name: 'Rob Ops', handle: 'rob', github: 'robops', emails: ['rob@example.com'],
      provenance: 'external', status: 'proposed', visibility: 'internal'});
  });

  it('matches a display name against the roster, and leaves a name nobody knows unresolved and listed', async () => {
    const {rows, report, manifest} = await graph;
    expect(edges(rows('edges/reporter'), 'GROK-100').map((e) => e.to)).toEqual(['P:rob']);
    expect(edges(rows('edges/assignee'), 'GROK-101')).toEqual([]);
    expect(manifest.sources.people).toBe('partial');
    expect(report('unresolved-people')).toEqual([{key: 'ghost', count: 1}]);
  });
});

describe('releases and their picks (build-plan.md WO-5)', () => {
  it('reads the record into a release with its branch and base, a dry run being a proposal', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/release'), 'Rel:1.0.1')).toMatchObject({
      type: 'release', name: '1.0.1', version: '1.0.1', kind: 'patch', state: 'testing', branch: 'release/1.0.1', base: 'Rel:1.0.0',
      record: RECORD, status: 'proposed', description: 'Dry run record generated 2026-01-15; not the live record.',
      provenance: 'annotation', source_layer: 'process', visibility: 'public',
    });
    expect(byId(rows('nodes/release'), 'Rel:1.0.0')).toMatchObject({branch: 'release/1.0.0', kind: 'minor', status: 'proposed'});
    expect(edges(rows('edges/targets-release'), 'GROK-101')).toEqual([expect.objectContaining({to: 'Rel:1.0.1', kind: 'fix-version', derived_by: 'annotation', evidence: [RECORD]})]);
  });

  it('skips a pick this checkout has no commit for and leaves the git source partial', async () => {
    const {rows, manifest} = await graph;
    expect(rows('nodes/commit')).toEqual([]);
    expect(rows('edges/includes')).toEqual([]);
    expect(manifest.sources.git).toBe('partial');
  });

  it('resolves a pick to the full sha, the release including it, the ticket it closes and the one it only names', async () => {
    const {rows, manifest} = await build({}, commits);
    expect(manifest.sources.git).toBe('ok');
    const picked = rows('nodes/commit');
    expect(picked).toHaveLength(2);
    expect(picked.map((c) => c.repo).sort()).toEqual(['public', 'reddata']);
    expect(picked[0]).toMatchObject({type: 'commit', sha: expect.stringMatching(/^[0-9a-f]{40}$/), subject: 'fixture', date: expect.any(String),
      provenance: 'git', source_layer: 'process', visibility: 'dev'});
    expect(rows('edges/includes').map((e) => [e.from, e.to])).toEqual(picked.map((c) => ['Rel:1.0.1', c.id]).sort());
    // a pick's ticket is the record's word: resolution only where the pick says `claim: resolves`, and a dry run lowers every edge it makes
    expect(rows('edges/resolves')).toEqual([expect.objectContaining({to: 'GROK-101', derived_by: 'annotation', confidence: 0.7, evidence: [RECORD]})]);
    expect(rows('edges/mentions').filter((e) => e.from.startsWith('commit:')).map((e) => [e.to, e.confidence]).sort())
      .toEqual([['GROK-100', 0.7], ['GROK-101', 0.7]]);
    // targets-release is keyed by its kind as well, so a picked ticket that also carries the fix version keeps both assertions
    expect(edges(rows('edges/targets-release'), 'GROK-100').map((e) => [e.kind, e.confidence]).sort()).toEqual([['fix-version', 1], ['picked', 0.7]]);
  }, 60_000);
});
