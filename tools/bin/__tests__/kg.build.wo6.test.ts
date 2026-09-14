/// The Dart batch consumer (build-plan.md WO-6): the node, edge and claim lines of
/// `.kg/batches/kg-dart.jsonl` and the four states the manifest reports for them.
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {spawnSync} from 'child_process';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const BATCH = path.join('.kg', 'batches', 'kg-dart.jsonl');
const NOW = new Date().toISOString();

/** The fixture as a git repository, so the batch revision can match HEAD or not. */
function makeRepo(): string {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-dart-'));
  fs.cpSync(fixture, repo, {recursive: true});
  spawnSync('git', ['init', '-q'], {cwd: repo});
  spawnSync('git', ['-c', 'user.email=kg@test', '-c', 'user.name=kg', 'commit', '--allow-empty', '-q', '-m', 'fixture'], {cwd: repo});
  return repo;
}

function head(repo: string): string {
  return spawnSync('git', ['rev-parse', '--verify', 'HEAD'], {cwd: repo, encoding: 'utf8'}).stdout.trim();
}

/** Rewrites the batch envelope; `null` deletes the file. */
function batch(repo: string, envelope: Record<string, unknown> | null, extra: string[] = []): void {
  const file = path.join(repo, BATCH);
  if (!envelope) {
    fs.rmSync(file);
    return;
  }
  const lines = fs.readFileSync(file, 'utf8').split('\n').filter(Boolean);
  lines[0] = JSON.stringify({...JSON.parse(lines[0]), ...envelope});
  fs.writeFileSync(file, [...lines, ...extra].join('\n') + '\n');
}

async function build(repo: string): Promise<{manifest: any, rows: (file: string) => any[]}> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'dart', db: false, output: 'json'});
    expect(error.mock.calls).toEqual([]);
    const rows = (file: string) => {
      const p = path.join(repo, '.kg', file.startsWith('reports/') ? file : `data/${file}.jsonl`);
      return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
    };
    return {manifest: JSON.parse(String(log.mock.calls[0][0])), rows};
  }
  finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

describe('the Dart batch (build-plan.md WO-6)', () => {
  it('loads nodes, edges and claims from a batch of the current revision', async () => {
    const repo = makeRepo();
    batch(repo, {built_at: NOW, revision: head(repo)});
    const {manifest, rows} = await build(repo);
    expect(manifest.sources.dart).toBe('ok');
    expect(manifest.problems.invalid_rows).toBe(0);
    expect(rows('nodes/source-file')).toMatchObject([{
      id: 'file:core/server/datlas/lib/src/services/bio_service.dart', type: 'source-file', name: 'bio_service.dart',
      language: 'dart', loc: 214, provenance: 'ast', source_layer: 'core', batch: 'dart:deadbeef',
    }]);
    expect(rows('nodes/declaration').map((r) => r.id)).toEqual([
      'decl:core/server/datlas/lib/src/routers/bio.dart#BioRouter',
      'decl:core/server/datlas/lib/src/services/bio_service.dart#BioService',
      'decl:core/server/datlas/lib/src/services/bio_service.dart#BioService.getSequence',
    ]);
    expect(rows('nodes/endpoint')).toMatchObject([{id: 'ep:GET /bio/sequences/{id}', method: 'GET', route: '/bio/sequences/{id}', path_params: ['id']}]);
    expect(rows('nodes/db-table')).toMatchObject([{id: 'table:public.sequences', schema: 'public'}]);
    expect(rows('edges/declares')).toHaveLength(2);
    // a reference property of a batch node becomes its own edge, like any other
    expect(rows('edges/router')).toMatchObject([{type: 'ref', name: 'router', from: 'ep:GET /bio/sequences/{id}', to: 'decl:core/server/datlas/lib/src/routers/bio.dart#BioRouter'}]);
    expect(rows('reports/claims.jsonl')).toMatchObject([
      {file: 'core/server/datlas/lib/src/routers/bio.dart', feature: 'domains/bio', rung: 1, source: 'marker', mode: 'participates', line: 9},
      {file: 'core/server/datlas/lib/src/services/bio_service.dart', feature: 'domains/bio', rung: 1, source: 'marker', props: {role: 'definition'}, line: 1},
    ]);
  });

  it('reports a batch from another revision or older than a week as stale, and still loads it', async () => {
    const other = makeRepo();
    batch(other, {built_at: NOW, revision: 'deadbeefdeadbeefdeadbeefdeadbeefdeadbeef'});
    const stale = await build(other);
    expect(stale.manifest.sources.dart).toBe('stale');
    expect(stale.rows('nodes/source-file')).toHaveLength(1);

    const old = makeRepo();
    batch(old, {built_at: '2026-01-02T03:04:05Z', revision: head(old)});
    const aged = await build(old);
    expect(aged.manifest.sources.dart).toBe('stale');
    expect(aged.rows('nodes/endpoint')).toHaveLength(1);
  });

  it('reports an absent batch as missing, without a problem', async () => {
    const repo = makeRepo();
    batch(repo, null);
    const {manifest, rows} = await build(repo);
    expect(manifest.sources.dart).toBe('missing');
    expect(manifest.problems.invalid_rows).toBe(0);
    expect(rows('nodes/source-file')).toEqual([]);
  });

  it('reports unreadable lines as partial and keeps the rest', async () => {
    const repo = makeRepo();
    batch(repo, {built_at: NOW, revision: head(repo)}, [
      '{"record":"node","type":',
      '{"record":"verse","type":"source-file","id":"file:core/x.dart"}',
      '{"record":"claim","feature":"domains/bio"}',
    ]);
    const {manifest, rows} = await build(repo);
    expect(manifest.sources.dart).toBe('partial');
    expect(manifest.problems.invalid_rows).toBe(3);
    expect(rows('nodes/source-file')).toHaveLength(1);
    expect(rows('reports/claims.jsonl')).toHaveLength(2);
  });
});
