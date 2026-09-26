/// `grok kg enrich media` (conventions.md §11.1): the work list over a built generation, the prompt, the answer parsed and
/// validated, the proposal merged into media.yaml, the cache, and what it refuses; the model and ffmpeg are seams here.
import {describe, it, expect} from 'vitest';
import fs from 'fs';
import path from 'path';
import * as yaml from 'js-yaml';
import {enrichMedia, matches, parseAnswer, validate, Prepared, Describer, FrameSampler} from '../utils/kg/enrich/media';
import {loadTypeSystem} from '../utils/kg/types';
import {copyFixture, git, write, buildFixture, runKg, kgRoot, KG_DIR} from './kg-fixture';

const IMG = 'public/help/visualize/viewers/img';
const BINS = `media:${IMG}/bins.png`;

function makeRepo(): string {
  const repo = copyFixture('build');
  // bins.png loses its description, so it is undescribed and shown; histogram.gif stays reviewed
  const record = path.join(repo, IMG, 'media.yaml');
  fs.writeFileSync(record, fs.readFileSync(record, 'utf8').replace(/^bins\.png:[\s\S]*$/m, 'bins.png:\n  reviewed: false\n'));
  git(repo, 'init', '-q');
  git(repo, 'add', '-A');
  git(repo, 'commit', '-q', '-m', 'fixture');
  return repo;
}

const answer = (overrides: Record<string, unknown> = {}) => ({proposal: {caption: 'The Bins slider', description: 'The context panel with the Bins slider at 12.', actions: [],
  kind: 'screenshot', quality: 'answer', quality_notes: '', illustrates: ['visualize/viewers/histogram', 'nope/feature'], ...overrides}});
const stills: FrameSampler = ({local}) => ({frames: [{file: local}], probe: {width: 1, height: 1}, stillOnly: false});

async function run(repo: string, extra: Record<string, unknown> = {}, describe: Describer = () => answer(), frames: FrameSampler = stills) {
  const built = await buildFixture(repo, 'homes,docs,landing,media');
  const system = loadTypeSystem(path.join(repo, KG_DIR));
  return enrichMedia({kgRoot: kgRoot(repo), repoRoot: repo, outRoot: built.out, system, limit: 20, stale: false, dryRun: false, model: 'claude-sonnet-5', describe, frames, ...extra});
}

describe('grok kg enrich media', () => {
  it('lists the shown, undescribed, unreviewed files most-shown first, with the page context and the feature candidates; a dry run writes nothing', async () => {
    const repo = makeRepo();
    const {items, outcomes} = await run(repo, {dryRun: true});
    expect(outcomes).toEqual([{id: BINS, pages: 1, frames: 0, status: 'dry-run'}]);
    expect(items[0]).toMatchObject({path: `${IMG}/bins.png`, format: 'png', blob: git(repo, 'hash-object', `${IMG}/bins.png`), bytes: 70});
    expect(items[0].pages).toEqual([{page: 'public/help/visualize/viewers/histogram.md', title: 'Histogram', anchor: 'bins', alt: 'The Bins slider', paragraph: expect.stringContaining('Watch')}]);
    expect(items[0].candidates.map((c) => c.id)).toEqual(['visualize/viewers', 'visualize/viewers/histogram']);
    expect(items[0].prompt).toContain('Feature ids that may appear in illustrates');
    expect(items[0].prompt).toContain('- visualize/viewers/histogram (Histogram)');
    expect(fs.readFileSync(path.join(repo, IMG, 'media.yaml'), 'utf8')).toContain('bins.png:\n  reviewed: false\n');
  });

  it('merges the validated answer into media.yaml as a proposal with the tool fields, keeps illustrates to the candidates, and reuses the cache', async () => {
    const repo = makeRepo();
    let calls = 0;
    const describe: Describer = () => { calls++; return answer(); };
    const {outcomes} = await run(repo, {}, describe);
    expect(outcomes).toEqual([{id: BINS, pages: 1, frames: 1, status: 'described', caption: 'The Bins slider', detail: undefined, record: `${IMG}/media.yaml`}]);
    const record = yaml.load(fs.readFileSync(path.join(repo, IMG, 'media.yaml'), 'utf8')) as Record<string, Record<string, unknown>>;
    expect(record['bins.png']).toEqual({caption: 'The Bins slider', description: 'The context panel with the Bins slider at 12.', kind: 'screenshot', quality: 'answer',
      illustrates: ['visualize/viewers/histogram'], reviewed: false, described_by: 'claude-sonnet-5', described_blob: git(repo, 'hash-object', `${IMG}/bins.png`), width: 1, height: 1});
    expect(record['histogram.gif'].reviewed).toBe(true);
    // the blob is quoted, so YAML never reads a digit-only one as a number
    expect(fs.readFileSync(path.join(repo, IMG, 'media.yaml'), 'utf8')).toMatch(/described_blob: '?[0-9a-f]{40}'?\n/);
    // the record now describes the file: a rebuild has it, check accepts it, and a second run finds nothing to do
    const rebuilt = await buildFixture(repo, 'homes,docs,landing,media');
    expect(rebuilt.rows('nodes/media').find((m) => m.id === BINS)).toMatchObject({caption: 'The Bins slider', reviewed: false, described_by: 'claude-sonnet-5'});
    expect(rebuilt.rows('edges/illustrates').find((e) => e.from === BINS)).toMatchObject({derived_by: 'llm', confidence: 0.6, proposed: true});
    expect((await runKg({_: ['kg', 'check'], kg: kgRoot(repo), landing: false})).ok).toBe(true);
    const again = await enrichMedia({kgRoot: kgRoot(repo), repoRoot: repo, outRoot: rebuilt.out, system: loadTypeSystem(path.join(repo, KG_DIR)), limit: 20, stale: false, dryRun: false,
      model: 'claude-sonnet-5', describe, frames: stills});
    expect(again.outcomes).toEqual([]);
    expect(calls).toBe(1);
    // --stale after the file changed: the cached answer is not for this blob, the model is asked again and the record refreshed
    fs.appendFileSync(path.join(repo, IMG, 'bins.png'), 'x');
    const changed = await buildFixture(repo, 'homes,docs,landing,media');
    const stale = await enrichMedia({kgRoot: kgRoot(repo), repoRoot: repo, outRoot: changed.out, system: loadTypeSystem(path.join(repo, KG_DIR)), limit: 20, stale: true, dryRun: false,
      model: 'claude-sonnet-5', describe, frames: stills});
    expect(stale.outcomes.map((o) => o.status)).toEqual(['described']);
    expect(calls).toBe(2);
    expect((yaml.load(fs.readFileSync(path.join(repo, IMG, 'media.yaml'), 'utf8')) as any)['bins.png'].described_blob).toBe(git(repo, 'hash-object', `${IMG}/bins.png`));
  });

  it('never overwrites a reviewed record, reports a failed or malformed answer, and refuses a frontier model', async () => {
    const repo = makeRepo();
    const record = path.join(repo, IMG, 'media.yaml');
    fs.writeFileSync(record, fs.readFileSync(record, 'utf8').replace('bins.png:\n  reviewed: false\n', 'bins.png:\n  reviewed: true\n'));
    const kept = await run(repo);
    expect(kept.outcomes).toEqual([]);
    fs.writeFileSync(record, fs.readFileSync(record, 'utf8').replace('bins.png:\n  reviewed: true\n', 'bins.png:\n  reviewed: false\n'));
    const failed = await run(repo, {}, () => ({error: 'claude exited 1'}));
    expect(failed.outcomes).toEqual([{id: BINS, pages: 1, frames: 1, status: 'failed', detail: 'claude exited 1'}]);
    const bad = await run(repo, {}, () => ({proposal: {caption: 'x', quality: 'superb'} as any}));
    expect(bad.outcomes[0]).toMatchObject({status: 'failed', detail: expect.stringMatching(/quality: "superb" is not one of/)});
    await expect(run(repo, {model: 'claude-fable-5-1'})).rejects.toThrow(/frontier model is refused/);
    await expect(run(repo, {model: 'claude-opus-5'})).rejects.toThrow(/refused/);
  });

  it('parses the fenced block of an answer and validates it: only proposal keys, typed members, quality capped without real frames', () => {
    expect(parseAnswer('Here you go:\n```json\n{"caption": "c", "description": "d"}\n```\n').proposal).toEqual({caption: 'c', description: 'd'});
    expect(parseAnswer('{"caption": "c"}').proposal).toEqual({caption: 'c'});
    expect(parseAnswer('no json here').error).toBe('no JSON block in the answer');
    expect(parseAnswer('```json\n{"caption": \n```').error).toMatch(/malformed JSON/);
    const system = loadTypeSystem(path.join(copyFixture('build'), KG_DIR));
    const item = {id: BINS, format: 'gif', candidates: [{id: 'visualize/viewers', name: 'Viewers'}], stillOnly: true} as Prepared;
    const checked = validate(system, {caption: 'c', description: 'd', quality: 'marketing', reviewed: true, blob: 'x', path: 'y', illustrates: ['visualize/viewers', 'other'], actions: [], ui_text: []} as any, item);
    expect(checked.proposal).toEqual({caption: 'c', description: 'd', quality: 'docs', illustrates: ['visualize/viewers']});
    expect(validate(system, {caption: 'c'}, item).error).toBe('no description in the answer');
    expect(validate(system, {description: 'd', kind: 'poster'} as any, item).error).toMatch(/kind: "poster" is not one of/);
  });

  it('--only takes a prefix, a glob, or a folder with a trailing **', () => {
    const file = `${IMG}/bins.png`;
    expect(matches(file, 'public/help/visualize')).toBe(true);
    expect(matches(file, 'public/help/visualize/')).toBe(true);
    expect(matches(file, 'public/help/compute')).toBe(false);
    expect(matches(file, 'public/help/**/*.png')).toBe(true);
    expect(matches(file, 'public/help/*/img/*.png')).toBe(false);
    expect(matches(file, 'public/help/visualize/**')).toBe(true);
    expect(matches(file, 'public/help/visualize/viewers/img/**')).toBe(true);
    expect(matches(file, 'public/help/compute/**')).toBe(false);
  });
});
