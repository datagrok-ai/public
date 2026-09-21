/// Media records (conventions.md §5.7, notation.md §3.7): the media.yaml and videos.yaml entries `grok kg check` reads
/// with the homes, the errors and the stale-description warning they can raise, and the blob ids the build keys them by.
import {describe, it, expect} from 'vitest';
import fs from 'fs';
import path from 'path';
import {loadTypeSystem, TypeSystem} from '../utils/kg/types';
import {loadHomes, makeReport, HomeSet} from '../utils/kg/homes';
import {blobIds, hashBlob, RECORD_KEYS} from '../utils/kg/media';
import {copyFixture, git, write, runKg, kgRoot, KG_DIR} from './kg-fixture';

const DIR = 'public/help/visualize/viewers/img';
const RECORD = `${DIR}/media.yaml`;
const PNG = `${DIR}/scatter-plot.png`;

function load(repo: string): {system: TypeSystem, homes: HomeSet} {
  const system = loadTypeSystem(path.join(repo, KG_DIR));
  return {system, homes: loadHomes(system, repo)};
}

function issues(repo: string): string[] {
  const {homes} = load(repo);
  return [...homes.errors, ...homes.warnings].filter((e) => /(media|videos).yaml$/.test(e.file)).map((e) => `${e.code} ${e.file}:${e.line}: ${e.message}`);
}

describe('media records (conventions.md §5.7)', () => {
  it('reads the good record: the members it sets, the illustrates targets, one record counted by check', async () => {
    const repo = copyFixture('good');
    const {system, homes} = load(repo);
    expect(homes.errors).toEqual([]);
    expect(homes.warnings.filter((w) => w.file.endsWith(".yaml"))).toEqual([]);
    expect(homes.media.get(`media:${PNG}`)).toEqual({id: `media:${PNG}`, path: PNG, file: RECORD, line: 1,
      illustrates: ['visualize/viewers/scatter-plot', 'C:dataframe'],
      data: {kind: 'screenshot', caption: 'A scatter plot of two columns', description: 'Points of the demo table plotted by age and weight, coloured by sex.',
        quality: 'answer', reviewed: true, described_by: 'person'}});
    expect(makeReport(system, homes).media).toBe(1);
    const run = await runKg({_: ['kg', 'check'], kg: kgRoot(repo)});
    expect(run.ok).toBe(true);
    expect(run.out.join('\n')).toMatch(/1 media files checked, 1 media record;/);
  });

  it('reports an entry naming no file, a key the record may not set, a bad illustrates target and a wrong value, each at its line', () => {
    const repo = copyFixture('good');
    write(repo, RECORD, [
      'gone.png:', '  caption: x',
      'scatter-plot.png:', '  path: elsewhere.png', '  caption: ok', '  illustrates: [visualize/nothing]', '  width: wide',
      '../scatter-plot.md:', '  caption: y',
      '',
    ].join('\n'));
    expect(issues(repo)).toEqual([
      `unknown-media ${RECORD}:1: 'gone.png' names no media file beside ${RECORD}`,
      `unknown-key ${RECORD}:3: 'scatter-plot.png': unknown record key 'path' (allowed: ${RECORD_KEYS.join(', ')})`,
      expect.stringMatching(new RegExp(`^unresolved-ref ${RECORD}:3: 'scatter-plot.png': illustrates\\[0\\]: 'visualize/nothing' does not resolve`)),
      expect.stringMatching(new RegExp(`^bad-value ${RECORD}:3: 'scatter-plot.png': width: expected a finite number, got \\"wide\\"$`)),
      `unknown-media ${RECORD}:8: '../scatter-plot.md' names no media file beside ${RECORD}`,
    ]);
    expect(load(repo).homes.media.size).toBe(0);
    expect(RECORD_KEYS).not.toContain('path');
    expect(RECORD_KEYS).not.toContain('blob');
    expect(RECORD_KEYS).not.toContain('url');
    expect(RECORD_KEYS).not.toContain('visibility');
    expect(RECORD_KEYS).not.toContain('provenance');
  });

  it('reads a hosted video from videos.yaml by its youtube key, with a title for its name, and refuses any other key shape', () => {
    const repo = copyFixture('good');
    write(repo, 'public/help/videos.yaml', 'youtube:abcdefghijk:\n  title: Scatter plots\n  caption: A lesson\n  illustrates: [visualize/viewers/scatter-plot]\n  reviewed: true\nvimeo:123:\n  caption: no\n');
    const {homes} = load(repo);
    expect(homes.errors.map((e) => `${e.code} ${e.line}: ${e.message}`)).toEqual(["unknown-media 6: 'vimeo:123' is not a hosted video key (youtube:<id>)"]);
    expect(homes.media.get('video:youtube:abcdefghijk')).toEqual({id: 'video:youtube:abcdefghijk', provider: 'youtube', externalId: 'abcdefghijk',
      file: 'public/help/videos.yaml', line: 1, illustrates: ['visualize/viewers/scatter-plot'], data: {name: 'Scatter plots', caption: 'A lesson', reviewed: true}});
  });

  it('refuses two records for one file and a record that is not a map', () => {
    const repo = copyFixture('good');
    write(repo, 'public/help/visualize/viewers/media.yaml', 'img/scatter-plot.png:\n  caption: again\nimg/scatter-plot.png#2: 5\n');
    // record files are read in path order, so the folder's own media.yaml holds the record and the parent's entry is the duplicate
    expect(issues(repo)).toEqual([
      `duplicate-id public/help/visualize/viewers/media.yaml:1: 'img/scatter-plot.png': media:${PNG} already has a record`,
      expect.stringMatching(/^unknown-media public\/help\/visualize\/viewers\/media.yaml:3: 'img\/scatter-plot.png#2' names no media file/),
    ]);
  });

  it('warns when the file changed since it was described, keyed by the git blob of the file', () => {
    const repo = copyFixture('good');
    git(repo, 'init', '-q');
    git(repo, 'add', '.');
    git(repo, 'commit', '-q', '-m', 'fixture');
    const blob = git(repo, 'hash-object', PNG);
    expect(blobIds(repo).get(PNG)).toBe(blob);
    expect(hashBlob(fs.readFileSync(path.join(repo, PNG)))).toBe(blob);
    const record = fs.readFileSync(path.join(repo, RECORD), 'utf8');
    write(repo, RECORD, `${record}  described_blob: ${blob}\n`);
    expect(issues(repo)).toEqual([]);
    // the working tree changed the file: the blob is what git would hash now, and the description is stale
    fs.appendFileSync(path.join(repo, PNG), 'x');
    const changed = git(repo, 'hash-object', PNG);
    expect(blobIds(repo).get(PNG)).toBe(changed);
    expect(issues(repo)).toEqual([`stale-description ${RECORD}:1: 'scatter-plot.png': described at blob ${blob.slice(0, 12)}, the file is now ${changed.slice(0, 12)}; run grok kg enrich media --stale`]);
    // a tree that is no repository: nothing from git, the raw hash stands in
    expect(blobIds(copyFixture('good')).size).toBe(0);
  });
});
