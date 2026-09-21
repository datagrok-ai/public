/// The media extractor (notation.md §3.7) over the build fixture as a git repository: media nodes from the files
/// under the roots and the records, hosted videos, the embeds of a page with their occurrence properties, illustrates
/// with the provenance the record earns, thumbnails, the broken embed, the public projection and the batch identity.
import {describe, it, expect} from 'vitest';
import fs from 'fs';
import path from 'path';
import {extractEmbeds, resolveTarget} from '../utils/kg/embeds';
import {headings} from '../utils/kg/citations';
import {copyFixture, git, write, buildFixture, Built} from './kg-fixture';

const PAGE = 'public/help/visualize/viewers/histogram.md';
const IMG = 'public/help/visualize/viewers/img';
const GIF = `media:${IMG}/histogram.gif`;
const VIDEO = 'video:youtube:abc123def45';

function makeRepo(): string {
  const repo = copyFixture('build');
  git(repo, 'init', '-q');
  git(repo, 'add', '-A');
  git(repo, 'commit', '-q', '-m', 'fixture');
  return repo;
}

const ONLY = 'homes,docs,media';
const graph: Promise<Built> = buildFixture(makeRepo(), ONLY);

describe('media extractor (notation.md §3.7)', () => {
  it('makes a node per file under the roots with its git blob, bytes, delivery url, thumbnail and the record fields, and per hosted video', async () => {
    const {rows, repo} = await graph;
    const media = rows('nodes/media');
    expect(media.map((m) => m.id)).toEqual([`media:${IMG}/bins.png`, `media:${IMG}/histogram-thumb.png`, GIF, `media:${IMG}/unused.png`, VIDEO]);
    const gif = media.find((m) => m.id === GIF);
    expect(gif).toEqual({id: GIF, type: 'media', name: 'histogram.gif', batch: expect.any(String), blob: git(repo, 'hash-object', `${IMG}/histogram.gif`), bytes: 43,
      format: 'gif', kind: 'animation', path: `${IMG}/histogram.gif`, url: `https://datagrok.ai/help/visualize/viewers/img/histogram.gif`,
      thumbnail: `media:${IMG}/histogram-thumb.png`, caption: 'Changing the bin count', description: expect.stringContaining('drags the Bins slider'),
      actions: ['opens the histogram', 'drags the Bins slider left', 'drags it back'], ui_text: ['Bins', 'Context Panel'], quality: 'answer', reviewed: true,
      described_by: 'person', described_blob: gif.blob, width: 800, height: 500, seconds: 6, provenance: 'annotation', source_layer: 'public', status: 'active', visibility: 'public'});
    const unused = media.find((m) => m.id === `media:${IMG}/unused.png`);
    expect(unused).toMatchObject({provenance: 'filesystem', format: 'png', bytes: 70});
    expect(Object.keys(unused)).not.toContain('reviewed');
    expect(Object.keys(unused)).not.toContain('thumbnail');
    expect(media.find((m) => m.id === VIDEO)).toMatchObject({name: 'Histograms in five minutes', format: 'youtube', provider: 'youtube', external_id: 'abc123def45',
      url: 'https://www.youtube.com/watch?v=abc123def45', thumbnail: `media:${IMG}/histogram-thumb.png`, quality: 'marketing', reviewed: true, seconds: 300, provenance: 'annotation'});
    expect(rows('edges/thumbnail').map((e) => [e.from, e.to])).toEqual([[GIF, `media:${IMG}/histogram-thumb.png`], [VIDEO, `media:${IMG}/histogram-thumb.png`]]);
  });

  it('draws one embeds edge per occurrence, numbered, with the form, the nearest heading, alt, title and the start of a timestamped link', async () => {
    const {rows} = await graph;
    const embeds = rows('edges/embeds').filter((e) => e.from === `doc:${PAGE}`).sort((a, b) => a.position - b.position);
    expect(embeds.map(({batch, evidence, derived_by, confidence, type, from, ...rest}) => rest)).toEqual([
      {to: GIF, position: 1, line: 9, form: 'image', alt: 'Changing bins', title: 'Bins'},
      {to: `media:${IMG}/bins.png`, position: 2, line: 13, form: 'tag', anchor: 'bins', alt: 'The Bins slider'},
      {to: VIDEO, position: 3, line: 15, form: 'link', anchor: 'bins', title: 'the lesson', start_seconds: 90},
      {to: VIDEO, position: 4, line: 16, form: 'link', anchor: 'bins', alt: 'Histogram lesson'},
      {to: VIDEO, position: 5, line: 18, form: 'iframe', anchor: 'bins'},
    ]);
    expect(embeds[0]).toMatchObject({derived_by: 'ast', confidence: 1, evidence: [PAGE]});
  });

  it('asserts illustrates from the records only: annotation 1.0 when reviewed, an llm proposal at 0.6 otherwise', async () => {
    const {rows} = await graph;
    expect(rows('edges/illustrates').map((e) => [e.from, e.to, e.derived_by, e.confidence, e.evidence])).toEqual([
      [`media:${IMG}/bins.png`, 'visualize/viewers/histogram', 'llm', 0.6, [`${IMG}/media.yaml`]],
      [GIF, 'visualize/viewers/histogram', 'annotation', 1, [`${IMG}/media.yaml`]],
      [VIDEO, 'visualize/viewers/histogram', 'annotation', 1, ['public/help/videos.yaml']],
    ]);
  });

  it('counts an embed of a file that does not exist as broken, makes no node for it and marks the media source partial', async () => {
    const {manifest, problems} = await graph;
    expect(manifest.sources.media).toBe('partial');
    expect(manifest.problems.broken_embeds).toBe(1);
    expect(problems.broken_embeds).toEqual([`${PAGE}:20: ![gone](img/gone.png)`]);
    expect(manifest.problems.untracked_media).toBeUndefined();
  });

  it('reports a local addition no page shows instead of indexing it, and a new record file moves the batch', async () => {
    const {manifest: before, repo} = await graph;
    write(repo, 'public/help/domains/bio/img/new.png', 'not really a png');
    write(repo, 'public/help/domains/bio/img/media.yaml', 'new.png:\n  caption: a local addition\n');
    const {manifest, rows, problems} = await buildFixture(repo, ONLY);
    expect(manifest.batch).not.toBe(before.batch);
    // described by a record, so it is indexed with the hash git would give it; the report still names it
    expect(rows('nodes/media').find((m) => m.path === 'public/help/domains/bio/img/new.png')).toMatchObject({caption: 'a local addition', blob: git(repo, 'hash-object', 'public/help/domains/bio/img/new.png')});
    expect(problems.untracked_media).toEqual(['public/help/domains/bio/img/new.png']);
    fs.rmSync(path.join(repo, 'public/help/domains/bio/img/media.yaml'));
    const again = await buildFixture(repo, ONLY);
    expect(again.rows('nodes/media').some((m) => m.path === 'public/help/domains/bio/img/new.png')).toBe(false);
    expect(again.problems.untracked_media).toEqual(['public/help/domains/bio/img/new.png']);
  });

  it('keeps media in the public projection with their thumbnails and illustrates, and drops the record evidence paths that are not public', async () => {
    const {rows, manifest} = await buildFixture(makeRepo(), ONLY, {public: true});
    expect(manifest.counts.nodes.media).toBe(5);
    const gif = rows('nodes/media').find((m) => m.id === GIF);
    expect(gif).toMatchObject({thumbnail: `media:${IMG}/histogram-thumb.png`, description: expect.any(String)});
    expect(Object.keys(gif)).not.toContain('home');
    expect(rows('edges/illustrates').map((e) => e.from)).toEqual([`media:${IMG}/bins.png`, GIF, VIDEO]);
    expect(rows('edges/embeds')).toHaveLength(5);
  });
});

describe('report media', () => {
  it('lists the backlog: undescribed, stale, awaiting review, unreferenced, broken, no alt, duplicates, largest, unfit', async () => {
    const {report} = await graph;
    const media = report('media');
    expect(media.summary).toBe('5 media (0 MB): 0 shown but undescribed, 1 described before the file changed, 1 awaiting review, 1 shown nowhere, 1 broken embeds, ' +
      '0 images without alt text, 1 files committed more than once, 0 rated unfit and still shown.');
    const rows = (title: string) => media.sections.find((s: any) => s.title === title).rows;
    expect(rows('stale')).toEqual([{id: `media:${IMG}/bins.png`, format: 'png', bytes: 70, pages: 1, described_blob: 'deadbeefdead', blob: expect.any(String)}]);
    expect(rows('awaiting review')).toEqual([{id: `media:${IMG}/bins.png`, format: 'png', bytes: 70, pages: 1, described_by: 'claude-sonnet-5', quality: 'docs'}]);
    // the thumbnail some gif names is not unreferenced; the png nothing shows is
    expect(rows('unreferenced')).toEqual([{id: `media:${IMG}/unused.png`, format: 'png', bytes: 70, path: `${IMG}/unused.png`}]);
    expect(rows('broken embeds')).toEqual([{embed: `${PAGE}:20: ![gone](img/gone.png)`}]);
    expect(rows('duplicates')).toEqual([{blob: expect.any(String), bytes: 70, copies: 3, paths: `${IMG}/bins.png, ${IMG}/histogram-thumb.png, ${IMG}/unused.png`}]);
    expect(rows('largest')[0]).toEqual({id: `media:${IMG}/bins.png`, format: 'png', bytes: 70, pages: 1});
    expect(rows('untracked')).toEqual([]);
  });
});

describe('embeds of a page (embeds.ts)', () => {
  it('reads every embed form with its target, alt, caption and start, skips fenced code and external images, and resolves site-absolute paths', () => {
    const body = [
      "import Poster from './img/poster.png';",
      '',
      '# Top',
      '',
      '![Inline](img/a.gif "A title") and ![Ref][pic] and ![Short]',
      '',
      '## Deeper {#deep}',
      '',
      '<img src={Poster} alt="Poster" width="200"/> <img src={require("./img/b.png").default} />',
      '<video poster="img/frame.png" loop><source src="img/clip.mp4" type="video/mp4"></video>',
      '<Image src="/docusaurus_img/slides/c.png" />',
      '<Iframe url="https://www.youtube.com/embed/ABCDEFGHIJK?rel=0" />',
      '<iframe src="https://www.youtube.com/embed/ABCDEFGHIJK?start=30"',
      '  width="512"></iframe>',
      'See https://youtu.be/ABCDEFGHIJK?t=5 and [watch](https://www.youtube.com/watch?v=ABCDEFGHIJK&t=45s).',
      '![Site](/help/uploads/d.png) ![Landing](/img/slides/e.gif) ![Ext](https://example.com/x.png) ![Esc](img/my%20file.png)',
      '{text: "Forty connectors", image: "/docusaurus_img/slides/access.png"}',
      '',
      '```',
      '![Fenced](img/fenced.png)',
      '```',
      '',
      '[pic]: img/ref.png "Ref title"',
    ].join('\n');
    const embeds = extractEmbeds('public/help/access/page.md', body, 4, headings(body));
    expect(embeds.map((e) => [e.position, e.line, e.form, e.target.kind === 'file' ? e.target.path : `${e.target.provider}:${e.target.id}`, e.anchor ?? '', e.alt ?? '', e.title ?? e.caption ?? '', e.start_seconds ?? '', e.poster ?? ''])).toEqual([
      [1, 8, 'image', 'public/help/access/img/a.gif', 'top', 'Inline', 'A title', '', ''],
      [2, 8, 'image', 'public/help/access/img/ref.png', 'top', 'Ref', 'Ref title', '', ''],
      [3, 12, 'tag', 'public/help/access/img/poster.png', 'deep', 'Poster', '', '', ''],
      [4, 12, 'tag', 'public/help/access/img/b.png', 'deep', '', '', '', ''],
      [5, 13, 'tag', 'public/help/access/img/clip.mp4', 'deep', '', '', '', 'public/help/access/img/frame.png'],
      [6, 14, 'tag', 'public/docusaurus/static/docusaurus_img/slides/c.png', 'deep', '', '', '', ''],
      [7, 15, 'iframe', 'youtube:ABCDEFGHIJK', 'deep', '', '', '', ''],
      [8, 16, 'iframe', 'youtube:ABCDEFGHIJK', 'deep', '', '', 30, ''],
      [9, 18, 'link', 'youtube:ABCDEFGHIJK', 'deep', '', '', 5, ''],
      [10, 18, 'link', 'youtube:ABCDEFGHIJK', 'deep', '', 'watch', 45, ''],
      [11, 19, 'image', 'public/help/uploads/d.png', 'deep', 'Site', '', '', ''],
      [12, 19, 'image', 'landing:web/img/slides/e.gif', 'deep', 'Landing', '', '', ''],
      [13, 19, 'image', 'public/help/access/img/my file.png', 'deep', 'Esc', '', '', ''],
      [14, 20, 'slide', 'public/docusaurus/static/docusaurus_img/slides/access.png', 'deep', '', 'Forty connectors', '', ''],
    ]);
    expect(resolveTarget('../../core/docs/x.png', 'public/help')).toEqual({kind: 'file', path: 'core/docs/x.png'});
    expect(resolveTarget('../../../etc/x.png', 'public/help')).toBeNull();
    expect(resolveTarget('notes.md', 'public/help')).toBeNull();
    expect(resolveTarget('https://www.youtube.com/watch?v=abc123def45&t=1048s', 'public/help')).toEqual({kind: 'hosted', provider: 'youtube', id: 'abc123def45'});
  });
});
