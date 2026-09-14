/// Documents (build-plan.md WO-3c): every markdown page the home discovery walks as a `doc-page` with a `doc-anchor`
/// per heading, its `~id` and ticket mentions, the legacy Test Track scenarios (a `feature:` area tag and no `id:`) as
/// `scenario` nodes automated by their sibling specs, and the tutorials of the Tutorials package.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {discoverHomeFiles, firstHeading} from '../../homes';
import {splitFrontmatter, Frontmatter} from '../../frontmatter';
import {proseLines, slugify} from '../../citations';
import {Emitter} from '../emitter';
import {Row} from '../normalize';
import {BuildContext, Extractor} from '../registry';
import {pkgId, docId, docKind, ticketId, sourceLayerOf} from '../ids';
import {HomeIndex, homesOf, emitMentions, ticketStub, kebab} from './markers';
import {firstParagraph} from './homes';
import {parsePlaywrightTests, playwrightTestId} from './ts/tests';

const TEST_TRACK = 'public/packages/UsageAnalysis/files/TestTrack/';
const TUTORIALS = 'public/packages/Tutorials/src/tracks';
const HEADING = /^(#{1,4})\s+(.+?)\s*(?:\{#([^}]+)\})?\s*#*\s*$/;
const NUMBERED_STEP = /^\s*\d+\.\s+\S/;
const PRIORITIES = ['p0', 'p1', 'p2', 'p3'];
const QUOTED = /(['"`])((?:\\.|(?!\1).)*)\1/g;
/** A changelog is the ts-changelog extractor's source; as a page it would add a heading anchor per released version. */
const PAGE_IGNORE = /(^|\/)CHANGELOG\.mdx?$/i;

export const docsExtractor: Extractor = {
  name: 'docs',
  layer: 'docs',
  modes: ['full', 'public'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const layer = new DocLayer(ctx.repoRoot, emitter, new HomeIndex(homesOf(ctx)));
    for (const file of discoverHomeFiles(ctx.repoRoot))
      if (/\.mdx?$/i.test(file) && !PAGE_IGNORE.test(file)) layer.emitPage(file);
    layer.emitTutorials();
    emitter.source('docs', layer.unresolved ? 'partial' : 'ok');
  },
};

class DocLayer {
  unresolved = 0;

  constructor(private repoRoot: string, private emitter: Emitter, private index: HomeIndex) {}

  emitPage(file: string): void {
    const fm = splitFrontmatter(this.read(file));
    const data = fm.data ?? {};
    if (fm.error) this.emitter.problem('invalid_rows', `${file}:${fm.errorLine ?? 1}: ${fm.error}`);
    const id = docId(file);
    const title = typeof data.title === 'string' ? data.title : undefined;
    const name = title ?? firstHeading(fm.body) ?? path.posix.basename(file);
    const keywords = (Array.isArray(data.keywords) ? data.keywords : typeof data.keywords === 'string' ? [data.keywords] : []).filter((k): k is string => typeof k === 'string');
    this.emitter.node({type: 'doc-page', id, name, description: typeof data.description === 'string' ? data.description : firstParagraph(fm.body), path: file, kind: docKind(file),
      title, mdx: /\.mdx$/i.test(file) || !!data.mdx ? true : undefined, unlisted: data.unlisted === true ? true : undefined, keywords: keywords.length ? keywords : undefined,
      provenance: Object.keys(data).length ? 'annotation' : 'filesystem', source_layer: sourceLayerOf(file)});
    const prose = proseLines(fm.body);
    const seen = new Map<string, number>();
    for (const {text} of prose) {
      const h = HEADING.exec(text);
      if (!h) continue;
      const base = h[3] ?? slugify(h[2]);
      if (!base) continue; // a heading of non-Latin words slugifies to nothing, so nothing can link to it
      const n = seen.get(base) ?? 0;
      seen.set(base, n + 1);
      const slug = n ? `${base}-${n}` : base;
      this.emitter.node({type: 'doc-anchor', id: docId(file, slug), name: h[2], path: file, page: id, slug, level: h[1].length, heading: h[2], provenance: 'annotation', source_layer: sourceLayerOf(file)});
    }
    this.unresolved += emitMentions(this.emitter, id, prose.map((l) => l.text).join('\n'), file, this.index).unresolved.length;
    if (file.startsWith(TEST_TRACK) && data.feature !== undefined && data.id === undefined) this.emitScenario(file, fm, name, prose.map((l) => l.text));
  }

  /** A legacy scenario: `TS:<folder segments>/<stem>` in kebab case; `feature:` is the area and never a covers. */
  private emitScenario(file: string, fm: Frontmatter, name: string, prose: string[]): void {
    const data = fm.data!;
    const rel = file.slice(TEST_TRACK.length).replace(/\.md$/i, '');
    const id = `TS:${rel.split('/').map(kebab).join('/')}`;
    const str = (key: string) => typeof data[key] === 'string' ? data[key] as string : undefined;
    const steps = prose.filter((l) => NUMBERED_STEP.test(l)).length;
    const row: Row = {type: 'scenario', id, name, description: firstParagraph(fm.body), path: file, priority: PRIORITIES.includes(str('priority') ?? '') ? str('priority') : undefined,
      target_layer: str('target_layer'), coverage_type: str('coverage_type'), manual_only: data.manual_only === true || data.target_layer === 'manual-only' ? true : undefined,
      manual_only_reason: str('manual_only_reason'), steps: steps || undefined, provenance: 'annotation', source_layer: 'public'};
    this.emitter.node(row);
    const dir = path.posix.dirname(file);
    for (const spec of Array.isArray(data.realized_as) ? data.realized_as : []) {
      if (typeof spec !== 'string') continue;
      const specFile = `${dir}/${spec.trim()}`;
      if (!fs.existsSync(path.join(this.repoRoot, specFile))) {
        this.unresolved++;
        this.emitter.problem('unresolved_ids', `${file}: realized_as ${spec} is not beside the scenario`);
        continue;
      }
      for (const t of parsePlaywrightTests(this.read(specFile)))
        this.emitter.edge({type: 'automates', from: playwrightTestId(specFile, t), to: id, derived_by: 'annotation', confidence: 1, evidence: [file]});
    }
    for (const bug of Array.isArray(data.related_bugs) ? data.related_bugs : []) {
      const key = typeof bug === 'string' ? bug : bug && typeof bug === 'object' ? (bug as Record<string, unknown>).id : undefined;
      if (typeof key !== 'string' || !/^(GROK-\d+|#\d+)$/.test(key.trim())) continue;
      const ticket = ticketId(key.trim());
      ticketStub(this.emitter, ticket);
      this.emitter.edge({type: 'mentions', from: id, to: ticket, derived_by: 'annotation', confidence: 1, evidence: [file]});
    }
  }

  /** `tracks/<dir>/index.ts` names the track; every `class X extends Tutorial` under `tracks/<dir>/tutorials/` is a tutorial. */
  emitTutorials(): void {
    for (const indexFile of this.glob(`${TUTORIALS}/*/index.ts`)) {
      const dir = path.posix.basename(path.posix.dirname(indexFile));
      const track = /new\s+Track\s*\(\s*(['"`])((?:\\.|(?!\1).)*)\1/.exec(this.read(indexFile))?.[2];
      if (!track) continue;
      for (const file of this.glob(`${TUTORIALS}/${dir}/tutorials/*.ts`)) {
        const text = this.read(file);
        const name = /class\s+\w+\s+extends\s+Tutorial\b[\s\S]*?get\s+name\s*\(\)\s*\{\s*return\s+(['"`])((?:\\.|(?!\1).)*)\1/.exec(text)?.[2];
        if (!name) continue;
        const description = /get\s+description\s*\(\)\s*\{\s*return\s+([\s\S]*?);\s*\}/.exec(text)?.[1];
        this.emitter.node({type: 'tutorial', id: `tutorial:${dir}/${path.posix.basename(file, '.ts')}`, name, path: file, track, package: pkgId('Tutorials'),
          description: description ? [...description.matchAll(QUOTED)].map((m) => m[2]).join('') || undefined : undefined, provenance: 'ast', source_layer: 'public'});
      }
    }
  }

  private glob(pattern: string): string[] {
    return globSync(pattern, {cwd: this.repoRoot, ignore: ['**/node_modules/**'], nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
  }

  private read(file: string): string {
    return fs.readFileSync(path.join(this.repoRoot, file), 'utf8');
  }
}
