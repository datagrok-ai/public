/// The home-document layer (build-plan.md WO-2): one node per home with its members, the edges its
/// frontmatter keys spell, `code:` roots and body citations as ownership claims for WO-4, and the
/// stubs the references need (tickets, declarations, cited pages).
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {TypeSystem, EdgeType, Member, NodeType, isSubtype, concreteAuthored} from '../../types';
import {loadHomes, Home, HomeSet, AnnotatedPage, HOME_IGNORE, GLOB_MAGIC, REPO_PREFIX, isHomeFileMember, splitRefAnchor} from '../../homes';
import {extractCitations, proseLines, Citation} from '../../citations';
import {keyLine, Frontmatter} from '../../frontmatter';
import {Emitter} from '../emitter';
import {Row} from '../normalize';
import {BuildContext, Extractor} from '../registry';
import {PREFIXED_ID, SCHEMED_ID, JIRA_KEY, parseId, ticketId, declId, docId, fileId, languageOf, sourceLayerOf, docKind} from '../ids';

export const homesExtractor: Extractor = {
  name: 'homes',
  layer: 'home-document',
  modes: ['full', 'public'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const homes = ctx.homes ?? loadHomes(ctx.system, ctx.repoRoot);
    const layer = new HomeLayer(ctx.system, ctx.repoRoot, homes, emitter);
    for (const page of homes.pages) layer.emitPage(page);
    for (const home of homes.homes) layer.emitHome(home);
    for (const e of homes.errors) emitter.problem('home_issues', `${e.file}${e.line ? `:${e.line}` : ''}: ${e.code}: ${e.message}`);
    if (homes.errors.length) emitter.source('homes', 'partial');
  },
};

/** What an edge key hangs off: a home, or a page that only annotates. */
interface Subject {
  id: string;
  type: NodeType;
  file: string;
  fm: Frontmatter;
}

class HomeLayer {
  private byId = new Map<string, Home>();
  private byAlias = new Map<string, Home>();

  constructor(private system: TypeSystem, private repoRoot: string, homes: HomeSet, private emitter: Emitter) {
    for (const home of homes.homes) {
      this.byId.set(home.id, home);
      for (const alias of home.aliases) this.byAlias.set(alias, home);
    }
  }

  emitHome(home: Home): void {
    const {type, data} = home;
    const row: Row = {type: type.name, id: home.id, name: home.name, home: home.file, source_layer: sourceLayerOf(home.file), provenance: 'annotation'};
    for (const [key, value] of Object.entries(data)) {
      if (value === null || key === 'feature' || key === 'id' || key === 'type' || key === 'aliases') continue;
      const member = type.members[key];
      if (!member || (isHomeFileMember(member) && !home.yaml)) continue;
      row[key] = member.kind === 'ref' ? this.resolveMember(value, member, `${home.file}: ${key}`) : value;
    }
    for (const member of Object.values(type.members))
      if (isHomeFileMember(member) && !home.yaml) row[member.name] = home.file;
    if (home.aliases.length) row.aliases = home.aliases;
    if (type.members.manual_only && row.manual_only === undefined && data.target_layer === 'manual-only') row.manual_only = true;
    if (row.description === undefined && !home.yaml) row.description = firstParagraph(home.body);
    if (!this.emitter.node(row).accepted) return;
    const subject: Subject = {id: home.id, type, file: home.file, fm: home.fm};
    const claimed = this.emitKeys(subject, data);
    if (!home.yaml) this.emitCitations(subject, home.body, claimed);
  }

  emitPage(page: AnnotatedPage): void {
    const type = this.system.nodes.get('doc-page');
    if (!type) return;
    const data = page.fm.data!;
    const id = docId(page.file);
    this.emitter.stub(id, type.name, typeof data.title === 'string' ? data.title : path.posix.basename(page.file), 'annotation', {path: page.file, kind: docKind(page.file)});
    this.emitKeys({id, type, file: page.file, fm: page.fm}, data);
  }

  /** Every edge key of the subject; returns the files its `code:` roots claim. */
  private emitKeys(subject: Subject, data: Record<string, unknown>): Set<string> {
    const claimed = new Set<string>();
    for (const [key, value] of Object.entries(data)) {
      if (value === null || !Array.isArray(value)) continue;
      if (key === 'edges') {
        for (const item of value) {
          if (!item || typeof item !== 'object') continue;
          const {type, to, ...props} = item as Record<string, unknown>;
          const edge = this.system.edges.get(String(type));
          if (edge && typeof to === 'string') this.emitEdge(subject, edge, 'from', to, props);
        }
        continue;
      }
      const edge = this.system.keys.get(key);
      if (!edge || edge.abstract) continue;
      const targetKey = key === 'code' ? 'path' : 'to';
      for (const item of value) {
        const target = typeof item === 'string' ? item : item && typeof item === 'object' ? (item as Record<string, unknown>)[targetKey] : undefined;
        if (typeof target !== 'string' || !target.trim()) continue;
        const props = typeof item === 'object' ? Object.fromEntries(Object.entries(item as Record<string, unknown>).filter(([k]) => k !== targetKey)) : {};
        if (key === 'code') this.emitCode(subject, edge, target.trim(), props, claimed);
        else this.emitEdge(subject, edge, edge.keySide, target, props);
      }
    }
    return claimed;
  }

  private emitEdge(subject: Subject, edge: EdgeType, subjectSide: 'from' | 'to', target: string, props: Record<string, unknown>): void {
    const otherSide = subjectSide === 'from' ? 'to' : 'from';
    const other = this.resolve(target, edge[otherSide], `${subject.file}: ${edge.key ?? edge.name}`);
    if (!other) return;
    this.stubFor(other);
    const [from, to] = subjectSide === 'from' ? [subject.id, other] : [other, subject.id];
    this.emitter.edge({type: edge.name, from, to, derived_by: 'annotation', confidence: 1, evidence: [subject.file], ...props});
  }

  /** A `code:` root: `path#Anchor` is the declaration itself; anything else expands to files and claims them (rung 2). */
  private emitCode(subject: Subject, edge: EdgeType, target: string, props: Record<string, unknown>, claimed: Set<string>): void {
    const hash = target.indexOf('#');
    const file = repoPath(hash < 0 ? target : target.slice(0, hash));
    if (hash >= 0) {
      const id = declId(file, target.slice(hash + 1));
      this.emitter.stub(id, 'declaration', target.slice(hash + 1), 'annotation', {language: languageOf(file), path: file});
      this.emitter.edge({type: edge.name, from: subject.id, to: id, derived_by: 'annotation', confidence: 1, evidence: [subject.file], ...props});
      return;
    }
    const line = keyLine(subject.fm, 'code');
    for (const p of this.expandRoot(file)) {
      if (!this.emitter.node({type: 'source-file', id: fileId(p), name: path.posix.basename(p), path: p, loc: countLines(path.join(this.repoRoot, p)),
        language: languageOf(p), provenance: 'filesystem', source_layer: sourceLayerOf(p)}).accepted) continue;
      this.emitter.claim({file: p, feature: subject.id, rung: 2, source: 'home', props, line});
      claimed.add(p);
    }
  }

  private expandRoot(p: string): string[] {
    const clean = p.replace(/\/+$/, '');
    if (GLOB_MAGIC.test(clean))
      return globSync(clean, {cwd: this.repoRoot, ignore: HOME_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
    const full = path.join(this.repoRoot, clean);
    if (!fs.existsSync(full)) return [];
    if (fs.statSync(full).isDirectory())
      return globSync(`${clean}/**`, {cwd: this.repoRoot, ignore: HOME_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
    return [clean];
  }

  /** Body citations: implementation files become claims (rung 3), documents become mentions of `doc:` stubs. */
  private emitCitations(subject: Subject, body: string, claimed: Set<string>): void {
    const isFeature = subject.type.root === 'feature';
    const mentions = this.system.edges.get('mentions');
    for (const c of extractCitations(subject.file, body, subject.fm.bodyLine)) {
      if (c.resolved === null || c.resolved === subject.file) continue;
      const resolved = this.resolveDocLink(c);
      if (!fs.existsSync(path.join(this.repoRoot, resolved))) continue;
      if (/\.mdx?$/i.test(resolved)) {
        if (!mentions || !this.system.nodes.has('doc-page')) continue;
        const from = mentions.from.some((t) => isSubtype(this.system, subject.type.name, t)) ? subject.id : docId(subject.file);
        if (from !== subject.id) this.emitter.stub(from, 'doc-page', path.posix.basename(subject.file), 'annotation', {path: subject.file, kind: docKind(subject.file)});
        this.emitter.stub(docId(resolved), 'doc-page', path.posix.basename(resolved), 'annotation', {path: resolved, kind: docKind(resolved)});
        this.emitter.edge({type: 'mentions', from, to: docId(resolved), derived_by: 'annotation', confidence: 1, evidence: [subject.file]});
        continue;
      }
      if (!isFeature || claimed.has(resolved) || [...claimed].some((f) => f.startsWith(`${resolved}/`))) continue;
      this.emitter.claim({file: resolved, feature: subject.id, rung: 3, source: 'home', props: {}, line: c.line});
      claimed.add(resolved);
    }
  }

  /** A Docusaurus link may drop the extension: `tile-viewer` means `tile-viewer.md` beside the page. */
  private resolveDocLink(c: Citation): string {
    const p = c.resolved!;
    if (c.kind === 'backtick' || path.posix.extname(p) || fs.existsSync(path.join(this.repoRoot, p))) return p;
    return [`${p}.md`, `${p}.mdx`, `${p}/index.md`].find((a) => fs.existsSync(path.join(this.repoRoot, a))) ?? p;
  }

  private resolveMember(value: unknown, member: Member, where: string): unknown {
    const one = (v: unknown) => typeof v === 'string' ? this.resolve(v, member.refs!, where) : v;
    if (member.list) return (Array.isArray(value) ? value : [value]).map(one).filter((v) => v !== undefined);
    return one(value);
  }

  /**
   * The canonical id of a reference. An alias resolves to the home's id. A prefixed id is exactly its prefix's type. A bare
   * id is the prefix-less type when the expected union admits one (a feature); otherwise it takes the single prefixed
   * candidate, or the one that has a home, and with several candidates and no home it is ambiguous: no edge, a problem.
   */
  private resolve(value: string, expected: string[], where: string): string | undefined {
    const raw = value.trim().replace(/^~/, '');
    if (JIRA_KEY.test(raw) || /^#\d+$/.test(raw)) return ticketId(raw);
    if (PREFIXED_ID.test(raw)) {
      const {id} = splitRefAnchor(raw);
      return this.lookup([id]) ?? id;
    }
    if (SCHEMED_ID.test(raw)) return raw;
    const {id} = splitRefAnchor(raw);
    const concrete = concreteAuthored(this.system, expected);
    if (concrete.some((t) => !t.prefix)) return this.lookup([id]) ?? id;
    const candidates = [...new Set(concrete.map((t) => `${t.prefix}:${id}`))];
    const found = this.lookup(candidates);
    if (found || candidates.length === 1) return found ?? candidates[0];
    this.emitter.problem('ambiguous_refs', `${where}: '${value}' could be ${candidates.join(' or ')}; write the prefix`);
    return undefined;
  }

  private lookup(candidates: string[]): string | undefined {
    for (const c of candidates) {
      const home = this.byId.get(c) ?? this.byAlias.get(c);
      if (home) return home.id;
    }
    return undefined;
  }

  /** The stub an extracted reference target needs when no extractor has produced the node. */
  private stubFor(id: string): void {
    const parsed = parseId(id);
    if (parsed.form === 'tracker')
      this.emitter.stub(id, 'ticket', id, 'annotation', {tracker: parsed.tracker, key: id, kind: 'unknown', state: 'open'});
    else if (parsed.scheme === 'decl' && parsed.path)
      this.emitter.stub(id, 'declaration', parsed.anchor ?? path.posix.basename(parsed.path), 'annotation', {language: languageOf(parsed.path), path: parsed.path});
    else if (parsed.scheme === 'doc' && parsed.path)
      this.emitter.stub(id, 'doc-page', path.posix.basename(parsed.path), 'annotation', {path: parsed.path, kind: docKind(parsed.path)});
  }
}

/** The first paragraph of prose after the frontmatter: not a heading, table, comment or import line. */
export function firstParagraph(body: string): string | undefined {
  const lines = proseLines(body).map((l) => l.text.trim());
  const start = lines.findIndex((l) => l && !/^(#|\||<!--|---|import\s|:::)/.test(l));
  if (start < 0) return undefined;
  let end = start;
  while (end < lines.length && lines[end] && !/^(#|\||<!--|:::)/.test(lines[end])) end++;
  return lines.slice(start, end).join(' ');
}

/** `landing:` and `infra:` prefixed roots name those repos' folders. */
function repoPath(p: string): string {
  const repo = REPO_PREFIX.exec(p.trim());
  return repo ? `${repo[1]}/${repo[2]}` : p.trim();
}

export function countLines(file: string): number {
  const buffer = fs.readFileSync(file);
  if (!buffer.length) return 0;
  let n = 0;
  for (let at = buffer.indexOf(0x0A); at >= 0; at = buffer.indexOf(0x0A, at + 1)) n++;
  return buffer[buffer.length - 1] === 0x0A ? n : n + 1;
}
