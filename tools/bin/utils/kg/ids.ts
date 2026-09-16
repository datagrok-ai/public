/// Ids of extracted nodes (conventions.md §2.2, build-plan.md "Common contracts"): one constructor per
/// scheme, the inverse `parseId`, and the path-derived conventions the extractors share (language,
/// source layer, location visibility, doc kind), and the `source-file` row every extractor emits for a file.
import * as fs from 'fs';
import * as path from 'path';
import type {Row} from './normalize';

/** Extracted id schemes and the node types they name. */
export const SCHEME_TYPES: Record<string, string[]> = {
  pkg: ['package'], lib: ['library'], func: ['function'], decl: ['declaration'], file: ['source-file'],
  ep: ['endpoint'], table: ['db-table'], semtype: ['semantic-type'], doc: ['doc-page'], test: ['test'],
  suite: ['test-suite'], sample: ['sample'], pr: ['pull-request'], gh: ['ticket'], commit: ['commit'],
  report: ['report'], img: ['image'], chg: ['changelog-entry'], conn: ['connection'], env: ['script-environment'],
  container: ['container'], mig: ['migration'], tutorial: ['tutorial'], query: ['query'], script: ['script'],
};
export const PREFIXED_ID = /^([A-Z][A-Za-z]{0,5}):(.+)$/;
export const SCHEMED_ID = /^([a-z][a-z0-9-]*):(.+)$/;
export const JIRA_KEY = /^GROK-\d+$/;
/** A DocPage reference spelled as the page's repo path (`public/help/.../bar-chart.md`), conventions.md §5.3. */
export const PAGE_PATH = /^[^\s:]+\/[^\s:]+\.mdx?$/i;
/** The user help tree; a page under it documents whatever cites it. */
export const HELP_DIR = 'public/help';
/** A `//help-url:` or `HelpUrl` value: the site-relative or absolute help URL, with any extension, slash and fragment. */
const HELP_URL = /^(?:https?:\/\/(?:[\w-]+\.)*datagrok\.ai)?\/help\/([^#?]+?)(?:\.mdx?)?\/?(?:[#?].*)?$/;
/** Where a folder has to be for a home in it to own it, and where a file has to be to count as an orphan. */
export const CODE_ROOTS = ['core/client/', 'core/server/', 'core/shared/', 'public/packages/', 'public/libraries/', 'public/js-api/'];
export const GITHUB_KEY = /^gh:public#\d+$/;
const PATH_SCHEMES = ['file', 'decl', 'doc', 'mig', 'sample'];
const LANGUAGES: Record<string, string> = {
  '.dart': 'dart', '.ts': 'ts', '.tsx': 'ts', '.js': 'js', '.mjs': 'js', '.cjs': 'js', '.jsx': 'js', '.java': 'java',
  '.py': 'python', '.r': 'r', '.sql': 'sql', '.jl': 'julia', '.m': 'octave', '.grok': 'grok',
};

export function posix(p: string): string {
  return p.replace(/\\/g, '/');
}

export function pkgId(folder: string): string {
  return `pkg:${folder}`;
}

/** `libraries/<folder>` and the JS API itself (`js-api`). */
export function libId(folder: string): string {
  return `lib:${folder}`;
}

export function fileId(file: string): string {
  return `file:${posix(file)}`;
}

/** [name] is `Name` or `Owner.member`; a getter/setter pair shares the name and differs in [accessor]; overloads share one id. */
export function declId(file: string, name: string, accessor?: 'get' | 'set'): string {
  return `decl:${posix(file)}#${name}${accessor ? `:${accessor}` : ''}`;
}

/** Scripts, queries and roles of a package; Dart commands use the package `core`. */
export function funcId(pkg: string, name: string): string {
  return `func:${pkg}:${name}`;
}

export function connId(pkg: string, name: string): string {
  return `conn:${pkg}:${name}`;
}

export function envId(pkg: string, name: string): string {
  return `env:${pkg}:${name}`;
}

export function containerId(pkg: string, folder: string): string {
  return `container:${pkg}:${folder}`;
}

export function semtypeId(name: string): string {
  return `semtype:${name}`;
}

/** `test:<framework>:<path>#<category>/<name>`; a test in no category is `#<name>`. */
export function testId(framework: string, file: string, category: string, name: string): string {
  return `test:${framework}:${posix(file)}#${category ? `${category}/` : ''}${name}`;
}

/** `suite:dg:<Pkg>:<category>`, or `suite:<framework>:<path>` for a framework whose suite is a file. */
export function suiteId(framework: 'dg' | 'playwright' | 'dart', pkgOrFile: string, category?: string): string {
  return framework === 'dg' ? `suite:dg:${pkgOrFile}:${category}` : `suite:${framework}:${posix(pkgOrFile)}`;
}

/** [file] relative to `packages/ApiSamples/scripts`; the extension is dropped. */
export function sampleId(file: string): string {
  return `sample:${posix(file).replace(/\.[^./]+$/, '')}`;
}

export function chgId(pkg: string, version: string, n: number): string {
  return `chg:${pkg}:${version}:${n}`;
}

export function docId(file: string, slug?: string): string {
  return `doc:${posix(file)}${slug ? `#${slug}` : ''}`;
}

export function prId(repo: 'reddata' | 'public', n: number): string {
  return `pr:${repo}#${n}`;
}

export function commitId(repo: 'reddata' | 'public', sha: string): string {
  return `commit:${repo}:${sha}`;
}

export function epId(method: string, route: string): string {
  return `ep:${method.toUpperCase()} ${route}`;
}

export function tableId(schema: string, name: string): string {
  return `table:${schema}.${name}`;
}

export function migId(file: string, cls?: string): string {
  return `mig:${posix(file)}${cls ? `#${cls}` : ''}`;
}

export function custId(slug: string): string {
  return `Cust:${slug}`;
}

export function teamId(slug: string): string {
  return `Team:${slug}`;
}

export function relId(version: string): string {
  return `Rel:${version}`;
}

/** Tracker keys as the graph writes them: `GROK-n` as is, `#n` / `public-n` as `gh:public#n`. */
export function ticketId(key: string): string {
  const gh = /^(?:#|public-)(\d+)$/.exec(key);
  return gh ? `gh:public#${gh[1]}` : key;
}

export interface ParsedId {
  form: 'prefix' | 'scheme' | 'tracker' | 'bare';
  prefix?: string;
  scheme?: string;
  tracker?: 'jira' | 'github';
  /** After the prefix or scheme; the whole id when bare. */
  local: string;
  /** For file-addressed schemes (file, decl, doc, mig, sample): the posix path before any `#`. */
  path?: string;
  /** The `#` part: a declaration name, a heading slug, a test name, an issue number. */
  anchor?: string;
  /** `:`-separated parts of a package-scoped id (func, conn, env, container, chg, suite, test, commit). */
  parts?: string[];
}

export function parseId(id: string): ParsedId {
  if (JIRA_KEY.test(id)) return {form: 'tracker', tracker: 'jira', local: id};
  if (GITHUB_KEY.test(id)) return {form: 'tracker', tracker: 'github', scheme: 'gh', local: id.slice(3), anchor: id.slice(id.indexOf('#') + 1)};
  const prefixed = PREFIXED_ID.exec(id);
  if (prefixed) return {form: 'prefix', prefix: prefixed[1], local: prefixed[2]};
  const schemed = SCHEMED_ID.exec(id);
  if (!schemed) return {form: 'bare', local: id};
  const [, scheme, local] = schemed;
  const hash = local.indexOf('#');
  const anchor = hash < 0 ? undefined : local.slice(hash + 1);
  const head = hash < 0 ? local : local.slice(0, hash);
  if (PATH_SCHEMES.includes(scheme)) return {form: 'scheme', scheme, local, path: head, anchor};
  return {form: 'scheme', scheme, local, anchor, parts: head.split(':')};
}

/** A readable name for a node that exists only as a reference target. */
export function stubName(id: string): string {
  const p = parseId(id);
  if (p.form === 'tracker') return id;
  if (p.anchor !== undefined && p.scheme !== 'gh') return p.anchor;
  if (p.path !== undefined) return path.posix.basename(p.path);
  if (p.parts) return p.parts[p.parts.length - 1];
  return titleCase(p.local.split('/').pop()!);
}

/** `ScatterPlot` -> `scatter-plot`, `MLMethods` -> `ml-methods`, `initial runs` -> `initial-runs`. */
export function kebab(segment: string): string {
  return segment.replace(/([a-z0-9])([A-Z])/g, '$1-$2').replace(/([A-Z]+)([A-Z][a-z])/g, '$1-$2').replace(/[\s_]+/g, '-').toLowerCase().replace(/-+/g, '-').replace(/^-|-$/g, '');
}

/** A Docusaurus link may drop the extension: `tile-viewer` is `tile-viewer.md`, `tile-viewer.mdx` or `tile-viewer/index.md`. */
export function docCandidates(p: string): string[] {
  return [`${p}.md`, `${p}.mdx`, `${p}/index.md`];
}

/** The help page a help URL points at, as a repo path: `.../help/datagrok/project` is one of `docCandidates` or
 * `project/project.md` under public/help, whichever exists. */
export function helpPage(repoRoot: string, url: string): string | undefined {
  const m = HELP_URL.exec(url.trim());
  if (!m) return undefined;
  const p = `${HELP_DIR}/${m[1]}`;
  return [...docCandidates(p), `${p}/${path.posix.basename(p)}.md`].find((c) => fs.existsSync(path.join(repoRoot, c)));
}

/** Lines of [text]: one per newline, plus one for an unterminated last line; empty is 0. */
export function countLines(text: string | Buffer): number {
  if (!text.length) return 0;
  const buffer = typeof text === 'string' ? Buffer.from(text, 'utf8') : text;
  let n = 0;
  for (let at = buffer.indexOf(0x0A); at >= 0; at = buffer.indexOf(0x0A, at + 1)) n++;
  return buffer[buffer.length - 1] === 0x0A ? n : n + 1;
}

/** The `source-file` row for [file]: what its path says, then [extra] (`loc`, `generated`, `package`). */
export function sourceFileRow(file: string, extra: Row): Row {
  return {type: 'source-file', id: fileId(file), name: path.posix.basename(file), path: file, language: languageOf(file), provenance: 'filesystem',
    source_layer: sourceLayerOf(file), ...extra};
}

/** `scatter-plot` -> `Scatter plot`. */
export function titleCase(segment: string): string {
  const words = segment.replace(/[-_]+/g, ' ');
  return words.charAt(0).toUpperCase() + words.slice(1);
}

export function languageOf(file: string): string {
  return LANGUAGES[path.posix.extname(posix(file)).toLowerCase()] ?? 'other';
}

/** `source_layer` of a node by the location of its file (build-plan.md Decisions "Layers"). */
export function sourceLayerOf(file: string): string {
  const p = posix(file);
  return p.startsWith('core/') ? 'core' : p.startsWith('infra/') ? 'infra' : 'public';
}

/** Visibility of a source or evidence path: public/ -> public, core/ -> dev, the internal folder -> internal. */
export function locationVisibility(file: string): string {
  const p = posix(file);
  if (p.startsWith('core/docs/knowledge-graph/internal/')) return 'internal';
  return p.startsWith('public/') || p.startsWith('landing/') ? 'public' : 'dev';
}

/** `doc-page.kind` by folder (build-plan.md WO-3c). */
export function docKind(file: string): string {
  const p = posix(file);
  const base = path.posix.basename(p);
  if (p.startsWith('public/help/')) return 'help';
  if (p.startsWith('core/docs/design/')) return 'design';
  if (p.startsWith('core/docs/runbooks/')) return 'runbook';
  if (p.startsWith('core/docs/reviews/')) return 'review';
  if (p.startsWith('core/docs/features/')) return 'record';
  if (p.startsWith('core/docs/plans/')) return 'plan';
  if (/^core\/docs\/[^/]+$/.test(p)) return 'core-doc';
  if (base === 'CLAUDE.md') return 'agent';
  if (/^readme\.mdx?$/i.test(base)) return 'readme';
  return 'other';
}
