"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.SCHEME_TYPES = exports.SCHEMED_ID = exports.PREFIXED_ID = exports.JIRA_KEY = exports.GITHUB_KEY = void 0;
exports.chgId = chgId;
exports.commitId = commitId;
exports.connId = connId;
exports.containerId = containerId;
exports.custId = custId;
exports.declId = declId;
exports.docId = docId;
exports.docKind = docKind;
exports.envId = envId;
exports.epId = epId;
exports.fileId = fileId;
exports.funcId = funcId;
exports.languageOf = languageOf;
exports.libId = libId;
exports.locationVisibility = locationVisibility;
exports.migId = migId;
exports.parseId = parseId;
exports.pkgId = pkgId;
exports.posix = posix;
exports.prId = prId;
exports.relId = relId;
exports.sampleId = sampleId;
exports.semtypeId = semtypeId;
exports.sourceLayerOf = sourceLayerOf;
exports.stubName = stubName;
exports.suiteId = suiteId;
exports.tableId = tableId;
exports.teamId = teamId;
exports.testId = testId;
exports.ticketId = ticketId;
exports.titleCase = titleCase;
var path = _interopRequireWildcard(require("path"));
function _interopRequireWildcard(e, t) { if ("function" == typeof WeakMap) var r = new WeakMap(), n = new WeakMap(); return (_interopRequireWildcard = function (e, t) { if (!t && e && e.__esModule) return e; var o, i, f = { __proto__: null, default: e }; if (null === e || "object" != typeof e && "function" != typeof e) return f; if (o = t ? n : r) { if (o.has(e)) return o.get(e); o.set(e, f); } for (const t in e) "default" !== t && {}.hasOwnProperty.call(e, t) && ((i = (o = Object.defineProperty) && Object.getOwnPropertyDescriptor(e, t)) && (i.get || i.set) ? o(f, t, i) : f[t] = e[t]); return f; })(e, t); }
/// Ids of extracted nodes (conventions.md §2.2, build-plan.md "Common contracts"): one constructor per
/// scheme, the inverse `parseId`, and the path-derived conventions the extractors share (language,
/// source layer, location visibility, doc kind).

/** Extracted id schemes and the node types they name. */
const SCHEME_TYPES = exports.SCHEME_TYPES = {
  pkg: ['package'],
  lib: ['library'],
  func: ['function'],
  decl: ['declaration'],
  file: ['source-file'],
  ep: ['endpoint'],
  table: ['db-table'],
  semtype: ['semantic-type'],
  doc: ['doc-page'],
  test: ['test'],
  suite: ['test-suite'],
  sample: ['sample'],
  pr: ['pull-request'],
  gh: ['ticket'],
  commit: ['commit'],
  report: ['report'],
  img: ['image'],
  chg: ['changelog-entry'],
  conn: ['connection'],
  env: ['script-environment'],
  container: ['container'],
  mig: ['migration']
};
const PREFIXED_ID = exports.PREFIXED_ID = /^([A-Z][A-Za-z]{0,5}):(.+)$/;
const SCHEMED_ID = exports.SCHEMED_ID = /^([a-z][a-z0-9-]*):(.+)$/;
const JIRA_KEY = exports.JIRA_KEY = /^GROK-\d+$/;
const GITHUB_KEY = exports.GITHUB_KEY = /^gh:public#\d+$/;
const PATH_SCHEMES = ['file', 'decl', 'doc', 'mig', 'sample'];
const LANGUAGES = {
  '.dart': 'dart',
  '.ts': 'ts',
  '.tsx': 'ts',
  '.js': 'js',
  '.mjs': 'js',
  '.cjs': 'js',
  '.jsx': 'js',
  '.java': 'java',
  '.py': 'python',
  '.r': 'r',
  '.sql': 'sql',
  '.jl': 'julia',
  '.m': 'octave',
  '.grok': 'grok'
};
function posix(p) {
  return p.replace(/\\/g, '/');
}
function pkgId(folder) {
  return `pkg:${folder}`;
}

/** `libraries/<folder>` and the JS API itself (`js-api`). */
function libId(folder) {
  return `lib:${folder}`;
}
function fileId(file) {
  return `file:${posix(file)}`;
}

/** [name] is `Name` or `Owner.member`; a getter/setter pair shares the name and differs in [accessor]; overloads share one id. */
function declId(file, name, accessor) {
  return `decl:${posix(file)}#${name}${accessor ? `:${accessor}` : ''}`;
}

/** Scripts, queries and roles of a package; Dart commands use the package `core`. */
function funcId(pkg, name) {
  return `func:${pkg}:${name}`;
}
function connId(pkg, name) {
  return `conn:${pkg}:${name}`;
}
function envId(pkg, name) {
  return `env:${pkg}:${name}`;
}
function containerId(pkg, folder) {
  return `container:${pkg}:${folder}`;
}
function semtypeId(name) {
  return `semtype:${name}`;
}
function testId(framework, file, category, name) {
  return `test:${framework}:${posix(file)}#${category}/${name}`;
}

/** `suite:dg:<Pkg>:<category>` or `suite:playwright:<path>`. */
function suiteId(framework, pkgOrFile, category) {
  return framework === 'dg' ? `suite:dg:${pkgOrFile}:${category}` : `suite:playwright:${posix(pkgOrFile)}`;
}

/** [file] relative to `packages/ApiSamples/scripts`; the extension is dropped. */
function sampleId(file) {
  return `sample:${posix(file).replace(/\.[^./]+$/, '')}`;
}
function chgId(pkg, version, n) {
  return `chg:${pkg}:${version}:${n}`;
}
function docId(file, slug) {
  return `doc:${posix(file)}${slug ? `#${slug}` : ''}`;
}
function prId(repo, n) {
  return `pr:${repo}#${n}`;
}
function commitId(repo, sha) {
  return `commit:${repo}:${sha}`;
}
function epId(method, route) {
  return `ep:${method.toUpperCase()} ${route}`;
}
function tableId(schema, name) {
  return `table:${schema}.${name}`;
}
function migId(file, cls) {
  return `mig:${posix(file)}${cls ? `#${cls}` : ''}`;
}
function custId(slug) {
  return `Cust:${slug}`;
}
function teamId(slug) {
  return `Team:${slug}`;
}
function relId(version) {
  return `Rel:${version}`;
}

/** Tracker keys as the graph writes them: `GROK-n` as is, `#n` / `public-n` as `gh:public#n`. */
function ticketId(key) {
  const gh = /^(?:#|public-)(\d+)$/.exec(key);
  return gh ? `gh:public#${gh[1]}` : key;
}
function parseId(id) {
  if (JIRA_KEY.test(id)) return {
    form: 'tracker',
    tracker: 'jira',
    local: id
  };
  if (GITHUB_KEY.test(id)) return {
    form: 'tracker',
    tracker: 'github',
    scheme: 'gh',
    local: id.slice(3),
    anchor: id.slice(id.indexOf('#') + 1)
  };
  const prefixed = PREFIXED_ID.exec(id);
  if (prefixed) return {
    form: 'prefix',
    prefix: prefixed[1],
    local: prefixed[2]
  };
  const schemed = SCHEMED_ID.exec(id);
  if (!schemed) return {
    form: 'bare',
    local: id
  };
  const [, scheme, local] = schemed;
  const hash = local.indexOf('#');
  const anchor = hash < 0 ? undefined : local.slice(hash + 1);
  const head = hash < 0 ? local : local.slice(0, hash);
  if (PATH_SCHEMES.includes(scheme)) return {
    form: 'scheme',
    scheme,
    local,
    path: head,
    anchor
  };
  return {
    form: 'scheme',
    scheme,
    local,
    anchor,
    parts: head.split(':')
  };
}

/** A readable name for a node that exists only as a reference target. */
function stubName(id) {
  const p = parseId(id);
  if (p.form === 'tracker') return id;
  if (p.anchor !== undefined && p.scheme !== 'gh') return p.anchor;
  if (p.path !== undefined) return path.posix.basename(p.path);
  if (p.parts) return p.parts[p.parts.length - 1];
  return titleCase(p.local.split('/').pop());
}

/** `scatter-plot` -> `Scatter plot`. */
function titleCase(segment) {
  const words = segment.replace(/[-_]+/g, ' ');
  return words.charAt(0).toUpperCase() + words.slice(1);
}
function languageOf(file) {
  return LANGUAGES[path.posix.extname(posix(file)).toLowerCase()] ?? 'other';
}

/** `source_layer` of a node by the location of its file (build-plan.md Decisions "Layers"). */
function sourceLayerOf(file) {
  const p = posix(file);
  return p.startsWith('core/') ? 'core' : p.startsWith('infra/') ? 'infra' : 'public';
}

/** Visibility of a source or evidence path: public/ -> public, core/ -> dev, the internal folder -> internal. */
function locationVisibility(file) {
  const p = posix(file);
  if (p.startsWith('core/docs/knowledge-graph/internal/')) return 'internal';
  return p.startsWith('public/') || p.startsWith('landing/') ? 'public' : 'dev';
}

/** `doc-page.kind` by folder (build-plan.md WO-3c). */
function docKind(file) {
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