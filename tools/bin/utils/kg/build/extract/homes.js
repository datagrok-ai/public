"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.firstParagraph = firstParagraph;
exports.homesExtractor = void 0;
var fs = _interopRequireWildcard(require("fs"));
var path = _interopRequireWildcard(require("path"));
var _glob = require("glob");
var _types = require("../../types");
var _homes = require("../../homes");
var _citations = require("../../citations");
var _frontmatter = require("../../frontmatter");
var _ids = require("../ids");
function _interopRequireWildcard(e, t) { if ("function" == typeof WeakMap) var r = new WeakMap(), n = new WeakMap(); return (_interopRequireWildcard = function (e, t) { if (!t && e && e.__esModule) return e; var o, i, f = { __proto__: null, default: e }; if (null === e || "object" != typeof e && "function" != typeof e) return f; if (o = t ? n : r) { if (o.has(e)) return o.get(e); o.set(e, f); } for (const t in e) "default" !== t && {}.hasOwnProperty.call(e, t) && ((i = (o = Object.defineProperty) && Object.getOwnPropertyDescriptor(e, t)) && (i.get || i.set) ? o(f, t, i) : f[t] = e[t]); return f; })(e, t); }
/// The home-document layer (build-plan.md WO-2): one node per home with its members, the edges its
/// frontmatter keys spell, `code:` roots and body citations as ownership claims for WO-4, and the
/// stubs the references need (tickets, declarations, cited pages).

const homesExtractor = exports.homesExtractor = {
  name: 'homes',
  layer: 'home-document',
  modes: ['full', 'public'],
  run(ctx, emitter) {
    const homes = ctx.homes ?? (0, _homes.loadHomes)(ctx.system, ctx.repoRoot);
    const layer = new HomeLayer(ctx.system, ctx.repoRoot, homes, emitter);
    for (const page of homes.pages) layer.emitPage(page);
    for (const home of homes.homes) layer.emitHome(home);
    if (homes.errors.length) emitter.source('homes', 'partial');
  }
};

/** What an edge key hangs off: a home, or a page that only annotates. */

class HomeLayer {
  byId = new Map();
  byAlias = new Map();
  constructor(system, repoRoot, homes, emitter) {
    this.system = system;
    this.repoRoot = repoRoot;
    this.emitter = emitter;
    for (const home of homes.homes) {
      this.byId.set(home.id, home);
      for (const alias of home.aliases) this.byAlias.set(alias, home);
    }
  }
  emitHome(home) {
    const {
      type,
      data
    } = home;
    const row = {
      type: type.name,
      id: home.id,
      name: home.name,
      home: home.file,
      source_layer: (0, _ids.sourceLayerOf)(home.file),
      provenance: 'annotation'
    };
    for (const [key, value] of Object.entries(data)) {
      if (value === null || key === 'feature' || key === 'id' || key === 'type' || key === 'aliases') continue;
      const member = type.members[key];
      if (!member || (0, _homes.isHomeFileMember)(member) && !home.yaml) continue;
      row[key] = member.kind === 'ref' ? this.resolveMember(value, member) : value;
    }
    for (const member of Object.values(type.members)) if ((0, _homes.isHomeFileMember)(member) && !home.yaml) row[member.name] = home.file;
    if (home.aliases.length) row.aliases = home.aliases;
    if (type.members.manual_only && row.manual_only === undefined && data.target_layer === 'manual-only') row.manual_only = true;
    if (row.description === undefined && !home.yaml) row.description = firstParagraph(home.body);
    this.emitter.node(row);
    const subject = {
      id: home.id,
      type,
      file: home.file,
      fm: home.fm
    };
    const claimed = this.emitKeys(subject, data);
    if (!home.yaml) this.emitCitations(subject, home.body, claimed);
  }
  emitPage(page) {
    const type = this.system.nodes.get('doc-page');
    if (!type) return;
    const data = page.fm.data;
    const id = (0, _ids.docId)(page.file);
    this.emitter.stub(id, type.name, typeof data.title === 'string' ? data.title : path.posix.basename(page.file), 'annotation', {
      path: page.file,
      kind: (0, _ids.docKind)(page.file)
    });
    this.emitKeys({
      id,
      type,
      file: page.file,
      fm: page.fm
    }, data);
  }

  /** Every edge key of the subject; returns the files its `code:` roots claim. */
  emitKeys(subject, data) {
    const claimed = new Set();
    for (const [key, value] of Object.entries(data)) {
      if (value === null || !Array.isArray(value)) continue;
      if (key === 'edges') {
        for (const item of value) {
          if (!item || typeof item !== 'object') continue;
          const {
            type,
            to,
            ...props
          } = item;
          const edge = this.system.edges.get(String(type));
          if (edge && typeof to === 'string') this.emitEdge(subject, edge, 'from', to, props);
        }
        continue;
      }
      const edge = this.system.keys.get(key);
      if (!edge || edge.abstract) continue;
      const targetKey = key === 'code' ? 'path' : 'to';
      for (const item of value) {
        const target = typeof item === 'string' ? item : item && typeof item === 'object' ? item[targetKey] : undefined;
        if (typeof target !== 'string' || !target.trim()) continue;
        const props = typeof item === 'object' ? Object.fromEntries(Object.entries(item).filter(([k]) => k !== targetKey)) : {};
        if (key === 'code') this.emitCode(subject, edge, target.trim(), props, claimed);else this.emitEdge(subject, edge, edge.keySide, target, props);
      }
    }
    return claimed;
  }
  emitEdge(subject, edge, subjectSide, target, props) {
    const otherSide = subjectSide === 'from' ? 'to' : 'from';
    const other = this.resolve(target, edge[otherSide]);
    this.stubFor(other);
    const [from, to] = subjectSide === 'from' ? [subject.id, other] : [other, subject.id];
    this.emitter.edge({
      type: edge.name,
      from,
      to,
      derived_by: 'annotation',
      confidence: 1,
      evidence: [subject.file],
      ...props
    });
  }

  /** A `code:` root: `path#Anchor` is the declaration itself; anything else expands to files and claims them (rung 2). */
  emitCode(subject, edge, target, props, claimed) {
    const hash = target.indexOf('#');
    const file = repoPath(hash < 0 ? target : target.slice(0, hash));
    if (hash >= 0) {
      const id = (0, _ids.declId)(file, target.slice(hash + 1));
      this.emitter.stub(id, 'declaration', target.slice(hash + 1), 'annotation', {
        language: (0, _ids.languageOf)(file),
        path: file
      });
      this.emitter.edge({
        type: edge.name,
        from: subject.id,
        to: id,
        derived_by: 'annotation',
        confidence: 1,
        evidence: [subject.file],
        ...props
      });
      return;
    }
    const line = (0, _frontmatter.keyLine)(subject.fm, 'code');
    for (const p of this.expandRoot(file)) {
      this.emitter.node({
        type: 'source-file',
        id: (0, _ids.fileId)(p),
        name: path.posix.basename(p),
        path: p,
        loc: countLines(path.join(this.repoRoot, p)),
        language: (0, _ids.languageOf)(p),
        provenance: 'filesystem',
        source_layer: (0, _ids.sourceLayerOf)(p)
      });
      this.emitter.claim({
        file: p,
        feature: subject.id,
        rung: 2,
        source: 'home',
        props,
        line
      });
      claimed.add(p);
    }
  }
  expandRoot(p) {
    const clean = p.replace(/\/+$/, '');
    if (_homes.GLOB_MAGIC.test(clean)) return (0, _glob.globSync)(clean, {
      cwd: this.repoRoot,
      ignore: _homes.HOME_IGNORE,
      nodir: true,
      posix: true,
      windowsPathsNoEscape: true
    }).sort();
    const full = path.join(this.repoRoot, clean);
    if (!fs.existsSync(full)) return [];
    if (fs.statSync(full).isDirectory()) return (0, _glob.globSync)(`${clean}/**`, {
      cwd: this.repoRoot,
      ignore: _homes.HOME_IGNORE,
      nodir: true,
      posix: true,
      windowsPathsNoEscape: true
    }).sort();
    return [clean];
  }

  /** Body citations: implementation files become claims (rung 3), documents become mentions of `doc:` stubs. */
  emitCitations(subject, body, claimed) {
    const isFeature = subject.type.root === 'feature';
    const mentions = this.system.edges.get('mentions');
    for (const c of (0, _citations.extractCitations)(subject.file, body, subject.fm.bodyLine)) {
      if (c.resolved === null || c.resolved === subject.file) continue;
      const resolved = this.resolveDocLink(c);
      if (!fs.existsSync(path.join(this.repoRoot, resolved))) continue;
      if (/\.mdx?$/i.test(resolved)) {
        if (!mentions || !this.system.nodes.has('doc-page')) continue;
        const from = mentions.from.some(t => (0, _types.isSubtype)(this.system, subject.type.name, t)) ? subject.id : (0, _ids.docId)(subject.file);
        if (from !== subject.id) this.emitter.stub(from, 'doc-page', path.posix.basename(subject.file), 'annotation', {
          path: subject.file,
          kind: (0, _ids.docKind)(subject.file)
        });
        this.emitter.stub((0, _ids.docId)(resolved), 'doc-page', path.posix.basename(resolved), 'annotation', {
          path: resolved,
          kind: (0, _ids.docKind)(resolved)
        });
        this.emitter.edge({
          type: 'mentions',
          from,
          to: (0, _ids.docId)(resolved),
          derived_by: 'annotation',
          confidence: 1,
          evidence: [subject.file]
        });
        continue;
      }
      if (!isFeature || claimed.has(resolved) || [...claimed].some(f => f.startsWith(`${resolved}/`))) continue;
      this.emitter.claim({
        file: resolved,
        feature: subject.id,
        rung: 3,
        source: 'home',
        props: {},
        line: c.line
      });
      claimed.add(resolved);
    }
  }

  /** A Docusaurus link may drop the extension: `tile-viewer` means `tile-viewer.md` beside the page. */
  resolveDocLink(c) {
    const p = c.resolved;
    if (c.kind === 'backtick' || path.posix.extname(p) || fs.existsSync(path.join(this.repoRoot, p))) return p;
    return [`${p}.md`, `${p}.mdx`, `${p}/index.md`].find(a => fs.existsSync(path.join(this.repoRoot, a))) ?? p;
  }
  resolveMember(value, member) {
    if (member.list) return (Array.isArray(value) ? value : [value]).map(v => typeof v === 'string' ? this.resolve(v, member.refs) : v);
    return typeof value === 'string' ? this.resolve(value, member.refs) : value;
  }

  /** The canonical id of a reference: aliases resolve to the home's id, a bare id gets the prefix of the type it lands in. */
  resolve(value, expected) {
    const raw = value.trim().replace(/^~/, '');
    if (_ids.JIRA_KEY.test(raw) || /^#\d+$/.test(raw)) return (0, _ids.ticketId)(raw);
    if (_ids.PREFIXED_ID.test(raw)) {
      const {
        id
      } = (0, _homes.splitRefAnchor)(raw);
      return this.lookup([id]) ?? id;
    }
    if (_ids.SCHEMED_ID.test(raw)) return raw;
    const {
      id
    } = (0, _homes.splitRefAnchor)(raw);
    const candidates = [...new Set((0, _types.concreteAuthored)(this.system, expected).map(t => t.prefix ? `${t.prefix}:${id}` : id))];
    return this.lookup(candidates) ?? candidates[0] ?? id;
  }
  lookup(candidates) {
    for (const c of candidates) {
      const home = this.byId.get(c) ?? this.byAlias.get(c);
      if (home) return home.id;
    }
    return undefined;
  }

  /** The stub an extracted reference target needs when no extractor has produced the node. */
  stubFor(id) {
    const parsed = (0, _ids.parseId)(id);
    if (parsed.form === 'tracker') this.emitter.stub(id, 'ticket', id, 'annotation', {
      tracker: parsed.tracker,
      key: id,
      kind: 'unknown',
      state: 'open'
    });else if (parsed.scheme === 'decl' && parsed.path) this.emitter.stub(id, 'declaration', parsed.anchor ?? path.posix.basename(parsed.path), 'annotation', {
      language: (0, _ids.languageOf)(parsed.path),
      path: parsed.path
    });else if (parsed.scheme === 'doc' && parsed.path) this.emitter.stub(id, 'doc-page', path.posix.basename(parsed.path), 'annotation', {
      path: parsed.path,
      kind: (0, _ids.docKind)(parsed.path)
    });
  }
}

/** The first paragraph of prose after the frontmatter: not a heading, table, comment or import line. */
function firstParagraph(body) {
  const lines = (0, _citations.proseLines)(body).map(l => l.text.trim());
  const start = lines.findIndex(l => l && !/^(#|\||<!--|---|import\s|:::)/.test(l));
  if (start < 0) return undefined;
  let end = start;
  while (end < lines.length && lines[end] && !/^(#|\||<!--|:::)/.test(lines[end])) end++;
  return lines.slice(start, end).join(' ');
}

/** `landing:` and `infra:` prefixed roots name those repos' folders. */
function repoPath(p) {
  const repo = _homes.REPO_PREFIX.exec(p.trim());
  return repo ? `${repo[1]}/${repo[2]}` : p.trim();
}
function countLines(file) {
  const buffer = fs.readFileSync(file);
  if (!buffer.length) return 0;
  let n = 0;
  for (let at = buffer.indexOf(0x0A); at >= 0; at = buffer.indexOf(0x0A, at + 1)) n++;
  return buffer[buffer.length - 1] === 0x0A ? n : n + 1;
}