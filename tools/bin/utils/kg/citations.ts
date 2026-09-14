/// Citations in a markdown body: repo paths in backticks, inline and reference-style links.
/// Extracted once into records so `check` (stale citations) and `build` (participation edges)
/// read the same evidence (conventions.md §5.4).
import * as path from 'path';

export interface Citation {
  /** The citing document, posix path relative to the monorepo root. */
  file: string;
  line: number;
  kind: 'backtick' | 'link' | 'reference-link';
  raw: string;
  /** Repo-relative posix path, or null when the target escapes the repository. */
  resolved: string | null;
  anchor?: string;
  /** A link to documentation (`*.md`, `*.mdx`) or a citation of an implementation file. */
  target: 'doc' | 'code';
}

export interface ProseLine {
  text: string;
  /** 0-based offset from the start of the body. */
  offset: number;
}

const REPO_PATH = /^(core|public|infra|landing)[\\/][A-Za-z0-9_.\\/*-]+$/;
const SCHEME = /^[a-z][a-z0-9+.-]*:/i;
const INLINE_LINK = /!?\[[^\]]*\]\((?:<([^>]+)>|([^)\s]+))(?:\s+"[^"]*")?\)/g;
const REFERENCE_DEFINITION = /^ {0,3}\[([^\]]+)\]:\s*(?:<([^>]+)>|(\S+))/;

/** Body lines outside fenced code blocks. A fence closes only on the same character with at least the opening length. */
export function proseLines(body: string): ProseLine[] {
  const out: ProseLine[] = [];
  let fence: {char: string, length: number} | null = null;
  body.split('\n').forEach((text, offset) => {
    const m = /^\s{0,3}(`{3,}|~{3,})/.exec(text);
    if (m) {
      const run = m[1];
      if (fence === null) {
        fence = {char: run[0], length: run.length};
        return;
      }
      if (run[0] === fence.char && run.length >= fence.length) {
        fence = null;
        return;
      }
    }
    if (fence === null) out.push({text, offset});
  });
  return out;
}

/** The GitHub / Docusaurus slug of a heading text. */
export function slugify(text: string): string {
  return text.toLowerCase().replace(/[^\w\s-]/g, '').trim().replace(/\s+/g, '-');
}

/** Anchors a markdown body answers to: heading slugs plus explicit `{#id}` suffixes. */
export function headingAnchors(body: string): Set<string> {
  const anchors = new Set<string>();
  for (const {text} of proseLines(body)) {
    const m = /^#{1,6}\s+(.+?)\s*(?:\{#([^}]+)\})?\s*#*\s*$/.exec(text);
    if (!m) continue;
    anchors.add(slugify(m[1]));
    if (m[2]) anchors.add(m[2]);
  }
  return anchors;
}

export function extractCitations(file: string, body: string, bodyLine: number): Citation[] {
  const out: Citation[] = [];
  const dir = path.posix.dirname(file);
  for (const {text, offset} of proseLines(body)) {
    const line = bodyLine + offset;
    for (const m of text.matchAll(/`([^`\n]+)`/g)) {
      const token = m[1].trim().replace(/:\d+(-\d+)?$/, '');
      const {target, anchor} = splitAnchor(token);
      if (REPO_PATH.test(target)) out.push(citation(file, line, 'backtick', m[1], target.replace(/\/+$/, ''), anchor));
    }
    for (const m of text.matchAll(INLINE_LINK))
      pushLink(out, file, dir, line, 'link', m[1] ?? m[2]);
    const def = REFERENCE_DEFINITION.exec(text);
    if (def) pushLink(out, file, dir, line, 'reference-link', def[2] ?? def[3]);
  }
  return out;
}

function pushLink(out: Citation[], file: string, dir: string, line: number, kind: Citation['kind'], raw: string): void {
  if (SCHEME.test(raw)) return;
  const {target, anchor} = splitAnchor(decodeURIComponent(raw.trim()));
  if (!target) {
    if (anchor) out.push(citation(file, line, kind, raw, file, anchor));
    return;
  }
  const resolved = target.startsWith('/') ? path.posix.normalize(target.slice(1)) :
    REPO_PATH.test(target) ? path.posix.normalize(target) : path.posix.normalize(path.posix.join(dir, target));
  out.push(citation(file, line, kind, raw, resolved.startsWith('..') || path.posix.isAbsolute(resolved) ? null : resolved.replace(/\/+$/, ''), anchor));
}

function splitAnchor(token: string): {target: string, anchor?: string} {
  const hash = token.indexOf('#');
  return hash < 0 ? {target: token} : {target: token.slice(0, hash), anchor: token.slice(hash + 1) || undefined};
}

function citation(file: string, line: number, kind: Citation['kind'], raw: string, resolved: string | null, anchor?: string): Citation {
  return {file, line, kind, raw, resolved, anchor, target: resolved !== null && /\.mdx?$/i.test(resolved) ? 'doc' : 'code'};
}
