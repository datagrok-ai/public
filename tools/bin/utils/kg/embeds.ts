/// Embeds: the media a page shows, as records the docs extractor hands to the media extractor. Markdown images
/// (inline, reference-style, a thumbnail wrapped in a YouTube link), <img>, <video>/<source>, <picture>, the
/// Docusaurus <Image> and <Iframe> components, YouTube iframes and watch links, and the {text, image} slides of
/// the marketing pages. Fenced code is skipped; the nearest heading above an embed is its anchor.
import * as path from 'path';
import {proseLines, Heading} from './citations';
import {isMedia} from './media';

export type EmbedForm = 'image' | 'tag' | 'iframe' | 'link' | 'slide';

export type EmbedTarget = {kind: 'file', path: string} | {kind: 'hosted', provider: 'youtube', id: string};

export interface Embed {
  line: number;
  /** Ordinal in the page, from 1. */
  position: number;
  form: EmbedForm;
  target: EmbedTarget;
  raw: string;
  anchor?: string;
  alt?: string;
  title?: string;
  caption?: string;
  start_seconds?: number;
  /** A poster or thumbnail the same occurrence shows for the target (a video's poster, a YouTube link's image). */
  poster?: string;
}

/** Where a site-absolute image path lives in the repositories. */
const ABSOLUTE_ROOTS: [RegExp, string][] = [
  [/^(?:\.\/)?docusaurus_img\//, 'public/docusaurus/static/docusaurus_img/'],
  [/^\/docusaurus_img\//, 'public/docusaurus/static/docusaurus_img/'],
  [/^\/help\//, 'public/help/'],
  [/^\/img\//, 'landing:web/img/'],
];
const YOUTUBE = /(?:https?:)?\/\/(?:www\.|m\.)?(?:youtube\.com\/(?:watch\?(?:[^"'\s)]*&)?v=|embed\/|v\/|shorts\/)|youtu\.be\/)([\w-]{11})(?:[?&][^"'\s)]*)?/;
const START = /[?&](?:t|start)=(\d+)s?/;
const INLINE_IMAGE = /!\[([^\]]*)\]\((?:<([^>]+)>|([^)\s]+))(?:\s+"([^"]*)")?\)/g;
const REFERENCE_IMAGE = /!\[([^\]]*)\](?:\[([^\]]*)\])?(?!\()/g;
const REFERENCE_DEFINITION = /^ {0,3}\[([^\]]+)\]:\s*(?:<([^>]+)>|(\S+))(?:\s+"([^"]*)")?/;
const LINK_AROUND_IMAGE = /\[(!\[[^\]]*\]\([^)]*\))\]\(([^)\s]+)\)/g;
const LINK = /(?<!!)\[([^\]]*)\]\(([^)\s]+)\)/g;
const BARE_URL = /(?<![("'<\]=])https?:\/\/[^\s<>"')\]]+/g;
const A_HREF = /<a\b[^>]*\bhref\s*=\s*"([^"]*)"/g;
const BACKGROUND = /background(?:-image)?\s*:\s*url\(\s*['"]?([^'")]+)['"]?\s*\)/g;
const TAG = /<(img|source|iframe|Image|Iframe)\b([^>]*?)\/?>/g;
const VIDEO = /<video\b([^>]*)>([\s\S]*?)<\/video>/g;
const SOURCE = /<source\b([^>]*?)\/?>/g;
const ATTRIBUTE = /([\w-]+)\s*=\s*(?:"([^"]*)"|'([^']*)'|\{([^}]*)\})/g;
const IMPORT = /^\s*import\s+(\w+)\s+from\s+['"]([^'"]+)['"]/;
const REQUIRE = /require\(\s*['"]([^'"]+)['"]\s*\)/;
const SLIDE = /\{\s*text:\s*(['"`])((?:\\.|(?!\1).)*)\1\s*,\s*image:\s*(['"`])((?:\\.|(?!\3).)*)\3\s*\}/g;
const MAX_ANCHOR_DEPTH = 4;

export function extractEmbeds(file: string, body: string, bodyLine: number, heads: Heading[]): Embed[] {
  const dir = path.posix.dirname(file);
  const lines = proseLines(body);
  const text = lines.map((l) => l.text).join('\n');
  const offsetOf = lineIndex(lines);
  const definitions = new Map<string, {target: string, title?: string}>();
  const imports = new Map<string, string>();
  for (const {text: t} of lines) {
    const def = REFERENCE_DEFINITION.exec(t);
    if (def) definitions.set(def[1].toLowerCase(), {target: def[2] ?? def[3], title: def[4]});
    const imp = IMPORT.exec(t);
    if (imp) imports.set(imp[1], imp[2]);
  }
  const found: (Omit<Embed, 'position' | 'anchor'> & {index: number})[] = [];
  const add = (index: number, form: EmbedForm, raw: string, target: string | undefined, props: Partial<Embed> = {}) => {
    if (!target) return;
    const resolved = resolveTarget(target, dir);
    if (!resolved) return;
    found.push({index, line: bodyLine + offsetOf(index), form, target: resolved, raw, ...props});
  };
  const linked = new Set<number>();
  for (const m of text.matchAll(LINK_AROUND_IMAGE)) {
    const video = YOUTUBE.exec(m[2]);
    if (!video) continue;
    const image = new RegExp(INLINE_IMAGE.source).exec(m[1]);
    const poster = image ? resolveTarget(image[2] ?? image[3], dir) : null;
    linked.add(m.index! + 1);
    add(m.index!, 'link', m[0], m[2], {alt: image?.[1] || undefined, title: image?.[4] || undefined, start_seconds: startOf(m[2]),
      poster: poster?.kind === 'file' ? poster.path : undefined});
  }
  for (const m of text.matchAll(INLINE_IMAGE))
    if (!linked.has(m.index!)) add(m.index!, 'image', m[0], m[2] ?? m[3], {alt: m[1] || undefined, title: m[4] || undefined});
  for (const m of text.matchAll(REFERENCE_IMAGE)) {
    const def = definitions.get((m[2] || m[1]).toLowerCase());
    if (def) add(m.index!, 'image', m[0], def.target, {alt: m[1] || undefined, title: def.title || undefined});
  }
  const consumed: [number, number][] = [];
  for (const m of text.matchAll(VIDEO)) {
    const attrs = attributes(m[1], imports);
    const poster = attrs.poster ? resolveTarget(attrs.poster, dir) : null;
    const sources = [...m[2].matchAll(SOURCE)].map((s) => attributes(s[1], imports).src).filter((s): s is string => !!s);
    for (const src of attrs.src ? [attrs.src, ...sources] : sources)
      add(m.index!, 'tag', m[0], src, {alt: attrs.alt || undefined, title: attrs.title || undefined, poster: poster?.kind === 'file' ? poster.path : undefined});
    consumed.push([m.index!, m.index! + m[0].length]);
  }
  for (const m of text.matchAll(TAG)) {
    if (consumed.some(([from, to]) => m.index! >= from && m.index! < to)) continue;
    const attrs = attributes(m[2], imports);
    const tag = m[1].toLowerCase();
    const src = tag === 'iframe' ? attrs.src ?? attrs.url : attrs.src;
    add(m.index!, tag === 'iframe' ? 'iframe' : 'tag', m[0], src, {alt: attrs.alt || undefined, title: attrs.title || undefined, caption: attrs.caption || undefined,
      start_seconds: src ? startOf(src) : undefined});
  }
  for (const m of text.matchAll(LINK))
    if (YOUTUBE.test(m[2]) && !linked.has(m.index!) && !text.slice(m.index! - 1, m.index!).startsWith('!'))
      add(m.index!, 'link', m[0], m[2], {title: m[1] || undefined, start_seconds: startOf(m[2])});
  for (const m of text.matchAll(A_HREF))
    if (YOUTUBE.test(m[1])) add(m.index!, 'link', m[0], m[1], {start_seconds: startOf(m[1])});
  for (const m of text.matchAll(BACKGROUND))
    add(m.index!, 'tag', m[0], m[1]);
  for (const m of text.matchAll(BARE_URL))
    if (YOUTUBE.test(m[0])) add(m.index!, 'link', m[0], m[0], {start_seconds: startOf(m[0])});
  for (const m of text.matchAll(SLIDE))
    add(m.index!, 'slide', m[0], m[4], {caption: m[2] || undefined});
  found.sort((a, b) => a.index - b.index);
  const anchors = heads.filter((h) => h.depth <= MAX_ANCHOR_DEPTH && h.line !== undefined);
  const out: Embed[] = [];
  const seen = new Set<string>();
  found.forEach((e) => {
    const key = `${e.index}:${e.form}:${targetKey(e.target)}`;
    if (seen.has(key)) return;
    seen.add(key);
    const {index, ...rest} = e;
    const heading = [...anchors].reverse().find((h) => h.line! < e.line - bodyLine);
    out.push({...rest, position: out.length + 1, anchor: heading?.slug});
  });
  return out;
}

export function targetKey(target: EmbedTarget): string {
  return target.kind === 'file' ? target.path : `${target.provider}:${target.id}`;
}

/** A repo path (posix, `landing:` prefixed for the marketing site), or a hosted video; null for anything else. */
export function resolveTarget(raw: string, dir: string): EmbedTarget | null {
  const value = raw.trim().replace(/^<|>$/g, '');
  const video = YOUTUBE.exec(value);
  if (video) return {kind: 'hosted', provider: 'youtube', id: video[1]};
  if (/^[a-z][a-z0-9+.-]*:/i.test(value) && !/^landing:/.test(value)) return null;
  let target = decodeURIComponent(value.split('#')[0].split('?')[0]);
  if (!target) return null;
  for (const [pattern, root] of ABSOLUTE_ROOTS)
    if (pattern.test(target)) {
      target = root + target.replace(pattern, '');
      break;
    }
  const prefixed = /^landing:(.+)$/.exec(target);
  if (prefixed) return isMedia(prefixed[1]) ? {kind: 'file', path: `landing:${path.posix.normalize(prefixed[1])}`} : null;
  const resolved = target.startsWith('/') ? path.posix.normalize(target.slice(1)) :
    /^(core|public|infra)\//.test(target) ? path.posix.normalize(target) : path.posix.normalize(path.posix.join(dir, target));
  if (resolved.startsWith('..') || path.posix.isAbsolute(resolved) || !isMedia(resolved)) return null;
  return {kind: 'file', path: resolved};
}

function startOf(url: string): number | undefined {
  const m = START.exec(url);
  return m ? Number(m[1]) : undefined;
}

/** Tag attributes; `src={Name}` through the page's imports, `src={require('x').default}` through the call. */
function attributes(text: string, imports: Map<string, string>): Record<string, string> {
  const out: Record<string, string> = {};
  for (const m of text.matchAll(ATTRIBUTE)) {
    const name = m[1].toLowerCase();
    if (m[4] !== undefined) {
      const req = REQUIRE.exec(m[4]);
      const value = req ? req[1] : imports.get(m[4].trim());
      if (value) out[name] = value;
    } else out[name] = m[2] ?? m[3];
  }
  return out;
}

/** Maps a character index in the joined prose back to the 0-based body offset of its line. */
function lineIndex(lines: {text: string, offset: number}[]): (index: number) => number {
  const starts: number[] = [];
  let at = 0;
  for (const l of lines) {
    starts.push(at);
    at += l.text.length + 1;
  }
  return (index) => {
    let lo = 0;
    let hi = starts.length - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (starts[mid] <= index) lo = mid;
      else hi = mid - 1;
    }
    return lines[lo]?.offset ?? 0;
  };
}
