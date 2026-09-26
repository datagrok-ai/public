/// The marketing site (github datagrok-ai/landing, `--landing <dir>`): every `web/**/*.html` page as a `web-page`
/// with the address nginx gives it, its h1–h3 headings as `doc-anchor` nodes, and the media it shows as embeds for
/// the media extractor. Read lexically, like the Dart pass: regexes over the HTML, no DOM.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {Emitter} from '../emitter';
import {BuildContext, Extractor} from '../context';
import {Heading, slugify} from '../../citations';
import {extractEmbeds} from '../../embeds';
import {docId} from '../../ids';
import {LANDING_PREFIX, pageUrl} from '../../roots';

const PAGES = 'web/**/*.html';
/** The site's copy of help uploads and the bare-html snippets are not pages. */
const IGNORE = ['web/help/**', 'web/bare-html/**'];
const TITLE = /<title>([^<]*)<\/title>/i;
const DESCRIPTION = /<meta\s+(?=[^>]*\bname\s*=\s*"description")[^>]*\bcontent\s*=\s*"([^"]*)"/i;
const HEADING = /<h([1-3])\b([^>]*)>([\s\S]*?)<\/h\1>/gi;
const ID_ATTRIBUTE = /\bid\s*=\s*"([^"]*)"/i;
const TAGS = /<[^>]+>/g;
/** A title that names the company rather than the page. */
const GENERIC_TITLE = /^datagrok$/i;

export const landingExtractor: Extractor = {
  name: 'landing',
  describes: {landing: 'the pages of the marketing site, their headings and the media they show'},
  modes: ['full', 'public'],
  run(ctx: BuildContext, emitter: Emitter): void {
    if (ctx.landingDir === undefined) {
      emitter.source('landing', 'missing');
      return;
    }
    for (const rel of globSync(PAGES, {cwd: ctx.landingDir, ignore: IGNORE, nodir: true, posix: true}).sort()) {
      const file = `${LANDING_PREFIX}${rel}`;
      const html = fs.readFileSync(path.join(ctx.landingDir, rel), 'utf8');
      const id = docId(file);
      const heads = htmlHeadings(html);
      const title = TITLE.exec(html)?.[1].trim();
      const name = title && !GENERIC_TITLE.test(title) ? title : heads.find((h) => h.depth === 1)?.text ?? path.posix.basename(rel, '.html');
      const description = DESCRIPTION.exec(html)?.[1].trim() || undefined;
      const admitted = emitter.node({type: 'web-page', id, name, description, path: file, kind: 'marketing', url: pageUrl(file), provenance: 'filesystem', source_layer: 'public'});
      if (!admitted.accepted) continue;
      for (const h of heads)
        emitter.node({type: 'doc-anchor', id: docId(file, h.slug), name: h.text, path: file, page: id, slug: h.slug, depth: h.depth, provenance: 'filesystem', source_layer: 'public'});
      for (const e of extractEmbeds(file, html, 1, heads)) emitter.embed(file, e);
    }
    emitter.source('landing', 'ok');
  },
};

/** The h1–h3 of an HTML page, in order, with the anchor each answers to: its `id` attribute, else the slug of its text, numbered per duplicate. */
export function htmlHeadings(html: string): Heading[] {
  const out: Heading[] = [];
  const seen = new Map<string, number>();
  for (const m of html.matchAll(HEADING)) {
    const text = m[3].replace(TAGS, '').replace(/\s+/g, ' ').trim();
    const explicit = ID_ATTRIBUTE.exec(m[2])?.[1];
    const base = explicit || slugify(text);
    if (!base) continue;
    const n = seen.get(base) ?? 0;
    seen.set(base, n + 1);
    out.push({depth: Number(m[1]), text, slug: n ? `${base}-${n}` : base, line: html.slice(0, m.index).split('\n').length - 1});
  }
  return out;
}
