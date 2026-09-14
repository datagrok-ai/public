/// API samples (build-plan.md WO-3c): every script under packages/ApiSamples/scripts as a `sample` with the API
/// members its `//api:` header names, else the `DG.` / `ui.` / `grok.` members its text uses (js-api/scripts/inventory.cjs),
/// the help page a `//help-url:` names as a mention, a `demonstrates` edge for a `~id` in the header, and the
/// `uses` edges of the JS API members it calls, resolved the way WO-3b resolves them for a source file.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {Emitter} from '../../emitter';
import {BuildContext, Extractor} from '../../registry';
import {sampleId, docId, docKind, languageOf} from '../../ids';
import {commentPrefix} from '../../annotations';
import {HomeIndex, homesOf, idTokens, resolveMention} from '../markers';
import {tsSources} from './declarations';
import {UsesLayer} from './uses';

const SCRIPTS_DIR = 'public/packages/ApiSamples/scripts';
const HELP_DIR = 'public/help';
const API_MEMBER = /\b((?:DG|ui|grok)(?:\.[A-Za-z_$][\w$]*){1,3})/g;
const HELP_URL = /^(?:https?:\/\/(?:[\w-]+\.)*datagrok\.ai)?\/help\/([^#?]+?)(?:\.mdx?)?\/?(?:[#?].*)?$/;

export const samplesExtractor: Extractor = {
  name: 'ts-samples',
  layer: 'public',
  modes: ['full', 'public'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const files = globSync(`${SCRIPTS_DIR}/**/*.{js,py,R,r}`, {cwd: ctx.repoRoot, ignore: ['**/node_modules/**'], nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
    const uses = new UsesLayer(emitter, tsSources(ctx, emitter));
    const index = new HomeIndex(homesOf(ctx));
    let missing = 0;
    for (const file of files) {
      const rel = file.slice(SCRIPTS_DIR.length + 1);
      const id = sampleId(rel);
      const text = fs.readFileSync(path.join(ctx.repoRoot, file), 'utf8');
      const language = languageOf(file);
      const header = parseSampleHeader(text, commentPrefix(language));
      const members = header.keys.api ? header.keys.api.split(',').map((s) => s.trim()).filter(Boolean) : [...new Set([...text.matchAll(API_MEMBER)].map((m) => m[1]))];
      const folder = path.posix.dirname(rel);
      if (!emitter.node({type: 'sample', id, name: path.posix.basename(rel).replace(/\.[^.]+$/, ''), description: header.description, path: file, language, folder: folder === '.' ? '' : folder,
        api_members: members.length ? members : undefined, provenance: Object.keys(header.keys).length ? 'annotation' : 'ast', source_layer: 'public'}).accepted) continue;
      if (header.keys['help-url']) {
        const page = helpPage(ctx.repoRoot, header.keys['help-url']);
        if (page) {
          emitter.stub(docId(page), 'doc-page', path.posix.basename(page), 'annotation', {path: page, kind: docKind(page)});
          emitter.edge({type: 'mentions', from: id, to: docId(page), derived_by: 'annotation', confidence: 1, evidence: [file]});
        }
        else {
          missing++;
          emitter.problem('unresolved_ids', `${file}: help-url ${header.keys['help-url']} names no page under ${HELP_DIR}`);
        }
      }
      for (const token of idTokens(header.text).keys()) {
        const feature = resolveMention(emitter, index, token, file);
        if (feature) emitter.edge({type: 'demonstrates', from: id, to: feature.id, derived_by: 'annotation', confidence: 1, evidence: [file]});
      }
      uses.emit(id, header.body, file);
    }
    uses.finish();
    emitter.source('ts-samples', missing ? 'partial' : 'ok');
  },
};

interface SampleHeader {
  /** `//key: value` lines of the leading comment block. */
  keys: Record<string, string>;
  /** The first prose line of the block. */
  description?: string;
  /** The whole block without comment markers. */
  text: string;
  /** The source after the block: what the sample runs. */
  body: string;
}

/** The leading comment block of a sample: header keys, the first prose line, the text a `~id` may sit in, and the code under it. */
export function parseSampleHeader(source: string, prefix: string): SampleHeader {
  const header: SampleHeader = {keys: {}, text: '', body: ''};
  const lines = source.split(/\r?\n/);
  let i = 0;
  while (i < lines.length && !lines[i].trim()) i++;
  const block: string[] = [];
  for (; i < lines.length && lines[i].trim().startsWith(prefix); i++) {
    const line = lines[i].trim().slice(prefix.length).trim();
    block.push(line);
    const key = /^([a-z][\w-]*):\s*(.*)$/.exec(line);
    if (key && !/^https?$/.test(key[1])) header.keys[key[1]] ??= key[2].trim();
    else if (header.description === undefined && line) header.description = line;
  }
  header.text = block.join('\n');
  header.body = lines.slice(i).join('\n');
  return header;
}

/** The help page a `//help-url:` points at, as a repo path: `.../help/datagrok/project` is `project.md`, `project.mdx`,
 * `project/index.md` or `project/project.md` under public/help, whichever exists. */
export function helpPage(repoRoot: string, url: string): string | undefined {
  const m = HELP_URL.exec(url.trim());
  if (!m) return undefined;
  const p = `${HELP_DIR}/${m[1]}`;
  return [`${p}.md`, `${p}.mdx`, `${p}/index.md`, `${p}/${path.posix.basename(p)}.md`].find((c) => fs.existsSync(path.join(repoRoot, c)));
}
