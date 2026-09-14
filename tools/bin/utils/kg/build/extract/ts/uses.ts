/// JS API usage per source file (build-plan.md WO-3b): the `DG.X`, `ui.x` and `grok.a.b` tokens of the
/// comment-stripped text, resolved through the exported JS API declarations into `uses` edges with a kind
/// and a count. `UsesLayer` is shared with the samples of WO-3c, where the sample is the `from`.
import * as fs from 'fs';
import * as path from 'path';
import {Emitter} from '../../emitter';
import {BuildContext, Extractor} from '../../registry';
import {fileId} from '../../ids';
import {tsSources, ApiEntry, TsSources} from './declarations';

const API_TOKEN = /\b(DG|ui|grok)((?:\.[A-Za-z_$][\w$]*){1,3})/g;
const COMMENT_OR_STRING = /\/\/.*|\/\*[\s\S]*?\*\/|(["'`])(?:\\[\s\S]|(?!\1)[^\\])*\1/g;
/** The declaration a `DG.X` token means when several exported ones share the name. */
const DG_PREFERENCE = ['class', 'enum', 'const', 'function', 'interface', 'type'];
const UPPER_CASE = /^[A-Z][A-Z0-9_]*$/;
const TOKEN_CONFIDENCE = 0.8;

interface Use {
  to: string;
  kind: string;
}

/** Qualified JS API tokens of [text] with their counts; comments are skipped, string contents are not. */
export function apiTokens(text: string): Map<string, number> {
  const stripped = text.replace(COMMENT_OR_STRING, (m, quote) => quote ? m : '');
  const counts = new Map<string, number>();
  for (const m of stripped.matchAll(API_TOKEN)) counts.set(m[0], (counts.get(m[0]) ?? 0) + 1);
  return counts;
}

export class UsesLayer {
  private api: Map<string, ApiEntry[]>;
  private unresolved = new Map<string, {files: number, first: string}>();

  constructor(private emitter: Emitter, sources: TsSources) {
    this.api = sources.apiIndex();
  }

  /** One `uses` edge per resolved declaration of [text], from [from] (a file or a sample). */
  emit(from: string, text: string, evidence: string): void {
    const uses = new Map<string, {kind: string, count: number}>();
    for (const [token, count] of apiTokens(text)) {
      const use = this.resolve(token);
      if (!use) {
        const seen = this.unresolved.get(token);
        if (seen) seen.files++;
        else this.unresolved.set(token, {files: 1, first: evidence});
        continue;
      }
      const cur = uses.get(use.to);
      if (cur) cur.count += count;
      else uses.set(use.to, {kind: use.kind, count});
    }
    // the scan counts API-looking text, so a shadowed `DG` or a mention inside a string reaches here too; the binder-aware pass raises this to 1
    for (const [to, {kind, count}] of uses)
      this.emitter.edge({type: 'uses', from, to, kind, count, derived_by: 'ast', confidence: TOKEN_CONFIDENCE, evidence: [evidence]});
  }

  /** The unresolved tokens as problems, one per token with the number of files it appears in. */
  finish(): void {
    for (const [token, {files, first}] of [...this.unresolved].sort(([a], [b]) => a < b ? -1 : 1))
      this.emitter.problem('unresolved_ids', `uses: ${token} names no exported JS API declaration (${files} file${files === 1 ? '' : 's'}, first ${first})`);
  }

  private resolve(token: string): Use | undefined {
    const [root, ...segments] = token.split('.');
    if (root === 'DG') {
      const top = this.entry(segments[0], (e) => !e.decl.container);
      const nested = segments.length > 1 && (!top || top.decl.kind === 'const') ? this.entry(`${segments[0]}.${segments[1]}`) : undefined;
      const entry = nested ?? top;
      return entry && {to: entry.id, kind: useKind(entry)};
    }
    if (root === 'ui') {
      const inUi = (e: ApiEntry) => e.file.path.endsWith('/ui.ts');
      const entry = (segments.length > 1 ? this.entry(`${segments[0]}.${segments[1]}`, inUi) : undefined) ?? this.entry(segments[0], inUi);
      return entry && {to: entry.id, kind: 'ui'};
    }
    const space = this.entry(segments[0], (e) => e.file.path.endsWith('/grok.ts') && !e.decl.container);
    const owner = space?.decl.alias;
    if (owner) {
      const entry = segments.length > 1 ? this.entry(`${owner}.${segments[1]}`) : this.entry(owner, (e) => e.decl.kind === 'class');
      return entry && {to: entry.id, kind: useKind(entry)};
    }
    const member = segments.length > 1 ? this.entry(`${segments[0]}.${segments[1]}`) : undefined;
    const entry = member ?? space ?? this.entry(segments[0], (e) => !e.decl.container);
    return entry && {to: entry.id, kind: useKind(entry)};
  }

  private entry(name: string, filter: (e: ApiEntry) => boolean = () => true): ApiEntry | undefined {
    const candidates = this.api.get(name)?.filter(filter);
    if (!candidates?.length) return undefined;
    const rank = (e: ApiEntry) => {
      const i = DG_PREFERENCE.indexOf(e.decl.kind);
      return i < 0 ? DG_PREFERENCE.length : i;
    };
    return candidates.reduce((best, e) => rank(e) < rank(best) ? e : best);
  }
}

/** uses.yaml `kind`: what sort of API member the target is. */
function useKind(entry: ApiEntry): string {
  const {kind, name} = entry.decl;
  if (kind === 'interface') return 'type';
  if (kind === 'class' || kind === 'enum' || kind === 'type') return kind;
  if (kind === 'const') return UPPER_CASE.test(name.slice(name.lastIndexOf('.') + 1)) ? 'enum' : 'function';
  return entry.file.path.endsWith('/ui.ts') ? 'ui' : 'function';
}

export const usesExtractor: Extractor = {
  name: 'ts-uses',
  layer: 'public',
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const sources = tsSources(ctx, emitter);
    const layer = new UsesLayer(emitter, sources);
    for (const file of sources.files)
      layer.emit(fileId(file.path), fs.readFileSync(path.join(ctx.repoRoot, file.path), 'utf8'), file.path);
    layer.finish();
    emitter.source('ts-uses', sources.failed ? 'partial' : 'ok');
  },
};
