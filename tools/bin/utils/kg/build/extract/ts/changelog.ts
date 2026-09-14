/// Package changelogs (build-plan.md WO-3c): one `changelog-entry` per bullet under a `## <version> (<date>)` or
/// `## v.next` heading, the tickets it names as mentions, and a `~id` that resolves to a feature as `changes`.
import * as fs from 'fs';
import * as path from 'path';
import {Emitter} from '../../emitter';
import {Row} from '../../normalize';
import {BuildContext, Extractor} from '../../registry';
import {pkgId, chgId} from '../../ids';
import {HomeIndex, homesOf, emitMentions} from '../markers';
import {listPackages} from './packages';

export interface ChangelogEntry {
  version: string;
  date?: string;
  /** 1-based position under the version heading. */
  n: number;
  text: string;
}

const HEADING = /^##\s+(.+?)\s*$/;
const BULLET = /^\s*[*+-]\s+(.*\S)\s*$/;
const VERB = /\b(add(?:ed|s)?|fix(?:e[ds])?|improve[ds]?|remove[ds]?)\b/i;
/** `changes.kind` by the first three letters of the verb. */
const KINDS: Record<string, string> = {add: 'added', fix: 'fixed', imp: 'improved', rem: 'removed'};
const NAME_CAP = 100;

export const changelogExtractor: Extractor = {
  name: 'ts-changelog',
  layer: 'public',
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const index = new HomeIndex(homesOf(ctx));
    let unresolved = 0;
    for (const pkg of listPackages(ctx.repoRoot)) {
      const file = `${pkg.dir}/CHANGELOG.md`;
      const full = path.join(ctx.repoRoot, file);
      if (!fs.existsSync(full)) continue;
      for (const e of parseChangelog(fs.readFileSync(full, 'utf8'))) {
        const id = chgId(pkg.folder, e.version, e.n);
        emitter.node({type: 'changelog-entry', id, name: e.text.length > NAME_CAP ? `${e.text.slice(0, NAME_CAP - 1)}…` : e.text, text: e.text, package: pkgId(pkg.folder),
          version: e.version, date: e.date, path: file, provenance: 'annotation', source_layer: 'public'});
        const kind = KINDS[VERB.exec(e.text)?.[1].slice(0, 3).toLowerCase() ?? ''];
        const changes = (target: {root?: string}): Row | undefined => target.root === 'feature' ? {type: 'changes', kind} : undefined;
        unresolved += emitMentions(emitter, id, e.text, `${file} (${e.version} #${e.n})`, index, changes).unresolved.length;
      }
    }
    emitter.source('ts-changelog', unresolved ? 'partial' : 'ok');
  },
};

/** The bullets of a changelog under their version headings: `## 1.2.3 (2026-01-31)` and `## v.next` (version `0`, no date);
 * a heading without a version number closes the section. Continuation lines of a bullet are joined with a space. */
export function parseChangelog(text: string): ChangelogEntry[] {
  const out: ChangelogEntry[] = [];
  const counters = new Map<string, number>();
  let current: {version: string, date?: string} | null = null;
  let last: ChangelogEntry | null = null;
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.replace(/\s+$/, '');
    const heading = HEADING.exec(line);
    if (heading || /^#/.test(line)) {
      last = null;
      if (!heading) continue;
      const version = /^v\.next\b/i.test(heading[1]) ? '0' : /\d+(?:\.\d+)*/.exec(heading[1])?.[0];
      current = version === undefined ? null : {version, date: calendarDate(/\((\d{4}-\d{2}-\d{2})/.exec(heading[1])?.[1])};
      continue;
    }
    if (!current) continue;
    const bullet = BULLET.exec(line);
    if (bullet) {
      const n = (counters.get(current.version) ?? 0) + 1;
      counters.set(current.version, n);
      last = {version: current.version, date: current.date, n, text: bullet[1]};
      out.push(last);
    }
    else if (last && /^\s+\S/.test(line)) last.text += ` ${line.trim()}`;
    else last = null;
  }
  return out;
}

/** Changelogs carry swapped days and months (`2025-27-01`); a heading date that is no calendar date is no date. */
function calendarDate(date: string | undefined): string | undefined {
  const parsed = date ? new Date(`${date}T00:00:00Z`) : undefined;
  return parsed && !Number.isNaN(parsed.getTime()) && parsed.toISOString().startsWith(date!) ? date : undefined;
}
