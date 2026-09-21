/// Media files: what a page shows rather than cites. The one place that says which extensions are
/// media, reads the records people and the enrichment pass write beside them (`media.yaml` per
/// folder, `videos.yaml` for hosted videos), and finds the git blob id of a file without reading it.
/// `check` (a missing image, a bad record), the home extractor (no ownership claim, no source-file)
/// and the media extractor read the same rules.
import * as fs from 'fs';
import * as path from 'path';
import {createHash} from 'crypto';
import {spawnSync} from 'child_process';
import {Issue, TypeSystem} from './types';
import {normalizeRow, Row} from './normalize';
import {parseYamlDocument, keyLine} from './frontmatter';
import {mediaId, hostedMediaId, posix} from './ids';

export type MediaFormat = 'png' | 'jpg' | 'gif' | 'svg' | 'webp' | 'mp4' | 'webm' | 'pdf' | 'youtube' | 'other';

const FORMATS: Record<string, MediaFormat> = {
  '.png': 'png', '.jpg': 'jpg', '.jpeg': 'jpg', '.gif': 'gif', '.svg': 'svg', '.webp': 'webp', '.mp4': 'mp4', '.webm': 'webm',
  '.pdf': 'pdf', '.ico': 'other', '.mov': 'other', '.avif': 'other', '.bmp': 'other',
};

export const MEDIA_EXT = Object.keys(FORMATS);
/** The record files: one per folder of media files, one catalog of hosted videos (any folder may carry either). */
export const RECORD_FILES = ['media.yaml', 'videos.yaml'];
/** What a record may say about a media node; everything else on the node comes from the file, the page or the build. */
export const RECORD_KEYS = ['kind', 'animated', 'width', 'height', 'seconds', 'caption', 'description', 'actions', 'ui_text',
  'quality', 'quality_notes', 'reviewed', 'described_by', 'described_blob', 'illustrates'];
/** A hosted video has no file name to take its name from, so its record may carry a title. */
const HOSTED_KEYS = [...RECORD_KEYS, 'title'];
const HOSTED_KEY = /^(youtube):([\w-]{11})$/;
export const ILLUSTRATES_TARGETS = ['feature', 'concept'];

export function isMedia(file: string): boolean {
  return path.posix.extname(file).toLowerCase() in FORMATS;
}

export function formatOf(file: string): MediaFormat {
  return FORMATS[path.posix.extname(file).toLowerCase()] ?? 'other';
}

/** One entry of a record file, validated: the node members it sets and the features it says the media illustrates. */
export interface MediaRecord {
  id: string;
  /** Repo path of the file; absent for a hosted video. */
  path?: string;
  provider?: string;
  externalId?: string;
  /** The record file and the line of the entry's key. */
  file: string;
  line: number;
  data: Row;
  illustrates: string[];
}

export interface MediaRecords {
  records: Map<string, MediaRecord>;
  errors: Issue[];
  warnings: Issue[];
}

export type ResolveRef = (value: string, expected: string[], ctx: {source: string, key: string}) => string | undefined;

/**
 * Reads the record files: a key is a file beside the record (`media.yaml`) or `youtube:<id>` (`videos.yaml`); the
 * entry's keys are the authorable members plus `illustrates`, checked against the media type through [system] and
 * the reference resolver of the homes; a `described_blob` that is not the file's current blob is a stale description.
 */
export function loadMediaRecords(system: TypeSystem, repoRoot: string, files: string[], resolve: ResolveRef): MediaRecords {
  const out: MediaRecords = {records: new Map(), errors: [], warnings: []};
  if (!system.nodes.has('media')) return out;
  let blobs: Map<string, string> | undefined;
  for (const file of files) {
    const hosted = path.posix.basename(file) === 'videos.yaml';
    const dir = path.posix.dirname(file);
    let fm;
    try {
      fm = parseYamlDocument(fs.readFileSync(path.join(repoRoot, file), 'utf8'));
    } catch (e: any) {
      out.errors.push({file, code: 'unreadable', message: `cannot read: ${e.message}`});
      continue;
    }
    if (fm.error || !fm.data) {
      out.errors.push({file, line: fm.errorLine, code: 'yaml-error', message: fm.error ?? 'empty record file'});
      continue;
    }
    for (const [key, value] of Object.entries(fm.data)) {
      const line = keyLine(fm, key) ?? 1;
      const error = (code: string, message: string, target?: string) => out.errors.push({file, line, code, message, target});
      const entry: Partial<MediaRecord> = {file, line};
      if (hosted) {
        const m = HOSTED_KEY.exec(key);
        if (!m) {
          error('unknown-media', `'${key}' is not a hosted video key (youtube:<id>)`, key);
          continue;
        }
        entry.provider = m[1];
        entry.externalId = m[2];
        entry.id = hostedMediaId(m[1], m[2]);
      } else {
        const target = path.posix.normalize(path.posix.join(dir, key));
        const exists = !key.startsWith('/') && !target.startsWith('..') && fs.existsSync(path.join(repoRoot, target)) && fs.statSync(path.join(repoRoot, target)).isFile();
        if (!exists || !isMedia(target)) {
          error('unknown-media', `'${key}' names no media file beside ${file}`, target);
          continue;
        }
        entry.path = target;
        entry.id = mediaId(target);
      }
      if (!value || typeof value !== 'object' || Array.isArray(value)) {
        error('bad-value', `'${key}': expected a map of record keys, got ${JSON.stringify(value)}`, key);
        continue;
      }
      const data = value as Record<string, unknown>;
      const allowed = hosted ? HOSTED_KEYS : RECORD_KEYS;
      let bad = false;
      for (const k of Object.keys(data))
        if (!allowed.includes(k)) {
          error('unknown-key', `'${key}': unknown record key '${k}' (allowed: ${allowed.join(', ')})`, key);
          bad = true;
        }
      const {illustrates, title, ...members} = data;
      const targets: string[] = [];
      if (illustrates !== undefined) {
        if (!Array.isArray(illustrates) || !illustrates.every((t) => typeof t === 'string' && t.trim())) {
          error('bad-edge', `'${key}': illustrates must be a list of feature or concept ids`, key);
          bad = true;
        } else {
          illustrates.forEach((t: string, i: number) => {
            const problem = resolve(t, ILLUSTRATES_TARGETS.filter((x) => system.nodes.has(x)), {source: file, key: `${key}.illustrates[${i}]`});
            if (problem) {
              error('unresolved-ref', `'${key}': illustrates[${i}]: ${problem}`, key);
              bad = true;
            } else targets.push(t.trim().replace(/^~/, ''));
          });
        }
      }
      const row: Row = {type: 'media', id: entry.id, format: entry.path ? formatOf(entry.path) : 'youtube', ...members};
      if (typeof title === 'string') row.name = title;
      const normalized = normalizeRow(system, row, {defaults: false});
      for (const p of normalized.problems)
        if (!bad || p.code !== 'unknown-key') {
          error(p.code, `'${key}': ${p.message}`, key);
          bad = true;
        }
      if (bad) continue;
      const {type, id, format, ...rest} = normalized.row;
      entry.data = rest;
      entry.illustrates = targets;
      if (entry.path && typeof rest.described_blob === 'string') {
        blobs ??= blobIds(repoRoot);
        const blob = blobs.get(entry.path) ?? hashBlob(fs.readFileSync(path.join(repoRoot, entry.path)));
        if (blob !== rest.described_blob)
          out.warnings.push({file, line, code: 'stale-description', message: `'${key}': described at blob ${String(rest.described_blob).slice(0, 12)}, the file is now ${blob.slice(0, 12)}; run grok kg enrich media --stale`, target: entry.path});
      }
      if (out.records.has(entry.id!)) error('duplicate-id', `'${key}': ${entry.id} already has a record`, key);
      else out.records.set(entry.id!, entry as MediaRecord);
    }
  }
  return out;
}

/** The git blob id of a buffer, as `git hash-object` computes it for a file no filter touches. */
export function hashBlob(content: Buffer): string {
  return createHash('sha1').update(`blob ${content.length}\0`).update(content).digest('hex');
}

/**
 * Blob ids of every tracked media file under [repoRoot], from the index of the monorepo and of the `public/`
 * submodule (and of any extra git root, keyed by the prefix its paths carry), without reading the files; a file
 * the working tree has changed is hashed the way git would. A tree that is no repository yields nothing.
 */
export function blobIds(repoRoot: string, roots: {dir: string, prefix: string}[] = [{dir: repoRoot, prefix: ''}, {dir: path.join(repoRoot, 'public'), prefix: 'public/'}]): Map<string, string> {
  const out = new Map<string, string>();
  const seen = new Set<string>();
  for (const {dir, prefix} of roots) {
    const top = spawnSync('git', ['rev-parse', '--show-toplevel'], {cwd: dir, encoding: 'utf8'});
    if (top.status !== 0) continue;
    const cwd = path.resolve(top.stdout.trim());
    if (seen.has(cwd) || path.resolve(dir) !== cwd) continue;
    seen.add(cwd);
    const ls = spawnSync('git', ['ls-files', '-s', '-z'], {cwd, encoding: 'utf8', maxBuffer: 256 * 1024 * 1024});
    if (ls.status !== 0) continue;
    for (const entry of ls.stdout.split('\0')) {
      // <mode> <blob> <stage>\t<path>
      const tab = entry.indexOf('\t');
      if (tab < 0) continue;
      const [, blob, stage] = entry.slice(0, tab).split(' ');
      const file = entry.slice(tab + 1);
      if (stage !== '0' || !isMedia(file)) continue;
      out.set(prefix + file, blob);
    }
    const status = spawnSync('git', ['status', '--porcelain', '-z', '--untracked-files=no'], {cwd, encoding: 'utf8', maxBuffer: 256 * 1024 * 1024});
    if (status.status !== 0) continue;
    const changed: string[] = [];
    const entries = status.stdout.split('\0');
    for (let i = 0; i < entries.length; i++) {
      const entry = entries[i];
      if (!entry) continue;
      const code = entry.slice(0, 2);
      const file = entry.slice(3);
      if (code.startsWith('R') || code.startsWith('C')) i++;
      if (code.includes('D')) out.delete(prefix + file);
      else if (isMedia(file) && fs.existsSync(path.join(cwd, file))) changed.push(file);
    }
    if (!changed.length) continue;
    const hashed = spawnSync('git', ['hash-object', '--', ...changed], {cwd, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024});
    if (hashed.status !== 0) continue;
    hashed.stdout.trim().split('\n').forEach((blob, i) => out.set(prefix + posix(changed[i]), blob.trim()));
  }
  return out;
}
