/// `grok kg enrich media` (conventions.md §11.1, the prepare-then-apply contract): for the media files the current
/// generation shows but no record describes (or describes for an older blob, `--stale`), sample frames with ffmpeg,
/// ask a cheaper model through `claude -p` what the file shows, validate the answer against the media type and merge
/// it into the folder's `media.yaml` as a proposal (`reviewed: false`, `described_by: <model>`, `described_blob`).
/// The build never runs this; a person reviews what it wrote.
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {createHash} from 'crypto';
import {spawnSync} from 'child_process';
import {globSync} from 'glob';
import * as yaml from 'js-yaml';
import {TypeSystem} from '../types';
import {Row, normalizeRow, compare} from '../normalize';
import {readJsonl, dataFile} from '../generation';
import {Roots, localPath} from '../roots';
import {proseLines} from '../citations';

/** The keys a proposal may set; `reviewed`, `described_by` and `described_blob` are the tool's, never the model's. */
export const PROPOSAL_KEYS = ['caption', 'description', 'actions', 'ui_text', 'kind', 'quality', 'quality_notes', 'illustrates'];
/** Bumped when the prompt or the answer shape changes, so cached answers to the old prompt are not reused. */
export const PROMPT_VERSION = 1;
const FRAMES = 6;
const FRAME_WIDTH = 1024;
const DEFAULT_MODEL = 'claude-sonnet-5';
/** The describer must be a cheaper model: a frontier id is refused whatever the flag says. */
const REFUSED_MODELS = /opus|fable/i;
const DESCRIBABLE = new Set(['png', 'jpg', 'gif', 'webp', 'mp4', 'webm']);
const CALL_TIMEOUT_MS = 300_000;
const JSON_BLOCK = /```(?:json)?\s*([\s\S]*?)```/;

export interface EnrichOptions {
  kgRoot: string;
  repoRoot: string;
  landingDir?: string;
  outRoot: string;
  system: TypeSystem;
  limit: number;
  only?: string;
  stale: boolean;
  dryRun: boolean;
  model: string;
  /** Test seams: the model call and the frame sampler. */
  describe?: Describer;
  frames?: FrameSampler;
}

export interface Prepared {
  id: string;
  path: string;
  format: string;
  blob: string;
  bytes: number;
  local: string;
  pages: Occurrence[];
  candidates: {id: string, name: string}[];
  frames: Frame[];
  probe?: Probe;
  /** ffmpeg absent: the frames are the file itself (a still) or its thumbnail (an animation), and the prompt says so. */
  stillOnly: boolean;
  prompt: string;
  cacheKey: string;
}

export interface Occurrence {
  page: string;
  title: string;
  anchor?: string;
  alt?: string;
  paragraph: string;
}

export interface Frame {
  file: string;
  seconds?: number;
}

export interface Probe {
  width?: number;
  height?: number;
  seconds?: number;
  frames?: number;
}

export interface Proposal {
  caption?: string;
  description?: string;
  actions?: string[];
  ui_text?: string[];
  kind?: string;
  quality?: string;
  quality_notes?: string;
  illustrates?: string[];
}

export type Describer = (item: Prepared, model: string) => {proposal?: Proposal, error?: string, raw?: string};
export type FrameSampler = (item: {local: string, format: string}, dir: string) => {frames: Frame[], probe?: Probe, stillOnly: boolean};

export interface Outcome {
  id: string;
  pages: number;
  frames: number;
  status: 'described' | 'cached' | 'failed' | 'dry-run' | 'skipped';
  caption?: string;
  detail?: string;
  record?: string;
}

export async function enrichMedia(o: EnrichOptions): Promise<{items: Prepared[], outcomes: Outcome[], notes: string[]}> {
  if (REFUSED_MODELS.test(o.model)) throw new Error(`--model ${o.model}: the describer is a cheaper model; a frontier model is refused`);
  const roots: Roots = {repoRoot: o.repoRoot, landingDir: o.landingDir};
  const genDir = o.outRoot;
  const rows = (kind: 'nodes' | 'edges', name: string): Row[] => [...readJsonl(dataFile(genDir, kind, name))] as Row[];
  const media = rows('nodes', 'media');
  const embeds = rows('edges', 'embeds');
  const pages = new Map(rows('nodes', 'doc-page').concat(rows('nodes', 'web-page')).map((p) => [String(p.id), p]));
  const features = new Map(rows('nodes', 'feature').map((f) => [String(f.id), f]));
  const documents = rows('edges', 'documents');
  const shownBy = new Map<string, Row[]>();
  for (const e of embeds) (shownBy.get(String(e.to)) ?? shownBy.set(String(e.to), []).get(String(e.to))!).push(e);
  const notes: string[] = [];
  const sampler = o.frames ?? ffmpegFrames;
  const work = media.filter((m) => typeof m.path === 'string' && DESCRIBABLE.has(String(m.format)) && m.reviewed !== true && shownBy.has(String(m.id)))
    .filter((m) => m.description === undefined || o.stale && typeof m.described_blob === 'string' && m.described_blob !== m.blob)
    .filter((m) => matches(String(m.path), o.only))
    .sort((a, b) => (shownBy.get(String(b.id))?.length ?? 0) - (shownBy.get(String(a.id))?.length ?? 0) || compare(String(a.id), String(b.id)))
    .slice(0, o.limit);
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-enrich-'));
  const items: Prepared[] = [];
  for (const m of work) {
    const local = localPath(roots, String(m.path));
    if (!local || !fs.existsSync(local)) continue;
    const occurrences = (shownBy.get(String(m.id)) ?? []).map((e) => occurrence(e, pages, roots)).filter((x): x is Occurrence => !!x);
    const candidates = candidateFeatures(occurrences, documents, features);
    const dir = path.join(scratch, createHash('sha1').update(String(m.id)).digest('hex').slice(0, 12));
    fs.mkdirSync(dir, {recursive: true});
    const sampled = o.dryRun ? {frames: [], stillOnly: false} : sampler({local, format: String(m.format)}, dir);
    const item: Prepared = {id: String(m.id), path: String(m.path), format: String(m.format), blob: String(m.blob), bytes: Number(m.bytes ?? 0), local,
      pages: occurrences, candidates, frames: sampled.frames, probe: sampled.probe, stillOnly: sampled.stillOnly, prompt: '', cacheKey: ''};
    item.prompt = prompt(item);
    item.cacheKey = createHash('sha1').update([item.blob, o.model, String(PROMPT_VERSION), `${FRAMES}x${FRAME_WIDTH}`, item.stillOnly ? 'still' : 'frames',
      occurrences.map((x) => `${x.page}#${x.anchor ?? ''}|${x.alt ?? ''}`).join(';'), candidates.map((c) => c.id).join(',')].join('\n')).digest('hex');
    items.push(item);
  }
  if (!items.length) notes.push('nothing to describe: every shown media file the filter admits has a description (use --stale to refresh stale ones)');
  if (o.dryRun) return {items, outcomes: items.map((i) => ({id: i.id, pages: i.pages.length, frames: 0, status: 'dry-run'})), notes};
  if (items.some((i) => i.stillOnly)) notes.push('ffmpeg not found: animations were judged by their first frame only, and never rated above docs');
  const describe = o.describe ?? claudeDescribe;
  const cacheDir = path.join(o.outRoot, '..', '..', 'enrich', 'media');
  fs.mkdirSync(cacheDir, {recursive: true});
  const outcomes: Outcome[] = [];
  for (const item of items) {
    const cached = path.join(cacheDir, `${item.cacheKey}.json`);
    let proposal: Proposal | undefined;
    let status: Outcome['status'] = 'described';
    let detail: string | undefined;
    if (fs.existsSync(cached)) {
      proposal = JSON.parse(fs.readFileSync(cached, 'utf8')).proposal;
      status = 'cached';
    } else {
      const answer = describe(item, o.model);
      if (!answer.proposal) {
        outcomes.push({id: item.id, pages: item.pages.length, frames: item.frames.length, status: 'failed', detail: answer.error ?? 'no answer'});
        continue;
      }
      proposal = answer.proposal;
    }
    const checked = validate(o.system, proposal, item);
    if (checked.error) {
      outcomes.push({id: item.id, pages: item.pages.length, frames: item.frames.length, status: 'failed', detail: checked.error});
      continue;
    }
    if (status === 'described') fs.writeFileSync(cached, `${JSON.stringify({key: item.cacheKey, id: item.id, model: o.model, at: new Date().toISOString(), proposal}, null, 2)}\n`);
    const applied = apply(roots, item, checked.proposal!, o.model);
    outcomes.push({id: item.id, pages: item.pages.length, frames: item.frames.length, status: applied.skipped ? 'skipped' : status, caption: checked.proposal!.caption,
      detail: applied.skipped, record: applied.record});
  }
  return {items, outcomes, notes};
}

/** A path prefix, or a glob with `*` and `**`. */
function matches(file: string, only?: string): boolean {
  if (!only) return true;
  if (!/[*?]/.test(only)) return file.startsWith(only.replace(/\/+$/, ''));
  const re = new RegExp(`^${only.replace(/[.+^${}()|[\]\\]/g, '\\$&').replace(/\*\*\/?/g, '\u0000').replace(/\*/g, '[^/]*').replace(/\u0000/g, '(?:.*/)?')}$`);
  return re.test(file);
}

/** Where the page shows the file: its title, the heading above, the alt text, and the paragraph around the embed. */
function occurrence(e: Row, pages: Map<string, Row>, roots: Roots): Occurrence | undefined {
  const page = pages.get(String(e.from));
  if (!page) return undefined;
  let paragraph = '';
  const local = localPath(roots, String(page.path));
  if (local && fs.existsSync(local) && /\.mdx?$/i.test(String(page.path))) {
    const lines = fs.readFileSync(local, 'utf8').split('\n');
    const at = Number(e.line ?? 0) - 1;
    const around = lines.slice(Math.max(0, at - 4), at + 5).filter((l) => l.trim() && !/^\s*(#|!\[|<|import\b|---)/.test(l));
    paragraph = around.join(' ').replace(/\s+/g, ' ').slice(0, 600);
  }
  return {page: String(page.path), title: String(page.name ?? page.id), anchor: e.anchor === undefined ? undefined : String(e.anchor), alt: e.alt === undefined ? undefined : String(e.alt), paragraph};
}

/** The features a proposal may name: those the showing pages document or are home of, and the parents of those. */
function candidateFeatures(occurrences: Occurrence[], documents: Row[], features: Map<string, Row>): {id: string, name: string}[] {
  const ids = new Set<string>();
  const pageIds = new Set(occurrences.map((x) => `doc:${x.page}`));
  for (const d of documents) if (pageIds.has(String(d.from))) ids.add(String(d.to));
  for (const f of features.values())
    if (typeof f.home === 'string' && occurrences.some((x) => x.page === f.home) || pageIds.has(String(f.user_help)) || pageIds.has(String(f.developer_help))) ids.add(String(f.id));
  for (const id of [...ids]) {
    const parts = id.split('/');
    for (let i = parts.length - 1; i >= 2; i--) ids.add(parts.slice(0, i).join('/'));
  }
  return [...ids].sort(compare).filter((id) => features.has(id)).map((id) => ({id, name: String(features.get(id)!.name ?? id)}));
}

export function prompt(item: Prepared): string {
  const probe = item.probe ? ` ${[item.probe.width && item.probe.height ? `${item.probe.width}x${item.probe.height}` : '', item.probe.seconds ? `${item.probe.seconds.toFixed(1)} s` : '', item.probe.frames ? `${item.probe.frames} frames` : ''].filter(Boolean).join(', ')}` : '';
  const frames = item.frames.length ? item.frames.map((f) => `- ${f.file}${f.seconds === undefined ? '' : ` (t=${f.seconds.toFixed(1)} s)`}`).join('\n') : `- ${item.local}`;
  const shown = item.pages.map((x) => `- "${x.title}" (${x.page})${x.anchor ? ` under the heading "${x.anchor}"` : ''}${x.alt ? `, alt text "${x.alt}"` : ''}${x.paragraph ? `; the text around it: "${x.paragraph}"` : ''}`).join('\n');
  const candidates = item.candidates.length ? item.candidates.map((c) => `- ${c.id} (${c.name})`).join('\n') : '- (none: leave illustrates empty)';
  return [
    'You describe one media file of the Datagrok documentation for a knowledge graph that marketing and support search. Read every file listed under "Frames" with the Read tool, then answer with exactly one fenced ```json block and nothing else.',
    '',
    `Asset: ${item.path} (${item.format}, ${item.bytes} bytes${probe})${item.stillOnly && /gif|mp4|webm/.test(item.format) ? ' — FIRST FRAME ONLY: describe what is visible, do not infer actions, do not rate above docs' : ''}`,
    '',
    'Shown on:',
    shown || '- (no page)',
    '',
    'Feature ids that may appear in illustrates (only these; the list may be empty):',
    candidates,
    '',
    'Frames:',
    frames,
    '',
    'Answer with this shape:',
    '{"caption": "at most 12 words, what it shows", "description": "2-4 sentences a reader would understand without the picture: what is on screen and, for an animation, what happens in order", "actions": ["for an animation: the user actions shown, in order; [] for a still"], "ui_text": ["labels, menu items, dialog titles legible in the picture"], "kind": "screenshot | animation | clip | diagram | icon | logo | thumbnail | badge | photo", "quality": "unfit | docs | answer | marketing", "quality_notes": "why not higher: cropping, clutter, low resolution, debug data, an old UI; empty for marketing", "illustrates": ["feature ids from the list above that the picture actually shows"]}',
    '',
    'Quality tiers: unfit = do not show anyone (placeholder, broken, unreadable, obviously outdated UI, test data); docs = fine inside its page, not on its own; answer = clear enough to send alone to a user asking how to do this; marketing = polished, fit for a presentation or the website.',
  ].join('\n');
}

/** The model call: `claude -p` reads the prompt on stdin, may only Read files, answers as JSON; the answer is the fenced block in its result. */
export const claudeDescribe: Describer = (item, model) => {
  const r = spawnSync('claude', ['-p', '--model', model, '--output-format', 'json', '--allowedTools', 'Read', '--max-turns', '12'],
    {input: item.prompt, encoding: 'utf8', shell: process.platform === 'win32', timeout: CALL_TIMEOUT_MS, maxBuffer: 16 * 1024 * 1024});
  if (r.error) return {error: `claude: ${r.error.message}`};
  if (r.status !== 0) return {error: `claude exited ${r.status}: ${(r.stderr || r.stdout || '').trim().slice(0, 300)}`, raw: r.stdout};
  let text = r.stdout;
  try {
    const envelope = JSON.parse(r.stdout);
    if (envelope.is_error) return {error: `claude: ${String(envelope.result ?? '').slice(0, 300)}`, raw: r.stdout};
    text = typeof envelope.result === 'string' ? envelope.result : JSON.stringify(envelope.structured_output ?? envelope.result ?? '');
  } catch {
    // not an envelope: the raw text may still carry the block
  }
  return parseAnswer(text, r.stdout);
};

export function parseAnswer(text: string, raw?: string): {proposal?: Proposal, error?: string, raw?: string} {
  const block = JSON_BLOCK.exec(text)?.[1] ?? (text.trim().startsWith('{') ? text.trim() : undefined);
  if (!block) return {error: 'no JSON block in the answer', raw};
  try {
    const parsed = JSON.parse(block);
    if (!parsed || typeof parsed !== 'object' || Array.isArray(parsed)) return {error: 'the answer is not an object', raw};
    return {proposal: parsed as Proposal, raw};
  } catch (e: any) {
    return {error: `malformed JSON in the answer: ${e.message}`, raw};
  }
}

/** Only the proposal keys, typed as the media members are, illustrates within the candidates, quality capped without real frames. */
export function validate(system: TypeSystem, proposal: Proposal, item: Prepared): {proposal?: Proposal, error?: string} {
  const out: Record<string, unknown> = {};
  for (const [k, v] of Object.entries(proposal)) {
    if (!PROPOSAL_KEYS.includes(k) || v === null || v === undefined || v === '') continue;
    out[k] = v;
  }
  const {illustrates, ...members} = out;
  const normalized = normalizeRow(system, {type: 'media', id: item.id, format: item.format, ...members}, {defaults: false});
  if (normalized.problems.length) return {error: normalized.problems.map((p) => p.message).join('; ')};
  const {type, id, format, ...rest} = normalized.row;
  const allowed = new Set(item.candidates.map((c) => c.id));
  const shows = Array.isArray(illustrates) ? illustrates.filter((x): x is string => typeof x === 'string' && allowed.has(x)) : [];
  if (Array.isArray(rest.actions) && !rest.actions.length) delete rest.actions;
  if (Array.isArray(rest.ui_text) && !rest.ui_text.length) delete rest.ui_text;
  if (item.stillOnly && /gif|mp4|webm/.test(item.format) && (rest.quality === 'answer' || rest.quality === 'marketing')) rest.quality = 'docs';
  if (typeof rest.description !== 'string') return {error: 'no description in the answer'};
  return {proposal: {...rest, ...(shows.length ? {illustrates: shows} : {})} as Proposal};
}

/** Merges the proposal into the folder's media.yaml: a reviewed entry is left alone, every other entry is replaced by the proposal
 * with the tool's own fields; keys sorted, the blob quoted so YAML never reads it as a number. */
export function apply(roots: Roots, item: Prepared, proposal: Proposal, model: string): {record: string, skipped?: string} {
  const dir = path.posix.dirname(item.path);
  const record = `${dir}/media.yaml`;
  const local = localPath(roots, record)!;
  const existing = fs.existsSync(local) ? (yaml.load(fs.readFileSync(local, 'utf8')) as Record<string, unknown> | null) ?? {} : {};
  const key = path.posix.basename(item.path);
  const current = existing[key];
  if (current && typeof current === 'object' && (current as Record<string, unknown>).reviewed === true) return {record, skipped: 'reviewed record kept'};
  const entry: Record<string, unknown> = {...proposal, reviewed: false, described_by: model, described_blob: item.blob};
  if (item.probe?.width) entry.width = item.probe.width;
  if (item.probe?.height) entry.height = item.probe.height;
  if (item.probe?.seconds) entry.seconds = Math.round(item.probe.seconds * 10) / 10;
  existing[key] = entry;
  const sorted = Object.fromEntries(Object.entries(existing).sort(([a], [b]) => compare(a, b)));
  fs.writeFileSync(local, yaml.dump(sorted, {lineWidth: 100, sortKeys: true, noRefs: true}));
  return {record};
}

/** Six frames spread over the clip (ffmpeg), or the still itself; without ffmpeg an animation falls back to its thumbnail or first frame. */
export const ffmpegFrames: FrameSampler = ({local, format}, dir) => {
  const animated = /gif|mp4|webm/.test(format);
  const ffmpeg = tool('ffmpeg');
  const ffprobe = tool('ffprobe');
  let probe: Probe | undefined;
  if (ffprobe) {
    const r = spawnSync(ffprobe, ['-v', 'error', '-select_streams', 'v:0', '-show_entries', 'stream=width,height,nb_frames,duration,r_frame_rate:format=duration', '-of', 'json', local], {encoding: 'utf8'});
    if (r.status === 0) {
      try {
        const info = JSON.parse(r.stdout);
        const s = info.streams?.[0] ?? {};
        const seconds = Number(s.duration ?? info.format?.duration);
        probe = {width: s.width, height: s.height, frames: s.nb_frames ? Number(s.nb_frames) : undefined, seconds: Number.isFinite(seconds) && seconds > 0 ? seconds : undefined};
      } catch {
        probe = undefined;
      }
    }
  }
  if (!animated) return {frames: [{file: local}], probe, stillOnly: false};
  if (!ffmpeg) {
    const thumb = local.replace(/\.gif$/i, '-thumb.png');
    return {frames: [{file: fs.existsSync(thumb) ? thumb : local, seconds: 0}], probe, stillOnly: true};
  }
  const seconds = probe?.seconds;
  const scale = `scale=min(${FRAME_WIDTH}\\,iw):-2`;
  const filter = seconds && seconds > 1 ? `fps=${FRAMES / seconds},${scale}` : scale;
  const r = spawnSync(ffmpeg, ['-loglevel', 'error', '-y', '-i', local, '-vf', filter, '-frames:v', String(FRAMES), '-fps_mode', 'vfr', path.join(dir, 'f%02d.png')], {encoding: 'utf8'});
  const files = fs.existsSync(dir) ? fs.readdirSync(dir).filter((f) => /^f\d+\.png$/.test(f)).sort() : [];
  if (r.status !== 0 || !files.length) return {frames: [{file: local, seconds: 0}], probe, stillOnly: true};
  const step = seconds && seconds > 1 ? seconds / FRAMES : undefined;
  return {frames: files.map((f, i) => ({file: path.join(dir, f), seconds: step === undefined ? undefined : Math.round(i * step * 10) / 10})), probe, stillOnly: false};
};

/** ffmpeg or ffprobe on the PATH, else where winget puts Gyan.FFmpeg. */
function tool(name: string): string | undefined {
  if (spawnSync(name, ['-version'], {encoding: 'utf8'}).status === 0) return name;
  const home = process.env.LOCALAPPDATA;
  if (!home) return undefined;
  const found = globSync('Microsoft/WinGet/Packages/Gyan.FFmpeg*/**/bin/' + name + '.exe', {cwd: home, absolute: true, windowsPathsNoEscape: true});
  return found.sort()[found.length - 1];
}
