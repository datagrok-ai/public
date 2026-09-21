/* Guide mode: with `BDD_GUIDE=<dir>` set, a run leaves a manifest per scenario — for every step
   the page before and after it, the element the step located, where the pointer went and what
   it clicked — that `tool/guide-render.py` turns into a how-to video. Without the variable
   nothing here runs; a step of a test costs nothing. */
import {mkdirSync, writeFileSync} from 'node:fs';
import {join} from 'node:path';
import type {Locator, Page, TestInfo} from '@playwright/test';

export interface GuideBox {
  x: number;
  y: number;
  width: number;
  height: number;
}

export interface GuidePoint {
  x: number;
  y: number;
}

export type GuideStepKind = 'action' | 'check' | 'setup';

export interface GuideStep {
  index: number;
  line: number;
  keyword: string;
  text: string;
  /** The step as an instruction to a person: "Click the Open local file icon in the Browse toolbar". */
  caption: string;
  kind: GuideStepKind;
  before: string;
  after: string;
  /** The element the step acted on or checked, in page pixels. */
  target?: GuideBox;
  /** Where the pointer went during the step, in order. */
  pointer: GuidePoint[];
  clicks: {x: number; y: number; button: string}[];
  /** Keys pressed during the step (`Enter`, `Control+A`), and text typed. */
  keys: string[];
  typed: string;
  ms: number;
}

export interface GuideManifest {
  feature: string;
  scenario: string;
  description: string;
  tags: string[];
  viewport: {width: number; height: number} | null;
  steps: GuideStep[];
}

interface Recording {
  dir: string;
  manifest: GuideManifest;
  lastType: string;
  pointer: GuidePoint[];
  clicks: {x: number; y: number; button: string}[];
  keys: string[];
  typed: string;
  located?: Locator;
  target?: GuideBox | null;
  open?: {index: number; line: number; keyword: string; text: string; started: number; before: string};
}

const recordings = new WeakMap<Page, Recording>();
const attached = new WeakSet<Page>();

export function guideDir(): string | undefined {
  return process.env.BDD_GUIDE || undefined;
}

/** The pause before the "after" screenshot: what the platform animates (a dialog sliding in, a
 * balloon) has to be on the page for the picture, while a test never waits for it. */
function settleMs(): number {
  return Number(process.env.BDD_GUIDE_SETTLE ?? 500);
}

export function slugOf(name: string): string {
  return name.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '').slice(0, 60) || 'untitled';
}

/** `When user clicks on browse tab` → keyword `When`, text `user clicks on browse tab`. */
export function parseTitle(title: string): {keyword: string; text: string} {
  const m = /^(Given|When|Then|And|But|\*)\s+(.*)$/s.exec(title.trim());
  return m ? {keyword: m[1], text: m[2]} : {keyword: '', text: title.trim()};
}

const VERBS: Record<string, string> = {
  clicks: 'Click', 'double-clicks': 'Double-click', 'right-clicks': 'Right-click', hovers: 'Hover', types: 'Type',
  presses: 'Press', uploads: 'Upload', opens: 'Open', closes: 'Close', chooses: 'Choose', selects: 'Select',
  picks: 'Pick', drags: 'Drag', checks: 'Check', unchecks: 'Uncheck', toggles: 'Toggle', enters: 'Enter',
  expands: 'Expand',
  collapses: 'Collapse', saves: 'Save', adds: 'Add', removes: 'Remove', sets: 'Set', switches: 'Switch',
  scrolls: 'Scroll', moves: 'Move', waits: 'Wait', resizes: 'Resize', focuses: 'Focus', pastes: 'Paste',
  clears: 'Clear', deletes: 'Delete', creates: 'Create', renames: 'Rename', shares: 'Share', runs: 'Run',
  submits: 'Submit', confirms: 'Confirm', cancels: 'Cancel', searches: 'Search', filters: 'Filter', sorts: 'Sort',
  zooms: 'Zoom', pans: 'Pan', reads: 'Read', navigates: 'Go', goes: 'Go', edits: 'Edit', applies: 'Apply',
};

function imperative(verb: string): string {
  const known = VERBS[verb.toLowerCase()];
  if (known)
    return known;
  let stem = verb;
  if (/ies$/.test(verb))
    stem = verb.slice(0, -3) + 'y';
  else if (/(ss|sh|ch|x|z)es$/.test(verb))
    stem = verb.slice(0, -2);
  else if (/s$/.test(verb) && !/ss$/.test(verb))
    stem = verb.slice(0, -1);
  return stem.charAt(0).toUpperCase() + stem.slice(1);
}

/** Element phrases as a reader sees them: `inside` is `in`, a tree path `Files---Demo` is
 * `Files › Demo`, `on browse tab` keeps its words. */
function readable(text: string): string {
  return text.replace(/\binside\b/g, 'in').replace(/---/g, ' › ').replace(/\s+/g, ' ').trim();
}

/** A step as the instruction it gives, or as the fact it checks. `user clicks on browse tab` →
 * `Click on browse tab`; `the table should have 5 rows` → `The table has 5 rows`. */
export function captionOf(text: string, kind: GuideStepKind): string {
  let t = readable(text).replace(/^(the )?user\s+/i, '');
  if (kind === 'check') {
    t = t.replace(/\bshould not be\b/g, 'is not').replace(/\bshould not have\b/g, 'does not have')
      .replace(/\bshould be\b/g, 'is').replace(/\bshould have\b/g, 'has')
      .replace(/\bshould (\w+)/g, (_, v: string) => v + 's').replace(/\bshould\b/g, '');
    return t.charAt(0).toUpperCase() + t.slice(1);
  }
  const m = /^(\S+)(.*)$/s.exec(t);
  if (!m)
    return t;
  // a state ("user is logged in", "the panel has a tab") is told as it is, not ordered
  if (/^(is|are|has|have|was|were)$/.test(m[1])) {
    const state = readable(text);
    return state.charAt(0).toUpperCase() + state.slice(1);
  }
  return imperative(m[1]) + m[2];
}

function typeOf(keyword: string, lastType: string): string {
  return keyword === 'And' || keyword === 'But' || keyword === '*' ? lastType : keyword;
}

/** The pointer is followed on the page's own `mouse`, whatever runtime path drives it (a gesture,
 * a menu walk, a viewer hit area); the methods are replaced in place, once per page. */
export function attach(page: Page): void {
  if (!guideDir() || attached.has(page))
    return;
  attached.add(page);
  const mouse = page.mouse as any;
  const move = mouse.move.bind(mouse);
  const click = mouse.click.bind(mouse);
  const dblclick = mouse.dblclick.bind(mouse);
  const down = mouse.down.bind(mouse);
  let at: GuidePoint = {x: 0, y: 0};
  const rec = (): Recording | undefined => recordings.get(page);
  mouse.move = async (x: number, y: number, options?: unknown) => {
    at = {x, y};
    rec()?.pointer.push(at);
    return move(x, y, options);
  };
  mouse.click = async (x: number, y: number, options?: {button?: string}) => {
    at = {x, y};
    rec()?.pointer.push(at);
    rec()?.clicks.push({x, y, button: options?.button ?? 'left'});
    return click(x, y, options);
  };
  mouse.dblclick = async (x: number, y: number, options?: {button?: string}) => {
    at = {x, y};
    rec()?.pointer.push(at);
    rec()?.clicks.push({x, y, button: options?.button ?? 'left'});
    return dblclick(x, y, options);
  };
  mouse.down = async (options?: {button?: string}) => {
    rec()?.clicks.push({x: at.x, y: at.y, button: options?.button ?? 'left'});
    return down(options);
  };
  const keyboard = page.keyboard as any;
  const press = keyboard.press.bind(keyboard);
  const type = keyboard.type.bind(keyboard);
  const insertText = keyboard.insertText.bind(keyboard);
  keyboard.press = async (key: string, options?: unknown) => {
    rec()?.keys.push(key);
    return press(key, options);
  };
  keyboard.type = async (text: string, options?: unknown) => {
    const r = rec();
    if (r)
      r.typed += text;
    return type(text, options);
  };
  keyboard.insertText = async (text: string) => {
    const r = rec();
    if (r)
      r.typed += text;
    return insertText(text);
  };
}

/** What the step found on the page: the last phrase resolved wins (a scope resolves before the
 * element inside it). Its rectangle is read at once when it is already there, else at the step's
 * end (a menu item that the click closed is read by then, an element that appears later is not). */
export async function located(page: Page, loc: Locator): Promise<void> {
  const r = recordings.get(page);
  if (!r || !r.open)
    return;
  r.located = loc;
  r.target = await loc.first().boundingBox({timeout: 200}).catch(() => null);
}

function recordingFor(page: Page, info: TestInfo): Recording {
  let r = recordings.get(page);
  const [, featureName = 'feature', scenarioName = 'scenario'] = info.titlePath;
  const dir = join(guideDir()!, slugOf(featureName), slugOf(scenarioName));
  if (r && r.dir === dir)
    return r;
  mkdirSync(dir, {recursive: true});
  const viewport = page.viewportSize();
  const manifest: GuideManifest =
    {feature: featureName, scenario: scenarioName, description: '', tags: info.tags, viewport, steps: []};
  r = {dir, manifest, lastType: 'Given', pointer: [], clicks: [], keys: [], typed: ''};
  recordings.set(page, r);
  return r;
}

async function shot(page: Page, dir: string, name: string): Promise<string> {
  await page.screenshot({path: join(dir, name), animations: 'disabled', caret: 'hide'}).catch(() => undefined);
  return name;
}

/** Opens a step: the page as it is before the step, for the picture the pointer moves over. A step
 * before the page exists (the login) is recorded without pictures. */
export async function begin(page: Page | undefined, info: TestInfo, line: number, title: string): Promise<void> {
  if (!guideDir() || !page || page.isClosed())
    return;
  const r = recordingFor(page, info);
  const {keyword, text} = parseTitle(title);
  const index = r.manifest.steps.length + 1;
  const stem = String(index).padStart(2, '0');
  r.pointer = [];
  r.clicks = [];
  r.keys = [];
  r.typed = '';
  r.located = undefined;
  r.target = undefined;
  r.open = {index, line, keyword, text, started: Date.now(), before: await shot(page, r.dir, `${stem}-before.png`)};
}

/** Closes the step: waits for the platform to settle, takes the "after" picture, resolves the
 * target rectangle and writes the manifest — after every step, so a run that stops halfway
 * leaves what it had. */
export async function end(page: Page | undefined): Promise<void> {
  if (!guideDir() || !page || page.isClosed())
    return;
  const r = recordings.get(page);
  if (!r || !r.open)
    return;
  const open = r.open;
  r.open = undefined;
  await page.waitForTimeout(settleMs());
  await page.evaluate(() => new Promise<void>((resolve) =>
    requestAnimationFrame(() => requestAnimationFrame(() => resolve())))).catch(() => undefined);
  const stem = String(open.index).padStart(2, '0');
  const after = await shot(page, r.dir, `${stem}-after.png`);
  let target = r.target ?? undefined;
  if (!target && r.located)
    target = await r.located.first().boundingBox({timeout: 200}).catch(() => null) ?? undefined;
  const type = typeOf(open.keyword, r.lastType);
  r.lastType = type;
  const acted = !!target || r.pointer.length > 0 || r.clicks.length > 0 || r.keys.length > 0 || r.typed.length > 0;
  const kind: GuideStepKind = type === 'Then' ? 'check' : acted ? 'action' : 'setup';
  if (!target && r.clicks.length > 0) {
    const last = r.clicks[r.clicks.length - 1];
    target = {x: last.x - 12, y: last.y - 12, width: 24, height: 24};
  }
  r.manifest.steps.push({index: open.index, line: open.line, keyword: open.keyword, text: open.text,
    caption: captionOf(open.text, kind), kind, before: open.before, after, target, pointer: r.pointer, clicks: r.clicks,
    keys: r.keys, typed: r.typed, ms: Date.now() - open.started});
  writeFileSync(join(r.dir, 'steps.json'), JSON.stringify(r.manifest, null, 2));
}
