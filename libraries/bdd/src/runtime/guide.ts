/* Guide mode: with `BDD_GUIDE=<dir>` set, a run leaves a manifest per scenario — for every step
   the page before and after it, the element the step located, where the pointer went and what
   it clicked — that `tool/guide-render.py` turns into a how-to video. Without the variable
   nothing here runs; a step of a test costs nothing. */
import {mkdirSync, readFileSync, writeFileSync} from 'node:fs';
import {join} from 'node:path';
import type {Locator, Page, TestInfo} from '@playwright/test';

export interface GuideBox {
  x: number;
  y: number;
  width: number;
  height: number;
}

/** What the pointer did, in the page's CSS pixels — the screenshots' pixels. */
export interface GuidePointer {
  /** `move`: where the pointer came to rest, or with `held` a point of a drag; `down` and `up`: a
   * button pressed or released there. */
  type: 'move' | 'down' | 'up';
  x: number;
  y: number;
  button?: 'left' | 'middle' | 'right';
  /** 2 on the second press of a double-click. */
  count?: number;
  held?: boolean;
  /** The element a press landed on: its platform name (or tag) and its rectangle. */
  on?: string;
  box?: GuideBox;
  /** The page as it was with the button held here — what a drag has drawn so far. */
  shot?: string;
}

export type GuideStepKind = 'action' | 'check' | 'setup';

const DRAG_SHOTS = 8;

/** A stop of a step: the page as the pointer set off for it, the element lit there, and what the
 * pointer did at it. A step has one stop, or one per stop of a path it walked (a menu's group, each
 * item, the leaf). */
export interface GuideLeg {
  shot: string;
  target?: GuideBox;
  pointer: GuidePointer[];
}

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
  /** The element the step acted on or checked, in page pixels — the last stop when it walked a path. */
  target?: GuideBox;
  legs: GuideLeg[];
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
  keys: string[];
  typed: string;
  located?: Locator;
  target?: GuideBox | null;
  legs: GuideLeg[];
  /** What the pointer did since the last stop. */
  pending: GuidePointer[];
  dragShots: number;
  silent: boolean;
  open?: {index: number; line: number; keyword: string; text: string; table?: string[][]; started: number; before: string};
}

const recordings = new WeakMap<Page, Recording>();
const attached = new WeakSet<Page>();

export function guideDir(): string | undefined {
  return process.env.BDD_GUIDE || undefined;
}

/** A guide films the full shell — the menus, view tabs and panels as a person has them; a test
 * page runs in simple mode. */
export function shellSimpleMode(): boolean {
  return !guideDir();
}

/** The pause before the "after" screenshot: what the platform animates (a dialog sliding in, a
 * balloon) has to be on the page for the picture, while a test never waits for it. */
function settleMs(): number {
  return Number(process.env.BDD_GUIDE_SETTLE ?? 500);
}

/** Checks a viewer of the video has no use for — what a test needs to know, not what a person
 * sees: error and balloon floors, server state, the readings and pixels a viewer reports, claims
 * against a snapshot ("than before"), property bags, widget counts, task-bar and command
 * bookkeeping. A check that names what is on the page (a dialog, a column, a row count, a value,
 * a legend item's color) stays. Matched against the lowercased text of a `Then` (an `And` after
 * one included) — never an action: "drags across the "view" area of …" is a step to show. */
const HIDDEN_CHECKS: RegExp[] = [
  /^no errors should have been logged$/, /^no error or warning balloon should have been shown$/,
  / on the server$/,
  /^the top menu command should have completed$/, /^the package autostarts have completed$/, /\btask bar\b/,
  /\breadings? of\b/, /\bas remembered\b/, /^user remembers /, /\blistens for\b/, /\bshould have fired\b/,
  /\brepainted\b/, /\bpainted\b/, /\bink than before\b/, /\bcontain the color\b/, /\bthan before$/, /\bthan the ".*" area$/,
  /\bpixels tall$/, /\bareas? of\b/, /\bshould (not )?have an? ".*" area$/, /\bproperty of\b/, /^properties of /,
  /\bvalue range of\b/, /\bcolor scale of\b/, /\bcells of\b.*\bwide\b/, /\bshould show (fewer|more) rows\b/,
  /^the (open tableview|current view) should (have|hold)/, /\bshould be added to the open tableview$/,
  /^the .* view should be current$/, /\bnever increase$/, /\bshould have finished updating$/, / should match "/, /\bmatching "/,
  /^the legend of .* should (be (wider|narrower|taller|shorter|placed|docked|in a corner|in the|collapsed|on the)|list (fewer|the same))/,
  /^every item in the legend of/,
];

export function hiddenInGuide(text: string): boolean {
  const t = text.trim().toLowerCase();
  return HIDDEN_CHECKS.some((re) => re.test(t));
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
 * `Files › Demo`, `on browse tab` keeps its words, and the lines of a pasted text (`\n`) are
 * listed with commas, the way a filter card's search box shows them. */
function readable(text: string): string {
  return text.replace(/\binside\b/g, 'in').replace(/---/g, ' › ').replace(/(\\n)+"/g, '"').replace(/\\n/g, ', ')
    .replace(/\s+/g, ' ').trim();
}

/** `pass` → `passes`, `show` → `shows`, `fly` → `flies`. */
function thirdPerson(verb: string): string {
  if (/(s|sh|ch|x|z)$/.test(verb))
    return verb + 'es';
  if (/[^aeiou]y$/.test(verb))
    return verb.slice(0, -1) + 'ies';
  return verb + 's';
}

/** A step as the instruction it gives, or as the fact it checks. `user clicks on browse tab` →
 * `Click on browse tab`; `the table should have 5 rows` → `The table has 5 rows`; `5 rows should
 * pass the filter` → `5 rows pass the filter` (a counted or quantified plural in the subject, or a
 * plural word right before `should`, decides the number). */
export function captionOf(text: string, kind: GuideStepKind, table?: string[][]): string {
  let t = readable(text).replace(/^(the )?user\s+/i, '');
  // a data table is told inline: "adds a scatter plot viewer with:" + | X | AGE | → "… with X = AGE"
  if (table?.length)
    t = t.replace(/:\s*$/, '') + ' ' + table.map((row) => row.length === 2 ? `${row[0]} = ${row[1]}` : row.join(' | ')).join(', ');
  if (kind === 'check') {
    const subject = t.split(/\bshould\b/)[0].trim();
    const last = subject.split(/\s+/).pop() ?? '';
    const plural = (/s$/.test(last) && !/(ss|us|is|'s)$/.test(last)) ||
      /^(only |all |no |some )?(the )?(rows|columns)\b/.test(subject) || /\b(\d+|all|no|some|following)\s+[a-z]+s\b/.test(subject);
    t = t.replace(/\bshould not be\b/g, plural ? 'are not' : 'is not')
      .replace(/\bshould not have\b/g, plural ? 'do not have' : 'does not have')
      .replace(/\bshould not (\w+)/g, (_, v: string) => (plural ? 'do not ' : 'does not ') + v)
      .replace(/\bshould be\b/g, plural ? 'are' : 'is').replace(/\bshould have\b/g, plural ? 'have' : 'has')
      .replace(/\bshould (\w+)/g, (_, v: string) => plural ? v : thirdPerson(v)).replace(/\bshould\b/g, '');
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

/** Runs in the page, on every document it loads: logs each real press, release and move of the
 * pointer, whatever sent it — the page's mouse, a locator's click, a drag — where the browser
 * delivered it, and what the press landed on. A move with no button held only says where the
 * pointer came to rest, so a run of them keeps the last; a synthetic `el.click()` has no place on
 * the page and is not a gesture. */
function capturePointer(): void {
  const w = window as any;
  if (w.__bddGuidePointer)
    return;
  const log: any[] = w.__bddGuidePointer = [];
  const buttons = ['left', 'middle', 'right'];
  // the listeners are inline arrows: tsx wraps a named function in a `__name` call the page lacks
  for (const type of ['down', 'up', 'move'] as const) {
    addEventListener(`pointer${type}`, (event) => {
      const e = event as PointerEvent;
      if (!e.isTrusted)
        return;
      const held = e.buttons !== 0;
      const prev = log[log.length - 1];
      if (type === 'move' && !held && prev?.type === 'move' && !prev.held)
        log.pop();
      const entry: any = {type, x: e.clientX, y: e.clientY};
      if (type === 'move' && held)
        entry.held = true;
      if (type !== 'move')
        entry.button = buttons[e.button] ?? 'left';
      if (type === 'down') {
        entry.count = 1;
        const el = e.target instanceof Element ? e.target : null;
        if (el) {
          const r = el.getBoundingClientRect();
          entry.box = {x: r.x, y: r.y, width: r.width, height: r.height};
          entry.on = (el.closest('[name]')?.getAttribute('name') ?? el.tagName.toLowerCase()).slice(0, 80);
        }
      }
      log.push(entry);
      if (log.length > 2000)
        log.splice(0, log.length - 2000);
    }, {capture: true, passive: true});
  }
  // the browser says which press was a double-click's second: two quick clicks of their own (a
  // property cell, then the select it opened a few pixels away) are not one
  addEventListener('dblclick', (event) => {
    const e = event as MouseEvent;
    const press = [...log].reverse().find((p) => p.type === 'down');
    if (e.isTrusted && press && Math.abs(press.x - e.clientX) < 1 && Math.abs(press.y - e.clientY) < 1)
      press.count = 2;
  }, {capture: true, passive: true});
}

/** What the pointer did since the last read. */
async function drain(page: Page): Promise<GuidePointer[]> {
  return page.evaluate(() => {
    const log = (window as any).__bddGuidePointer as unknown[] | undefined;
    return (log ? log.splice(0, log.length) : []) as GuidePointer[];
  }).catch(() => []);
}

/** The pointer is followed in the page, whatever runtime path drives it (a gesture, a menu walk, a
 * viewer hit area, a locator's own click). The page's `mouse` is wrapped only to picture drags. */
export async function attach(page: Page): Promise<void> {
  if (!guideDir() || attached.has(page))
    return;
  attached.add(page);
  await page.addInitScript(capturePointer);
  await page.evaluate(capturePointer).catch(() => undefined);
  const rec = (): Recording | undefined => recordings.get(page);
  const mouse = page.mouse as any;
  const move = mouse.move.bind(mouse);
  const down = mouse.down.bind(mouse);
  const up = mouse.up.bind(mouse);
  let held = false;
  // a move with a button held is a drag: the page is pictured along the way (at most DRAG_SHOTS
  // per step), so the video shows what the drag draws — a selection box, an annotation region
  mouse.move = async (x: number, y: number, options?: unknown) => {
    const r = rec();
    if (!held || !r?.open || r.dragShots >= DRAG_SHOTS)
      return move(x, y, options);
    r.pending.push(...await drain(page));
    const last = [...r.pending].reverse().find((p) => p.type !== 'up');
    const from = last ?? {x, y};
    const legs = Math.min(3, Math.max(1, Math.round(Math.hypot(x - from.x, y - from.y) / 120)));
    for (let i = 1; i <= legs; i++) {
      const at = {x: from.x + (x - from.x) * i / legs, y: from.y + (y - from.y) * i / legs};
      await move(at.x, at.y, {steps: 4});
      r.pending.push(...await drain(page));
      r.dragShots++;
      const name = await shot(page, r.dir, `${String(r.open.index).padStart(2, '0')}-drag${r.dragShots}.png`);
      const reached = r.pending[r.pending.length - 1];
      if (reached?.type === 'move' && reached.held && Math.hypot(reached.x - at.x, reached.y - at.y) < 1)
        reached.shot = name;
      else
        r.pending.push({type: 'move', x: at.x, y: at.y, held: true, shot: name});
    }
  };
  mouse.down = async (options?: unknown) => {
    held = true;
    return down(options);
  };
  mouse.up = async (options?: unknown) => {
    held = false;
    return up(options);
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

/** A stop on the path the step walks: the page as it is now (a menu opened by the stop before)
 * and the element's place on it; the guide moves the pointer stop by stop and lights each. What
 * the pointer did before the first stop — on an element located before it, or anywhere (the
 * right-click that opened the menu) — is a stop of its own, on the "before" picture. */
export async function hop(page: Page, loc: Locator): Promise<void> {
  if (!recordings.get(page)?.open)
    return;
  const box = await loc.first().boundingBox({timeout: 200}).catch(() => null);
  if (box)
    await hopAt(page, box);
}

/** The same stop for a place that is no element of its own — a row of a canvas grid. */
export async function hopAt(page: Page, box: GuideBox): Promise<void> {
  const r = recordings.get(page);
  if (!r || !r.open)
    return;
  closeLeg(r, r.open.before, r.target ?? undefined, await drain(page));
  const stem = String(r.open.index).padStart(2, '0');
  r.legs.push({shot: await shot(page, r.dir, `${stem}-hop${r.legs.length + 1}.png`), target: box, pointer: []});
  r.located = undefined;
  r.target = box;
}

/** What the pointer did since the last stop goes to that stop; before the first stop, to a stop on
 * the page the step began with (when it acted there, or located something there). */
function closeLeg(r: Recording, before: string, target: GuideBox | undefined, drained: GuidePointer[]): void {
  r.pending.push(...drained);
  if (r.legs.length > 0)
    r.legs[r.legs.length - 1].pointer.push(...r.pending);
  else if (target || r.pending.length > 0)
    r.legs.push({shot: before, target, pointer: r.pending});
  r.pending = [];
}

/** The open step is set-up a person does not take (the login, a pinned setting, a wait for a
 * package), never part of a guide. */
export function silent(page: Page): void {
  const r = recordings.get(page);
  if (r?.open)
    r.silent = true;
}

function sameFile(dir: string, a: string, b: string): boolean {
  try {
    return readFileSync(join(dir, a)).equals(readFileSync(join(dir, b)));
  }
  catch {
    return false;
  }
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
  r = {dir, manifest, lastType: 'Given', keys: [], typed: '', legs: [], pending: [], dragShots: 0, silent: false};
  recordings.set(page, r);
  return r;
}

async function shot(page: Page, dir: string, name: string): Promise<string> {
  await page.screenshot({path: join(dir, name), animations: 'disabled', caret: 'hide'}).catch(() => undefined);
  return name;
}

/** Opens a step: the page as it is before the step, for the picture the pointer moves over. A step
 * before the page exists (the login) is recorded without pictures. */
export async function begin(page: Page | undefined, info: TestInfo, line: number, title: string, table?: string[][]): Promise<void> {
  if (!guideDir() || !page || page.isClosed())
    return;
  const r = recordingFor(page, info);
  const {keyword, text} = parseTitle(title);
  const index = r.manifest.steps.length + 1;
  const stem = String(index).padStart(2, '0');
  // what the pointer did between steps is part of neither
  await drain(page);
  r.keys = [];
  r.typed = '';
  r.located = undefined;
  r.target = undefined;
  r.legs = [];
  r.pending = [];
  r.dragShots = 0;
  r.silent = false;
  r.open = {index, line, keyword, text, table, started: Date.now(), before: await shot(page, r.dir, `${stem}-before.png`)};
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
  r.pending.push(...await drain(page));
  // a press nothing was located for lights the element it landed on, when that is not half the page
  const press = [...r.legs.flatMap((l) => l.pointer), ...r.pending].reverse().find((p) => p.type === 'down');
  const viewport = page.viewportSize();
  if (!target && press?.box && (!viewport || press.box.width * press.box.height < viewport.width * viewport.height / 2))
    target = press.box;
  else if (!target && press)
    target = {x: press.x - 12, y: press.y - 12, width: 24, height: 24};
  closeLeg(r, open.before, target, []);
  const type = typeOf(open.keyword, r.lastType);
  r.lastType = type;
  const acted = !!target || r.legs.some((l) => l.pointer.length > 0) || r.keys.length > 0 || r.typed.length > 0;
  // every step a reader would take is in the guide, a table opened through the API included (the
  // page after it is the point); left out are the login, a step that changed nothing on the page,
  // and a check or a wait a person has no use for
  const hidden = r.silent || ((type === 'Then' || type === 'Given') && hiddenInGuide(open.text));
  const kind: GuideStepKind = hidden ? 'setup' : type === 'Then' ? 'check' :
    acted || !sameFile(r.dir, open.before, after) ? 'action' : 'setup';
  r.manifest.steps.push({index: open.index, line: open.line, keyword: open.keyword, text: open.text,
    caption: captionOf(open.text, kind, open.table), kind, before: open.before, after, target, legs: r.legs,
    keys: r.keys, typed: r.typed, ms: Date.now() - open.started});
  writeFileSync(join(r.dir, 'steps.json'), JSON.stringify(r.manifest, null, 2));
}
