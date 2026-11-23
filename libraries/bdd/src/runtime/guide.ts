import {mkdirSync, readFileSync, writeFileSync} from 'node:fs';
import {join} from 'node:path';
import type {Locator, Page, TestInfo} from '@playwright/test';

export interface GuideBox {
  x: number;
  y: number;
  width: number;
  height: number;
}

/** What the pointer did, in the page's CSS pixels (the screenshots' pixels times the manifest's `scale`). */
export interface GuidePointer {
  /** `move`: where the pointer came to rest, or with `held` a point of a drag; `down` and `up`: a
   * button pressed or released there; `still`: a picture of the page while the step went on (text
   * being typed, a wheel turn), the pointer where it was. */
  type: 'move' | 'down' | 'up' | 'still';
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
  /** A move kept as a corner of the way (`waypoint`): the next move does not replace it. */
  kept?: boolean;
  /** A press the video shows and the page never got (`pressShown`): on the box it names, whatever is lit. */
  shown?: boolean;
  /** The keys held through a press (`Shift`, `Control`): a demo shows them as key caps. */
  mods?: string[];
}

/** `note`: a line of text the video shows over the page (`the video says …`), no gesture of its own. */
export type GuideStepKind = 'action' | 'check' | 'setup' | 'note';

const DRAG_SHOTS = 8;
/** At most this many pictures of one typing. */
const BEATS = 6;

/** A stop of a step: the page as the pointer set off for it, the element lit there, and what the
 * pointer did at it. A step has one stop, or one per stop of a path it walked (a menu's group, each
 * item, the leaf). The next stop's picture is the page as the pointer's arrival here left it (a group
 * opened, an item lit): the video cuts to it on arrival. */
export interface GuideLeg {
  shot: string;
  target?: GuideBox;
  pointer: GuidePointer[];
  /** `menu`: an item of a menu or a list on the way, held long enough to read; `list`: the control a list opens
   * from, the list held a little longer. */
  kind?: GuideLegKind;
}

export type GuideLegKind = 'menu' | 'list';

export interface GuideStep {
  index: number;
  line: number;
  keyword: string;
  text: string;
  /** The step as an instruction to a person: "Click the Open local file icon in the Browse toolbar". */
  caption: string;
  /** The caption is the feature's own (`# caption:`), not one made from the step's text: a demo shows it. */
  captioned?: boolean;
  kind: GuideStepKind;
  before: string;
  after: string;
  /** The element the step acted on or checked, in page pixels — the last stop when it walked a path. */
  target?: GuideBox;
  legs: GuideLeg[];
  /** Keys pressed during the step (`Enter`, `Control+A`), and text typed. */
  keys: string[];
  typed: string;
  /** The keys the step pressed as its own gesture (`user presses Escape`): a demo shows them as key caps. */
  pressed?: string[];
  ms: number;
  /** The view in front when the step ended: a demo folds the steps after an action into its result only while
   * they stay in its view. */
  view?: string;
  /** A dialog the step opened, by the inputs it shows: a demo rests on it as long as they take to read. */
  dialog?: {inputs: number};
  /** Where a note goes when the feature pins it (`the video says "…" at the top`); else the renderer places it. */
  noteAt?: 'top' | 'bottom';
}

export interface GuideManifest {
  feature: string;
  scenario: string;
  description: string;
  tags: string[];
  viewport: {width: number; height: number} | null;
  /** Screenshot pixels per CSS pixel: a small window is filmed at twice its density. */
  scale: number;
  /** The part of the page the video shows (`frameTo`): every picture is cut to it. */
  frame?: GuideBox;
  steps: GuideStep[];
}

interface Recording {
  dir: string;
  manifest: GuideManifest;
  lastType: string;
  keys: string[];
  typed: string;
  pressed: string[];
  located?: Locator;
  target?: GuideBox | null;
  legs: GuideLeg[];
  /** What the pointer did since the last stop. */
  pending: GuidePointer[];
  dragShots: number;
  beats: number;
  silent: boolean;
  /** The open step is a line of text for the video (`say`), pinned to the top or the bottom when the feature says. */
  note?: string;
  noteAt?: 'top' | 'bottom';
  /** A note was the last step shown since the last action: a hover's tooltip is no part of what follows it. */
  afterNote?: boolean;
  /** The step aimed at a hit area (`locatedArea`): finer than the viewer a later read of it locates. */
  area?: boolean;
  open?: {index: number; line: number; keyword: string; text: string; table?: string[][]; caption?: string; started: number; before: string};
}

const recordings = new WeakMap<Page, Recording>();
const attached = new WeakSet<Page>();

export function guideDir(): string | undefined {
  return process.env.BDD_GUIDE || undefined;
}

/** A help-page GIF (`grok-bdd guide --help-mode`): the bare shell, no captions, no set-up. */
export function helpMode(): boolean {
  return !!guideDir() && !!process.env.BDD_GUIDE_HELP;
}

/** A demo clip (`@demo`, `grok-bdd guide --demo`): the shell with its sidebar only, no captions, no set-up. */
export function demoMode(): boolean {
  return !!guideDir() && !!process.env.BDD_GUIDE_DEMO;
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
 * delivered it, and what the press landed on (with the keys held through it). A move with no button
 * held only says where the pointer came to rest, so a run of them keeps the last — but for a corner
 * of its way (`waypoint`); a synthetic `el.click()` has no place on the page and is not a gesture. */
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
      if (type === 'move' && !held && prev?.type === 'move' && !prev.held && !prev.kept)
        log.pop();
      const entry: any = {type, x: e.clientX, y: e.clientY};
      if (type === 'move' && held)
        entry.held = true;
      if (type !== 'move')
        entry.button = buttons[e.button] ?? 'left';
      if (type === 'down') {
        entry.count = 1;
        const mods = [e.ctrlKey && 'Control', e.metaKey && 'Meta', e.altKey && 'Alt', e.shiftKey && 'Shift']
          .filter(Boolean);
        if (mods.length > 0)
          entry.mods = mods;
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

/** A picture of the page in the middle of the step — the text typed so far, the list after a wheel
 * turn — shown where the pointer is when it was taken. */
async function beat(page: Page, r: Recording): Promise<void> {
  r.pending.push(...await drain(page));
  await page.evaluate(() => new Promise<void>((resolve) =>
    requestAnimationFrame(() => requestAnimationFrame(() => resolve())))).catch(() => undefined);
  const last = [...r.pending, ...r.legs.flatMap((l) => l.pointer)].reverse().find((p) => p.type !== 'still');
  r.pending.push({type: 'still', x: last?.x ?? 0, y: last?.y ?? 0,
    shot: await shot(page, r.dir, `${String(r.open!.index).padStart(2, '0')}-beat${++r.beats}.png`)});
}

/** The text in a few pieces, a picture after each: typing shows as it goes, not all at once. */
function pieces(text: string): string[] {
  const size = Math.max(1, Math.ceil(text.length / BEATS));
  const out: string[] = [];
  for (let i = 0; i < text.length; i += size)
    out.push(text.slice(i, i + size));
  return out;
}

/** A locator's typing goes past `page.keyboard`: its `fill` and `pressSequentially` are wrapped once
 * for every page, and do what they did whenever the page is not being filmed. */
let locatorWrapped = false;
/** The locator's own `fill`, for a text pictured going in that a gesture then puts in at once (`pictureTyping`). */
let rawFill: ((this: Locator, value: string, options?: unknown) => Promise<void>) | undefined;
function wrapLocatorTyping(loc: Locator): void {
  if (locatorWrapped)
    return;
  locatorWrapped = true;
  const proto = Object.getPrototypeOf(loc);
  const fill = rawFill = proto.fill;
  const pressSequentially = proto.pressSequentially;
  const filming = (l: Locator): Recording | undefined => {
    const r = recordings.get(l.page());
    return r?.open && !r.silent ? r : undefined;
  };
  proto.fill = async function(this: Locator, value: string, options?: unknown) {
    const r = filming(this);
    if (!r || !value)
      return fill.call(this, value, options);
    let typed = '';
    for (const piece of pieces(value)) {
      typed += piece;
      await fill.call(this, typed, options);
      await beat(this.page(), r);
    }
  };
  proto.pressSequentially = async function(this: Locator, text: string, options?: unknown) {
    const r = filming(this);
    if (!r)
      return pressSequentially.call(this, text, options);
    for (const piece of pieces(text)) {
      await pressSequentially.call(this, piece, options);
      await beat(this.page(), r);
    }
  };
}

/** Rows as wide as their panel, pressed where a person presses them — on their caption: an accordion pane's header, a
 * tree node, a menu item. */
const CAPTIONED = '.d4-accordion-pane-header, .d4-tree-view-node, .d4-menu-item, .u2-menu-item, [role="treeitem"], ' +
  '[role="menuitem"]';

/** A locator's click, double-click and hover land on the middle of the element: a row as wide as its panel is pressed
 * far right of its caption (an accordion pane's header 250 px right of "Explore"). A guide aims a row (`CAPTIONED`) at
 * its caption — the middle of its first text — when that is well off the middle and nothing of its own (a link, a
 * button, a field) is there; a test, and a gesture given a position, are as they were. */
let clicksWrapped = false;
function wrapLocatorClicks(loc: Locator): void {
  if (clicksWrapped)
    return;
  clicksWrapped = true;
  const proto = Object.getPrototypeOf(loc);
  for (const name of ['click', 'dblclick', 'hover'] as const) {
    const gesture = proto[name];
    proto[name] = async function(this: Locator, options?: {position?: {x: number; y: number}}) {
      const r = recordings.get(this.page());
      if (!r?.open || r.silent || options?.position)
        return gesture.call(this, options);
      // looked up only on an element already there (the gesture keeps its own wait); no named functions in the page:
      // tsx wraps them in a `__name` call the page lacks
      const position = await this.first().evaluate((el, rows) => {
        const box = el.getBoundingClientRect();
        if (!el.matches(rows) || box.width < 3 * box.height)
          return null;
        const walker = document.createTreeWalker(el, NodeFilter.SHOW_TEXT);
        for (let n = walker.nextNode(); n; n = walker.nextNode()) {
          if (!(n.textContent ?? '').trim())
            continue;
          const range = document.createRange();
          range.selectNodeContents(n);
          const t = range.getBoundingClientRect();
          if (t.width === 0 || t.height === 0)
            continue;
          const x = Math.max(box.left + 2, Math.min(box.right - 2, t.left + t.width / 2));
          const y = box.top + box.height / 2;
          const hit = document.elementFromPoint(x, y);
          const own = hit?.closest('a, button, input, select, textarea, [role="button"], [role="link"], ' +
            '[role="checkbox"]');
          if (Math.abs(x - (box.left + box.width / 2)) < 40 || !hit || !el.contains(hit) ||
            (own && own !== el && el.contains(own)))
            return null;
          return {x: x - box.left, y: y - box.top};
        }
        return null;
      }, CAPTIONED, {timeout: 500}).catch(() => null);
      return gesture.call(this, position ? {...options, position} : options);
    };
  }
}

/** The pointer is followed in the page, whatever runtime path drives it (a gesture, a menu walk, a
 * viewer hit area, a locator's own click). The page's `mouse` is wrapped to picture drags and wheel
 * turns, the typing to picture the text as it goes in, and a locator's clicks and hover to press a row
 * on its caption. */
export async function attach(page: Page): Promise<void> {
  if (!guideDir() || attached.has(page))
    return;
  attached.add(page);
  wrapLocatorTyping(page.locator('body'));
  wrapLocatorClicks(page.locator('body'));
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
    const legs = Math.min(DRAG_SHOTS - r.dragShots, Math.max(2, Math.round(Math.hypot(x - from.x, y - from.y) / 60)));
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
  // one event, as a test sends it: split into notches, a handler that moves a step per event (Bio's MSA header)
  // went further on film than in the test
  const wheel = mouse.wheel.bind(mouse);
  mouse.wheel = async (dx: number, dy: number) => {
    await wheel(dx, dy);
    const r = rec();
    if (r?.open && !r.silent)
      await beat(page, r);
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
  // a native select finds an option by the keys typed within a second: only a text field is pictured as it fills
  const intoText = (): Promise<boolean> => page.evaluate(() => {
    const e = document.activeElement as HTMLElement | null;
    return !!e && (e.tagName === 'INPUT' || e.tagName === 'TEXTAREA' || e.isContentEditable);
  }).catch(() => false);
  keyboard.type = async (text: string, options?: unknown) => {
    const r = rec();
    if (r)
      r.typed += text;
    if (!r?.open || r.silent || !await intoText())
      return type(text, options);
    for (const piece of pieces(text)) {
      await type(piece, options);
      await beat(page, r);
    }
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
  // an area the step aims at stays lit: the viewer located after it (a settle) is its scope
  if (!r || !r.open || r.area)
    return;
  r.located = loc;
  r.target = await loc.first().boundingBox({timeout: 200}).catch(() => null);
  // an element below the fold (a property category low in the context panel) is scrolled into view by the gesture
  // that follows, so the box read here is not where the press lands: filmed, the element is scrolled first and the
  // step's first picture taken again. Only a small element off the window, before the step has any stop: a scope
  // (a whole panel) is never scrolled
  const view = page.viewportSize();
  const box = r.target;
  if (box && view && r.legs.length === 0 && box.height < view.height / 2 &&
    (box.y < 0 || box.y + box.height > view.height)) {
    await loc.first().scrollIntoViewIfNeeded({timeout: 1000}).catch(() => undefined);
    r.target = await loc.first().boundingBox({timeout: 200}).catch(() => null);
    await shot(page, r.dir, r.open.before);
  }
}

/** The step acts on a place inside the element it located — a hit area of a viewer or a widget (an
 * atom of a sketcher, a grid cell): that place is lit, not the element round it, and a small one is
 * zoomed into as an icon is. */
export function locatedArea(page: Page, box: GuideBox): void {
  const r = recordings.get(page);
  if (!r || !r.open)
    return;
  r.located = undefined;
  r.target = box;
  r.area = true;
}

/** The video shows this part of the page only (a sketcher's dialog, for a guide about the sketcher):
 * set by a set-up step, it holds for the scenario, and every picture is cut to it — the box grown by
 * `margin` and kept on the page. */
export function frameTo(page: Page, box: GuideBox, margin = 16): void {
  const r = recordings.get(page);
  if (!r || !r.open)
    return;
  const view = page.viewportSize() ?? {width: Infinity, height: Infinity};
  const x = Math.max(0, Math.floor(box.x - margin));
  const y = Math.max(0, Math.floor(box.y - margin));
  const right = Math.min(view.width, Math.ceil(box.x + box.width + margin));
  const bottom = Math.min(view.height, Math.ceil(box.y + box.height + margin));
  r.manifest.frame = {x, y, width: (right - x) & ~1, height: (bottom - y) & ~1};
}

/** A stop on the path the step walks: the page as it is now (a menu opened by the stop before)
 * and the element's place on it; the guide moves the pointer stop by stop and lights each. What
 * the pointer did before the first stop — on an element located before it, or anywhere (the
 * right-click that opened the menu) — is a stop of its own, on the "before" picture. A stop is
 * taken as the pointer sets off for it: the next stop's picture is what the pointer's arrival led
 * to, and the video cuts to it on arrival. `kind` paces it (`GuideLeg`). */
export async function hop(page: Page, loc: Locator, kind?: GuideLegKind): Promise<void> {
  if (!recordings.get(page)?.open)
    return;
  const box = await loc.first().boundingBox({timeout: 200}).catch(() => null);
  if (box)
    await hopAt(page, box, kind);
}

/** The leaf of a menu walk, before the click a caller makes: a stop as the pointer sets off for it, then — the
 * pointer on it — a stop with the item lit under the pointer, so the video shows the item the press lands on lit
 * before the press. An item the pointer already rests on is lit in its first stop. A test only clicks. */
export async function hopLeaf(page: Page, leaf: Locator, timeout?: number): Promise<void> {
  await hop(page, leaf, 'menu');
  if (!filming(page) || await leaf.first().evaluate((e) => e.matches(':hover')).catch(() => false))
    return;
  await leaf.first().hover({timeout});
  await hop(page, leaf, 'menu');
}

/** Whether the page's open step is being filmed (a guide, not its set-up). */
export function filming(page: Page): boolean {
  const r = recordings.get(page);
  return !!r?.open && !r.silent;
}

/** Keeps the pointer's last move as a corner of its way: a run of moves with no button held is logged as the one
 * point it came to rest at, and a way that goes round something (into a menu's flyout along the group's row, not
 * across the items below it) would be filmed as a straight line. */
export async function waypoint(page: Page): Promise<void> {
  if (!filming(page))
    return;
  await page.evaluate(() => {
    const log = (window as any).__bddGuidePointer as {type: string; held?: boolean; kept?: boolean}[] | undefined;
    const last = log?.[log.length - 1];
    if (last?.type === 'move' && !last.held)
      last.kept = true;
  }).catch(() => undefined);
}

/** Pointer events the video shows where the step sends none: logged as the page logs its own. */
async function logShown(page: Page, events: GuidePointer[]): Promise<void> {
  if (!filming(page))
    return;
  await page.evaluate((list) => {
    (window as any).__bddGuidePointer?.push(...list);
  }, events).catch(() => undefined);
}

/** The pointer moved on film to the middle of a place the step moves nothing to for real (a row of a list only the
 * video shows). */
export async function pointShown(page: Page, box: GuideBox): Promise<void> {
  await logShown(page, [{type: 'move', x: box.x + box.width / 2, y: box.y + box.height / 2, kept: true}]);
}

/** A press the video shows where the step presses nothing for real: the place is pressed on film, the page gets no
 * event (a native select's list, which headless Chromium does not draw, opened and an option of it taken). Its box
 * is the element pressed, whatever the stop lit. */
export async function pressShown(page: Page, box: GuideBox, on: string): Promise<void> {
  const x = box.x + box.width / 2;
  const y = box.y + box.height / 2;
  await logShown(page, [{type: 'move', x, y}, {type: 'down', x, y, button: 'left', count: 1, box, on, shown: true},
    {type: 'up', x, y, button: 'left'}]);
}

/** A move the video leaves out: the pointer goes there for real — off a canvas whose lens or tooltip follows it, which
 * would otherwise stay drawn where the pointer rested while the video's pointer goes on — and the video's pointer goes
 * on from where it was. */
export async function moveUnfilmed(page: Page, x: number, y: number): Promise<void> {
  const r = recordings.get(page);
  if (!r?.open || r.silent)
    return;
  r.pending.push(...await drain(page));
  await page.mouse.move(x, y);
  await drain(page);
}

/** A picture of the page in the middle of the step, the pointer where it is (a list scrolling under it). */
export async function still(page: Page): Promise<void> {
  const r = recordings.get(page);
  if (r?.open && !r.silent)
    await beat(page, r);
}

/** A key the step presses as its gesture (`user presses Escape`), whichever way it reaches the page: a demo shows it
 * as a key cap. */
export function keyPressed(page: Page, key: string): void {
  const r = recordings.get(page);
  if (r?.open && !r.silent)
    r.pressed.push(key);
}

/** A text a gesture puts in at once (a paste), pictured going in as a typing is: the editor filled piece by piece
 * short of the whole text, a picture after each; the gesture then puts the whole of it in itself. A single-line
 * field shows the lines of the text as values separated by commas. */
export async function pictureTyping(page: Page, loc: Locator, text: string): Promise<void> {
  const r = recordings.get(page);
  if (!r?.open || r.silent || !rawFill)
    return;
  const single = await loc.evaluate((e) => e.tagName === 'INPUT').catch(() => true);
  const shown = single ? text.replace(/\n+$/, '').replace(/\n+/g, ', ') : text;
  let typed = '';
  for (const piece of pieces(shown).slice(0, -1)) {
    typed += piece;
    if (!await rawFill.call(loc, typed).then(() => true, () => false))
      return;
    await beat(page, r);
  }
}

/** A guide films a column picker with no tooltip over it: one the pointer raised beside the picker (a column
 * header's) covers the row the picker is about to take. On while the picker is filmed, off after it. */
export async function pickerFilmed(page: Page, open: boolean): Promise<void> {
  if (!filming(page))
    return;
  await page.evaluate((on) => {
    if (!document.getElementById('bdd-picker-style')) {
      const style = document.createElement('style');
      style.id = 'bdd-picker-style';
      style.textContent = 'body.bdd-picker-open .d4-tooltip {display: none !important}';
      document.head.append(style);
    }
    document.body.classList.toggle('bdd-picker-open', on);
  }, open).catch(() => undefined);
}

/** The same stop for a place that is no element of its own — a row of a canvas grid. */
export async function hopAt(page: Page, box: GuideBox, kind?: GuideLegKind): Promise<void> {
  const r = recordings.get(page);
  if (!r || !r.open)
    return;
  // an element located before the walk is no stop of its own unless the pointer did something there
  closeLeg(r, r.open.before, r.target ?? undefined, await drain(page), false);
  const stem = String(r.open.index).padStart(2, '0');
  r.legs.push({shot: await shot(page, r.dir, `${stem}-hop${r.legs.length + 1}.png`), target: box, pointer: [],
    ...(kind ? {kind} : {})});
  r.located = undefined;
  r.target = box;
}

/** What the pointer did since the last stop goes to that stop; before the first stop, to a stop on
 * the page the step began with (when it acted there, or — `located` — located something there). */
function closeLeg(r: Recording, before: string, target: GuideBox | undefined, drained: GuidePointer[],
  located = true): void {
  r.pending.push(...drained);
  if (r.legs.length > 0)
    r.legs[r.legs.length - 1].pointer.push(...r.pending);
  else if ((located && target) || r.pending.length > 0)
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

/** The open step is a line of text the video shows over the page, whatever the step's own words; `at` pins it to the
 * top or the bottom of the page, where the renderer would otherwise choose. */
export function say(page: Page, text: string, at?: 'top' | 'bottom'): void {
  const r = recordings.get(page);
  if (r?.open) {
    r.note = text;
    r.noteAt = at;
  }
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
    {feature: featureName, scenario: scenarioName, description: '', tags: info.tags, viewport,
      scale: info.project.use.deviceScaleFactor ?? 1, steps: []};
  r = {dir, manifest, lastType: 'Given', keys: [], typed: '', pressed: [], legs: [], pending: [], dragShots: 0,
    beats: 0, silent: false};
  recordings.set(page, r);
  return r;
}

async function shot(page: Page, dir: string, name: string): Promise<string> {
  await page.screenshot({path: join(dir, name), animations: 'disabled', caret: 'hide'}).catch(() => undefined);
  return name;
}

/** Opens a step: the page as it is before the step, for the picture the pointer moves over. A step
 * before the page exists (the login) is recorded without pictures. */
export async function begin(page: Page | undefined, info: TestInfo, line: number, title: string, table?: string[][],
  caption?: string): Promise<void> {
  if (!guideDir() || !page || page.isClosed())
    return;
  const r = recordingFor(page, info);
  const {keyword, text} = parseTitle(title);
  const index = r.manifest.steps.length + 1;
  const stem = String(index).padStart(2, '0');
  // what the pointer did between steps is part of neither
  await drain(page);
  // the dialogs open before the step: one that is open after it and not among them is the step's (`end`)
  await page.evaluate((sel) => document.querySelectorAll(sel).forEach((d) => { (d as any).__bddSeen = true; }),
    DIALOGS).catch(() => undefined);
  r.keys = [];
  r.typed = '';
  r.pressed = [];
  r.located = undefined;
  r.target = undefined;
  r.legs = [];
  r.pending = [];
  r.dragShots = 0;
  r.beats = 0;
  r.silent = false;
  r.note = undefined;
  r.noteAt = undefined;
  r.area = false;
  // a demo's set-up after a note can change the page under a tooltip a hover left up (the rows switched under it): the
  // tooltip goes first. Right after a hover the set-up is part of what the hover led to, and the tooltip stays
  if (demoMode() && r.afterNote && typeOf(keyword, r.lastType) === 'Given')
    await page.evaluate(() => (window as any).ui?.tooltip?.hide?.()).catch(() => undefined);
  r.open = {index, line, keyword, text, table, caption, started: Date.now(), before: await shot(page, r.dir, `${stem}-before.png`)};
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
  const view = await page.evaluate(() => (window as any).grok?.shell?.v?.name as string | undefined).catch(() => undefined);
  const inputs = await newDialogInputs(page);
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
  // and a check or a wait a person has no use for — unless the feature gave it a caption to show
  const hidden = r.silent || (!open.caption && (type === 'Then' || type === 'Given') && hiddenInGuide(open.text));
  const kind: GuideStepKind = r.note !== undefined ? 'note' : hidden ? 'setup' : type === 'Then' ? 'check' :
    acted || !sameFile(r.dir, open.before, after) ? 'action' : 'setup';
  const caption = r.note ?? open.caption;
  r.afterNote = kind === 'note' ? true : kind === 'action' ? false : r.afterNote;
  r.manifest.steps.push({index: open.index, line: open.line, keyword: open.keyword, text: open.text,
    caption: caption ?? captionOf(open.text, kind, open.table), ...(caption ? {captioned: true} : {}), kind,
    before: open.before, after, target, legs: r.legs, keys: r.keys, typed: r.typed,
    ...(r.pressed.length > 0 ? {pressed: r.pressed} : {}), ms: Date.now() - open.started,
    ...(view ? {view} : {}), ...(inputs !== undefined ? {dialog: {inputs}} : {}),
    ...(r.noteAt ? {noteAt: r.noteAt} : {})});
  writeFileSync(join(r.dir, 'steps.json'), JSON.stringify(r.manifest, null, 2));
}

const DIALOGS = '.d4-dialog, [data-u2="dialog"]';

/** The inputs a dialog the step opened shows (its input rows; else its bare fields), or nothing when no dialog is
 * open that was not open before the step. */
function newDialogInputs(page: Page): Promise<number | undefined> {
  // no named functions in the page: tsx wraps them in a `__name` call the page lacks
  return page.evaluate((sel) => {
    const fresh = Array.from(document.querySelectorAll(sel)).filter((d) => !(d as any).__bddSeen &&
      d.getClientRects().length > 0 && getComputedStyle(d).visibility !== 'hidden');
    const dialog = fresh.filter((d) => !fresh.some((o) => o !== d && o.contains(d))).pop();
    if (!dialog)
      return undefined;
    for (const query of ['.ui-input-root, [data-u2-part="editor"]', 'input, select, textarea']) {
      const found = Array.from(dialog.querySelectorAll(query)).filter((e) => e.getClientRects().length > 0 &&
        getComputedStyle(e).visibility !== 'hidden');
      const rows = found.filter((e) => !found.some((o) => o !== e && o.contains(e))).length;
      if (rows > 0)
        return rows;
    }
    return 0;
  }, DIALOGS).catch(() => undefined);
}
