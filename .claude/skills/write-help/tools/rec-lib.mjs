// Help-GIF recorder: signs in with a session token from 'grok s token' (nothing secret is printed),
// records a CDP screencast, draws the pointer with click marks (yellow = left, green = right), and
// encodes a GIF, a -thumb.png and a contact sheet with ffmpeg. See ../gifs.md.
import {createRequire} from 'module';
import {fileURLToPath} from 'url';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {execFileSync, spawnSync} from 'child_process';

// Playwright comes from the bdd library of the public repo (this file sits in .claude/skills/write-help/tools);
// REC_PLAYWRIGHT_FROM names another package.json to resolve it from.
const HERE = path.dirname(fileURLToPath(import.meta.url));
const require = createRequire(process.env.REC_PLAYWRIGHT_FROM ?? path.join(HERE, '../../../../libraries/bdd/package.json'));
const {chromium} = (() => {
  try { return require('@playwright/test'); }
  catch (_) { return require('playwright'); }
})();

export const HOST = process.env.REC_HOST ?? 'dev';
export const W = +(process.env.REC_W ?? 800), H = +(process.env.REC_H ?? Math.round(W * 5 / 8));
// Hardware WebGL: software rendering stalls 3D viewers. The screencast ignores the context's deviceScaleFactor and
// sends 1x frames unless the scale is forced. REC_BROWSER_ARGS overrides (space-separated).
const BROWSER_ARGS = process.env.REC_BROWSER_ARGS?.split(' ').filter(Boolean) ??
  [...(process.platform === 'win32' ? ['--use-angle=d3d11'] : []), '--enable-gpu', '--ignore-gpu-blocklist',
    '--force-device-scale-factor=2'];

// A full ffmpeg build (Playwright's own writes WebM only): FFMPEG, else imageio-ffmpeg's binary, else ffmpeg on PATH.
function findFfmpeg() {
  if (process.env.FFMPEG) return process.env.FFMPEG;
  for (const py of ['python3', 'python', 'py']) {
    try { return execFileSync(py, ['-c', 'import imageio_ffmpeg;print(imageio_ffmpeg.get_ffmpeg_exe())'], {stdio: ['ignore', 'pipe', 'ignore']}).toString().trim(); }
    catch (_) {}
  }
  return 'ffmpeg';
}
const FFMPEG = findFfmpeg();

// The server URL comes from the grok CLI config; the session token from 'grok s token', so no key is read here.
function serverUrl(alias) {
  const lines = fs.readFileSync(path.join(os.homedir(), '.grok', 'config.yaml'), 'utf8').split(/\r?\n/);
  const i = lines.findIndex((l) => l.trim() === `${alias}:`);
  if (i < 0) throw new Error(`server alias '${alias}' is not in ~/.grok/config.yaml (set REC_HOST)`);
  const indent = lines[i].search(/\S/);
  for (let j = i + 1; j < lines.length && (lines[j].trim() === '' || lines[j].search(/\S/) > indent); j++) {
    const m = lines[j].trim().match(/^url:\s*['"]?([^'"\s]+)/);
    if (m) return m[1];
  }
  throw new Error(`no url for alias '${alias}' in ~/.grok/config.yaml`);
}

function sessionToken(alias) {
  const r = spawnSync(process.env.REC_GROK ?? 'grok', ['s', 'token', '--host', alias],
    {encoding: 'utf8', shell: process.platform === 'win32'});
  const token = (r.stdout ?? '').trim().split(/\r?\n/).pop();
  if (r.status !== 0 || !token) throw new Error(`'grok s token --host ${alias}' failed: ${(r.stderr ?? '').trim().slice(0, 200)}`);
  return token;
}

export async function open() {
  const api = serverUrl(HOST);
  const origin = api.replace(/\/api\/?$/, '');
  const token = sessionToken(HOST);
  const browser = await chromium.launch({headless: true, channel: 'chromium', args: BROWSER_ARGS});
  const context = await browser.newContext({viewport: {width: W, height: H}, deviceScaleFactor: 2});
  await context.addCookies([{name: 'auth', value: token, url: origin}]);
  await context.addInitScript((t) => { try { localStorage.setItem('auth', t); } catch (_) {} }, token);
  const page = await context.newPage();
  await page.goto(origin, {waitUntil: 'domcontentloaded'});
  await page.waitForFunction(() => { try { return window.grok?.shell?.user?.login != null && document.querySelector('.layout-workarea, .d4-root') != null; } catch (_) { return false; } }, null, {timeout: 120000, polling: 1000});
  await page.waitForTimeout(3000);
  return {browser, context, page, origin};
}

/** Closes the panels a GIF must not show and resets the workspace. */
export async function cleanShell(page) {
  await page.evaluate(() => {
    grok.shell.closeAll();
    const w = grok.shell.windows;
    w.showBrowse = false; w.showToolbox = false; w.showContextPanel = false; w.showHelp = false;
    w.showConsole = false; w.showVariables = false;
    try { w.simpleMode = false; } catch (_) {}
  });
  await page.waitForTimeout(800);
}

export async function installPointer(page) {
  await page.evaluate(() => {
    if (window.__recPointer) return;
    window.__recPointer = true;
    const arrow = document.createElement('div');
    arrow.innerHTML = '<svg width="16" height="24" viewBox="0 0 14 21"><path d="M1 1 L1 17 L5 13.5 L8 20.5 L11 19 L8 12.5 L13.5 12.5 Z" fill="white" stroke="black" stroke-width="1.2"/></svg>';
    Object.assign(arrow.style, {position: 'fixed', left: '0', top: '0', zIndex: 2147483647, pointerEvents: 'none',
      transition: 'transform 60ms linear', transform: 'translate(-100px,-100px)'});
    const disc = document.createElement('div');
    Object.assign(disc.style, {position: 'fixed', left: '0', top: '0', width: '26px', height: '26px', boxSizing: 'border-box',
      borderRadius: '50%', zIndex: 2147483646, pointerEvents: 'none', display: 'none',
      transition: 'transform 60ms linear', opacity: '0.75', border: '2px solid rgba(0,0,0,0.35)'});
    // Top-layer popovers: a select picker or another popover is drawn in the top layer above any z-index, so the
    // pointer is a popover too, raised over whatever top-layer element opens last.
    for (const el of [disc, arrow]) {
      el.popover = 'manual';
      Object.assign(el.style, {inset: 'auto', left: '0', top: '0', margin: '0', overflow: 'visible'});
    }
    Object.assign(arrow.style, {padding: '0', border: 'none', background: 'transparent'});
    Object.assign(disc.style, {padding: '0'});
    document.documentElement.append(disc, arrow);
    const raise = () => { for (const el of [disc, arrow]) { try { el.hidePopover(); } catch (_) {} el.showPopover(); } };
    raise();
    let topLayer = '';
    const watch = () => {
      const now = [...document.querySelectorAll(':popover-open, :modal')].filter((e) => e !== disc && e !== arrow).length +
        ':' + document.querySelectorAll('select:open').length;
      if (now !== topLayer) { topLayer = now; raise(); }
      requestAnimationFrame(watch);
    };
    requestAnimationFrame(watch);
    let last = 0, down = false, hideTimer = null, px = -100, py = -100;
    const move = (e) => {
      // a long step (Mouse.press) jumps, or the arrow trails the press mark for a frame
      const jump = Math.hypot(e.clientX - px, e.clientY - py) > 40;
      px = e.clientX; py = e.clientY;
      for (const el of [arrow, disc])
        el.style.transition = jump ? 'none' : 'transform 60ms linear';
      arrow.style.transform = `translate(${e.clientX - 1}px,${e.clientY - 1}px)`;
      disc.style.transform = `translate(${e.clientX - 13}px,${e.clientY - 13}px)`;
    };
    const press = (e) => {
      const now = performance.now();
      if (down && now - last < 30) return;
      last = now; down = true;
      disc.style.background = e.button === 2 ? 'rgb(34,197,94)' : 'rgb(255,214,0)';
      move(e);
      clearTimeout(hideTimer);
      disc.style.display = 'block';
    };
    const release = () => { down = false; clearTimeout(hideTimer); hideTimer = setTimeout(() => { if (!down) disc.style.display = 'none'; }, 250); };
    for (const t of ['pointermove', 'mousemove']) window.addEventListener(t, move, true);
    for (const t of ['pointerdown', 'mousedown']) window.addEventListener(t, press, true);
    for (const t of ['pointerup', 'mouseup']) window.addEventListener(t, release, true);
  });
}

export class Mouse {
  constructor(page) { this.page = page; this.x = W / 2; this.y = H / 2; }
  async moveTo(x, y, speed = 900) {
    const d = Math.hypot(x - this.x, y - this.y);
    const steps = Math.max(4, Math.round(d / speed * 30));
    for (let i = 1; i <= steps; i++) {
      const t = i / steps, e = t * t * (3 - 2 * t);
      await this.page.mouse.move(this.x + (x - this.x) * e, this.y + (y - this.y) * e);
      await this.page.waitForTimeout(1000 / 30);
    }
    this.x = x; this.y = y;
  }
  async click(x, y, {button = 'left', pause = 350, after = 500} = {}) {
    await this.moveTo(x, y);
    await this.page.waitForTimeout(pause);
    await this.page.mouse.down({button});
    await this.page.waitForTimeout(160);
    await this.page.mouse.up({button});
    await this.page.waitForTimeout(after);
  }
  /** Steps onto (x, y) in one move and presses at once: for targets that react on mouseenter, such as the groups of
   *  a horizontal menu bar, so the press mark appears in the same frame as the reaction. */
  async press(x, y, {button = 'left', after = 500} = {}) {
    await this.page.mouse.move(x, y);
    this.x = x; this.y = y;
    await this.page.mouse.down({button});
    await this.page.waitForTimeout(160);
    await this.page.mouse.up({button});
    await this.page.waitForTimeout(after);
  }
  async dblclick(x, y) {
    await this.moveTo(x, y);
    await this.page.waitForTimeout(300);
    for (let i = 0; i < 2; i++) {
      await this.page.mouse.down({clickCount: i + 1});
      await this.page.waitForTimeout(70);
      await this.page.mouse.up({clickCount: i + 1});
      await this.page.waitForTimeout(90);
    }
    await this.page.waitForTimeout(500);
  }
  async drag(x1, y1, x2, y2, {steps = 25, stepMs = 50} = {}) {
    await this.moveTo(x1, y1);
    await this.page.waitForTimeout(300);
    await this.page.mouse.down();
    await this.page.waitForTimeout(200);
    for (let i = 1; i <= steps; i++) {
      this.x = x1 + (x2 - x1) * i / steps; this.y = y1 + (y2 - y1) * i / steps;
      await this.page.mouse.move(this.x, this.y);
      await this.page.waitForTimeout(stepMs);
    }
    await this.page.waitForTimeout(250);
    await this.page.mouse.up();
    await this.page.waitForTimeout(500);
  }
  async center(locator) {
    const b = await locator.first().boundingBox();
    if (!b) throw new Error(`no box for ${locator}`);
    return [b.x + b.width / 2, b.y + b.height / 2];
  }
  async clickEl(locator, opts) { const [x, y] = await this.center(locator); await this.click(x, y, opts); }
  async hoverEl(locator) { const [x, y] = await this.center(locator); await this.moveTo(x, y); }
  async type(text, delay = 90) { await this.page.keyboard.type(text, {delay}); }
}

/** Width of a baseline or progressive JPEG, read from its SOF marker. */
function jpegWidth(file) {
  const b = fs.readFileSync(file);
  for (let i = 2; i < b.length;) {
    const marker = b[i + 1], len = b.readUInt16BE(i + 2);
    if (marker >= 0xc0 && marker <= 0xc2) return b.readUInt16BE(i + 7);
    i += 2 + len;
  }
  return W;
}

export class Recorder {
  constructor(page, dir) { this.page = page; this.dir = dir; this.frames = []; }
  async start(nudge = [W / 2 + 1, H / 2 + 1]) {
    fs.rmSync(this.dir, {recursive: true, force: true});
    fs.mkdirSync(this.dir, {recursive: true});
    this.cdp = await this.page.context().newCDPSession(this.page);
    this.cdp.on('Page.screencastFrame', async (f) => {
      const file = path.join(this.dir, `f${String(this.frames.length).padStart(5, '0')}.jpg`);
      fs.writeFileSync(file, Buffer.from(f.data, 'base64'));
      this.frames.push({file, t: f.metadata.timestamp});
      try { await this.cdp.send('Page.screencastFrameAck', {sessionId: f.sessionId}); } catch (_) {}
    });
    await this.cdp.send('Page.startScreencast', {format: 'jpeg', quality: 95, maxWidth: W * 2, maxHeight: H * 2});
    // a still page emits no frames: nudge a repaint so the first frame exists. Frames from before the nudge show the
    // pointer where the setup left it, so they are dropped
    await this.page.waitForTimeout(200);
    const before = this.frames.length;
    await this.page.mouse.move(nudge[0], nudge[1]);
    await this.page.waitForTimeout(600);
    if (this.frames.length > before)
      this.frames.splice(0, before);
  }
  async stop() {
    await this.page.waitForTimeout(300);
    await this.cdp.send('Page.stopScreencast');
    // a still page sends no frames, so the last frame lasts until the stop: a final wait in the scene is kept
    this.end = this.frames.length ? Math.max(this.frames[this.frames.length - 1].t + 0.1, Date.now() / 1000) : 0;
  }
  /** Writes <out>.gif and <out>-thumb.png; thumbAt is a fraction of the timeline. Without size: 800×500, or with a
   *  crop, at most 800 px wide with the crop's aspect ratio. */
  make(out, {fps = 12, thumbAt = 0.95, colors = 128, crop = null, size = null} = {}) {
    if (!size) {
      const w = crop ? Math.min(800, Math.round(crop[2])) : 800;
      size = crop ? [w & ~1, Math.round(w * crop[3] / crop[2]) & ~1] : [800, 500];
    }
    const list = path.join(this.dir, 'list.txt');
    const lines = [];
    this.frames.forEach((f, i) => {
      const next = i + 1 < this.frames.length ? this.frames[i + 1].t : this.end;
      lines.push(`file '${f.file.replace(/\\/g, '/')}'`, `duration ${Math.max(0.001, next - f.t).toFixed(3)}`);
    });
    lines.push(`file '${this.frames[this.frames.length - 1].file.replace(/\\/g, '/')}'`);
    fs.writeFileSync(list, lines.join('\n'));
    const k = jpegWidth(this.frames[0].file) / W;
    if (k < 2)
      console.warn(`frames are ${k}x the viewport, not 2x: the GIF is upscaled`);
    const cropF = crop ? `crop=${Math.round(crop[2] * k)}:${Math.round(crop[3] * k)}:${Math.round(crop[0] * k)}:${Math.round(crop[1] * k)},` : '';
    const vf = `fps=${fps},${cropF}scale=${size[0]}:${size[1]}:flags=lanczos,split[a][b];[a]palettegen=max_colors=${colors}:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle`;
    execFileSync(FFMPEG, ['-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', list, '-vf', vf, '-loop', '0', `${out}.gif`]);
    const tAt = this.frames[0].t + (this.end - this.frames[0].t) * thumbAt;
    const thumb = [...this.frames].reverse().find((f) => f.t <= tAt) ?? this.frames[0];
    execFileSync(FFMPEG, ['-y', '-loglevel', 'error', '-i', thumb.file, '-vf', `${cropF}scale=${size[0]}:${size[1]}:flags=lanczos`, `${out}-thumb.png`]);
    const kb = Math.round(fs.statSync(`${out}.gif`).size / 1024);
    const dur = (this.end - this.frames[0].t).toFixed(1);
    console.log(`gif ${out}.gif: ${size[0]}x${size[1]}, ${this.frames.length} frames, ${dur}s, ${kb} KB`);
  }
  /** A contact sheet of evenly spaced frames, for reviewing a take without watching it. */
  sheet(out, n = 12) {
    const pick = [];
    for (let i = 0; i < n; i++) pick.push(this.frames[Math.floor(i * (this.frames.length - 1) / (n - 1))].file);
    const args = ['-y', '-loglevel', 'error'];
    for (const p of pick) args.push('-i', p);
    const filt = pick.map((_, i) => `[${i}:v]scale=400:250[v${i}]`).join(';') + ';' +
      pick.map((_, i) => `[v${i}]`).join('') + `xstack=inputs=${n}:layout=` +
      pick.map((_, i) => `${(i % 4) * 400}_${Math.floor(i / 4) * 250}`).join('|') + '[out]';
    execFileSync(FFMPEG, [...args, '-filter_complex', filt, '-map', '[out]', out]);
  }
}

/** Visible menu item by its exact label. */
export function menuItem(page, label) {
  return page.locator('.d4-menu-item-label:visible').getByText(label, {exact: true}).first();
}

/** Walks a menu path with the mouse: hover each group, click the leaf. A group in a horizontal menu bar opens on
 *  mouseenter, so the pointer approaches it from below (never sliding along the bar over its neighbours), stops just
 *  outside it, and steps in with the press: the press mark and the opening menu land in the same frame. */
export async function walkMenu(mouse, labels) {
  let prevHorz = false;
  for (let i = 0; i < labels.length; i++) {
    const item = menuItem(mouse.page, labels[i]);
    const [x, y] = await mouse.center(item);
    const horz = await item.evaluate((e) => e.closest('.d4-menu-item')?.classList.contains('d4-menu-item-horz') ?? false);
    if (horz) {
      const bottom = await item.evaluate((e) => e.closest('.d4-menu-item').getBoundingClientRect().bottom);
      // no pause below the bar: it is often a grid header, whose tooltip would pop up
      await mouse.moveTo(x, bottom + 12);
      await mouse.press(x, y, {after: 650});
    }
    else if (i < labels.length - 1) {
      // move horizontally first so the pointer doesn't cross sibling groups diagonally
      if (i === 0) { await mouse.moveTo(mouse.x, y); await mouse.moveTo(x, y); }
      else if (prevHorz) { await mouse.moveTo(mouse.x, y); await mouse.moveTo(x, y); }
      else { await mouse.moveTo(x, mouse.y); await mouse.moveTo(x, y); }
      await mouse.page.waitForTimeout(650);
    }
    else {
      if (i === 0 || prevHorz) await mouse.moveTo(mouse.x, y);
      else await mouse.moveTo(x, mouse.y);
      await mouse.click(x, y, {after: 900});
    }
    prevHorz = horz;
  }
}

/** Runs setup code (not recorded), then records scene(page, mouse) into <out>.gif. The browser is closed on failure too. */
export async function film(name, setup, scene, {out, thumbAt, fps, pre, start, noTooltips, crop, size, colors} = {}) {
  const {browser, page} = await open();
  try {
    await cleanShell(page);
    await page.evaluate(`(async () => { ${setup} })()`);
    await page.waitForTimeout(2500);
    if (pre) await pre(page);
    if (noTooltips) await page.addStyleTag({content: '.d4-tooltip { display: none !important; }'});
    await installPointer(page);
    const mouse = new Mouse(page);
    if (start) { mouse.x = start[0]; mouse.y = start[1]; }
    await page.mouse.move(mouse.x, mouse.y);
    const rec = new Recorder(page, path.join(os.tmpdir(), 'rec-lib-frames', name));
    await rec.start([mouse.x + 1, mouse.y + 1]);
    try { await scene(page, mouse); }
    finally { await rec.stop(); }
    const cropBox = typeof crop === 'function' ? await crop(page) : crop;
    // Default output: ./out in the current directory. Keep it outside the repository.
    const target = out ?? path.join(process.cwd(), 'out', name);
    fs.mkdirSync(path.dirname(target), {recursive: true});
    rec.make(target, {thumbAt, fps, crop: cropBox, size, colors});
    rec.sheet(`${target}-sheet.png`);
    fs.rmSync(rec.dir, {recursive: true, force: true});
  }
  finally { await browser.close(); }
}
