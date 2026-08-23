/* ---
realizes: [chem.cp.molecule-cell-actions]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

declare const grok: any;
declare const DG: any;

const DATASET = 'System:DemoFiles/chem/smiles.csv';
const MOL_COL = 'canonical_smiles';

async function openChemTable(page: Page, path: string): Promise<void> {
  await page.evaluate(async (p) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(p);
    grok.shell.addTableView(df);
    (window as any).__df = df;
    await new Promise<void>((res) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
      setTimeout(res, 8000);
    });
  }, path);
  await waitForChemMenu(page);
  await waitForMolecule(page);
  for (let i = 0; i < 50; i++) {
    if (await page.locator('[name="viewer-Grid"] canvas').count() > 0) break;
    await page.waitForTimeout(200);
  }
  await page.locator('[name="viewer-Grid"] canvas').first().waitFor({timeout: 30_000, state: 'attached'});
  await page.waitForTimeout(4000);
}

async function molCellPoint(page: Page, visualRow: number): Promise<{x: number; y: number}> {
  return page.evaluate(({col, vr}) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(col, vr).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {col: MOL_COL, vr: visualRow});
}

async function openCellMenu(page: Page, visualRow: number): Promise<void> {
  await page.keyboard.press('Escape').catch(() => {});
  await page.mouse.click(2, 2).catch(() => {});
  await page.waitForTimeout(200);
  const p = await molCellPoint(page, visualRow);
  const overlay = page.locator('[name="viewer-Grid"] canvas[name="overlay"]');
  if (await overlay.count() > 0) await overlay.first().evaluate((el: HTMLElement) => el.focus());
  await page.mouse.click(p.x, p.y);
  await page.evaluate(async ({vr, col}) => {
    const df = grok.shell.tv.dataFrame;
    const grid = grok.shell.tv.grid;
    const tr = grid.gridRowToTable(vr);
    for (let i = 0; i < 40; i++) {
      if (df.currentRowIdx === tr && df.currentCol != null && df.currentCol.name === col) break;
      await new Promise((r) => setTimeout(r, 50));
    }
  }, {vr: visualRow, col: MOL_COL});
  await page.mouse.move(p.x, p.y);
  await page.mouse.click(p.x, p.y, {button: 'right'});
  await page.locator('.d4-menu-popup').waitFor({timeout: 5_000});
  await page.waitForTimeout(300);
}

async function menuItemBoxes(page: Page): Promise<{label: string; w: number; h: number}[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup'))
    .flatMap((p) => Array.from(p.querySelectorAll('.d4-menu-item')))
    .map((i) => {
      const label = (i.querySelector('.d4-menu-item-label')?.textContent ?? '').trim();
      const r = i.getBoundingClientRect();
      return {label, w: Math.round(r.width), h: Math.round(r.height)};
    }).filter((i) => i.label.length > 0));
}

// Selector recon-notes — see chem.md "Molecule-cell context menu — the Copy-as leaves are a submenu":
// - the Copy-as / Export / Sort leaves sit at 0x0 until the "Current Value" group is hovered open,
//   and the expanded group renders them in a SECOND .d4-menu-popup
// - Playwright text/actionability resolution does not reach these items; click the item's own rect
async function expandCurrentValueGroup(page: Page): Promise<void> {
  const groupPoint = await page.evaluate(() => {
    const popups = Array.from(document.querySelectorAll('.d4-menu-popup'));
    let best: {el: Element; r: DOMRect} | null = null;
    for (const p of popups) {
      for (const el of Array.from(p.querySelectorAll('*'))) {
        if (!/^current\s+value$/i.test((el.textContent ?? '').trim())) continue;
        const r = el.getBoundingClientRect();
        if (r.width === 0 || r.height === 0) continue;
        // technical: deepest match wins — ancestors carry the same text
        if (!best || el.compareDocumentPosition(best.el) & Node.DOCUMENT_POSITION_PRECEDING) best = {el, r};
      }
    }
    if (!best) return null;
    return {
      x: Math.round(best.r.x + best.r.width / 2),
      y: Math.round(best.r.y + best.r.height / 2),
      w: Math.round(best.r.width),
      h: Math.round(best.r.height),
      cls: best.el.className?.toString?.() ?? '',
      tag: best.el.tagName,
    };
  });
  console.log(`[cell-actions] Current Value group = ${JSON.stringify(groupPoint)}`);
  const laidOutLabels = await page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup'))
    .flatMap((p) => Array.from(p.querySelectorAll('.d4-menu-item')))
    .map((i) => {
      const t = (i.querySelector('.d4-menu-item-label')?.textContent ?? '').trim();
      const r = i.getBoundingClientRect();
      return t ? `${t}${r.width > 0 && r.height > 0 ? '' : ' [0x0]'}` : '';
    }).filter(Boolean).slice(0, 60));
  expect(groupPoint,
    'the cell context menu must expose a laid-out "Current Value" group — grid_popup_menu.dart#L75-85 ' +
    'registers the SemanticValue actions under it. Laid-out labels actually present: ' +
    JSON.stringify(laidOutLabels)).not.toBeNull();
  // technical: MenuItem._lastMouseMove is static and dedupes on the client point (menu.dart#L79,
  // #L207-208) — a single move to the same coordinates on a reopened menu is dropped before
  // _showSubMenu (#L566), so land on a second point inside the group first. See chem.md.
  const nudge = Math.max(1, Math.min(8, Math.floor(groupPoint!.w / 2) - 2));
  expect(nudge, `the pre-move offset must land inside the "Current Value" item, whose rect is ` +
    `${groupPoint!.w}x${groupPoint!.h}`).toBeLessThan(groupPoint!.w / 2);
  await page.mouse.move(groupPoint!.x - nudge, groupPoint!.y);
  await page.mouse.move(groupPoint!.x, groupPoint!.y);
  // technical: barrier on the group opening — a Copy-as leaf acquiring layout
  await page.waitForFunction(() => Array.from(document.querySelectorAll('.d4-menu-item-label'))
    .some((l) => {
      const t = (l.textContent ?? '').trim();
      if (!t.startsWith('Copy as')) return false;
      const r = (l.closest('.d4-menu-item') ?? l).getBoundingClientRect();
      return r.width > 0 && r.height > 0;
    }), null, {timeout: 10_000});
}

async function clickCellMenuItem(page: Page, visualRow: number, label: string): Promise<void> {
  await openCellMenu(page, visualRow);
  await expandCurrentValueGroup(page);
  const present = await page.evaluate((lbl) => Array.from(document.querySelectorAll('.d4-menu-popup'))
    .some((p) => Array.from(p.querySelectorAll('.d4-menu-item-label'))
      .some((l) => (l.textContent ?? '').trim() === lbl)), label);
  expect(present, `context-menu item "${label}" not found on the molecule cell`).toBe(true);
  const box = await page.evaluate((lbl) => {
    const popups = Array.from(document.querySelectorAll('.d4-menu-popup'));
    const laidOut = popups.filter((p) => {
      const r = p.getBoundingClientRect(); return r.width > 0 && r.height > 0;
    });
    const popup = laidOut.reverse().find((p) => Array.from(p.querySelectorAll('.d4-menu-item-label'))
      .some((l) => (l.textContent ?? '').trim() === lbl)) ?? laidOut[0];
    const diag = {
      popups: popups.length,
      laidOut: laidOut.length + (popup ? 1 : 0),
      popupRects: popups.map((p) => { const r = p.getBoundingClientRect(); return [Math.round(r.width), Math.round(r.height)]; }),
    };
    if (!popup) return {box: null, diag: {...diag, why: 'no laid-out popup'}};
    const item = Array.from(popup.querySelectorAll('.d4-menu-item-label'))
      .find((l) => (l.textContent ?? '').trim() === lbl);
    if (!item) return {box: null, diag: {...diag, why: 'label not in laid-out popup'}};
    const host = (item.closest('.d4-menu-item') ?? item) as HTMLElement;
    const r = host.getBoundingClientRect();
    const cs = getComputedStyle(host);
    const pcs = getComputedStyle(popup as HTMLElement);
    return {
      box: {x: Math.round(r.x + r.width / 2), y: Math.round(r.y + r.height / 2), w: r.width, h: r.height},
      diag: {
        ...diag,
        itemRect: [Math.round(r.width), Math.round(r.height)],
        itemStyle: {display: cs.display, visibility: cs.visibility, opacity: cs.opacity},
        popupStyle: {display: pcs.display, visibility: pcs.visibility, opacity: pcs.opacity, overflow: pcs.overflow},
        itemsInPopup: popup.querySelectorAll('.d4-menu-item').length,
        visibleLabels: Array.from(popup.querySelectorAll('.d4-menu-item'))
          .filter((i) => { const b = i.getBoundingClientRect(); return b.width > 0 && b.height > 0; })
          .map((i) => (i.querySelector('.d4-menu-item-label')?.textContent ?? '').trim())
          .filter(Boolean),
      },
    };
  }, label);
  console.log(`[cell-actions] menu diag for "${label}" = ${JSON.stringify(box.diag)}`);
  expect(box.box, `context-menu item "${label}" not reachable — ${JSON.stringify(box.diag)}`).not.toBeNull();
  expect(box.box!.w * box.box!.h,
    `context-menu item "${label}" measured zero area — ${JSON.stringify(box.diag)}`)
    .toBeGreaterThan(0);
  await page.mouse.move(box.box!.x, box.box!.y);
  await page.mouse.click(box.box!.x, box.box!.y);
}

async function clearClipboard(page: Page): Promise<void> {
  const left = await page.evaluate(async () => {
    try {
      await navigator.clipboard.writeText('');
    } catch (e) {
      return 'write rejected: ' + String(e);
    }
    return await navigator.clipboard.readText().catch((e) => 'read rejected: ' + String(e));
  });
  console.log(`[cell-actions] clipboard cleared, reads back ${JSON.stringify(left)}`);
  expect(left,
    'the clipboard must read back empty before the copy action — otherwise a stale value from the ' +
    'previous action can satisfy the format checks while this action never fires').toBe('');
}

async function readClipboardText(page: Page): Promise<string> {
  return page.evaluate(async () => {
    for (let i = 0; i < 20; i++) {
      try {
        const t = await navigator.clipboard.readText();
        if (t && t.length > 0) return t;
      } catch (e) {}
      await new Promise((r) => setTimeout(r, 100));
    }
    return await navigator.clipboard.readText().catch(() => '');
  });
}

async function readClipboardImageBytes(page: Page): Promise<number> {
  return page.evaluate(async () => {
    for (let i = 0; i < 25; i++) {
      try {
        const items = await navigator.clipboard.read();
        for (const it of items) {
          const pngType = it.types.find((t) => t === 'image/png');
          if (pngType) {
            const blob = await it.getType(pngType);
            if (blob.size > 0) return blob.size;
          }
        }
      } catch (e) {}
      await new Promise((r) => setTimeout(r, 120));
    }
    return -1;
  });
}

// technical: installed once, before any action — balloons auto-hide, see chem.md
async function installBalloonRecorder(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__balloonLog) return;
    w.__balloonLog = [];
    w.__consoleErrLog = [];
    const origErr = console.error;
    console.error = function (...args: any[]) {
      w.__consoleErrLog.push({t: Date.now(), text: args.map((a: any) => String(a)).join(' ')});
      origErr.apply(console, args as any);
    };
    new MutationObserver((muts: any) => {
      for (const m of muts) for (const n of m.addedNodes) {
        if (n.nodeType !== 1) continue;
        const el = n as Element;
        const hit = (el.classList?.contains('d4-balloon') ? el : el.querySelector?.('.d4-balloon')) as Element | null;
        if (hit) w.__balloonLog.push({t: Date.now(), text: (hit.textContent || '').trim()});
      }
    }).observe(document.body, {childList: true, subtree: true});
  });
}

async function readActionTrace(page: Page, sinceMs: number): Promise<{balloons: any[]; errors: any[]}> {
  return page.evaluate((since) => {
    const w = window as any;
    return {
      balloons: (w.__balloonLog ?? []).filter((b: any) => b.t >= since),
      errors: (w.__consoleErrLog ?? []).filter((e: any) => e.t >= since),
    };
  }, sinceMs);
}

async function canonicalize(page: Page, smiles: string): Promise<string> {
  return page.evaluate((smi) => {
    const mb = grok.chem.convert(smi, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    return grok.chem.convert(mb, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
  }, smiles);
}

test.use(specTestOptions);

test('Chem: Molecule cell Copy-as / Export / Sort-by-similarity actions', async ({page, context}) => {
  test.setTimeout(300_000);

  await context.grantPermissions(['clipboard-read', 'clipboard-write']);
  await loginToDatagrok(page);
  await installBalloonRecorder(page);

  const clip: Record<string, string> = {};
  const queryRow = 2;

  await softStep('Setup: open smiles.csv, canonical_smiles typed as Molecule', async () => {
    await openChemTable(page, DATASET);
    const meta = await page.evaluate((col) => {
      const df = (window as any).__df;
      const c = df.columns.byName(col);
      return {semType: c.semType, rows: df.rowCount, val0: df.get(col, 0)};
    }, MOL_COL);
    expect(meta.semType, 'canonical_smiles must detect as Molecule for cell actions to bind').toBe('Molecule');
    expect(meta.rows).toBeGreaterThan(5);
    (page as any).__cell0Smiles = meta.val0;
  });

  await softStep('Scenario 1 Step 1: cell context menu lists all seven molecule-cell actions', async () => {
    await openCellMenu(page, 0);
    await expandCurrentValueGroup(page);
    const items = await menuItemBoxes(page);
    console.log(`[cell-actions] cell menu items = ${JSON.stringify(items)}`);
    for (const expected of ['Copy as SMILES', 'Copy as MOLFILE V2000', 'Copy as MOLFILE V3000',
      'Copy as SMARTS', 'Copy as Image', 'Export as SVG', 'Sort by similarity']) {
      const hit = items.find((i) => i.label === expected);
      expect(hit, `cell context menu missing "${expected}" — saw ${JSON.stringify(items.map((i) => i.label))}`)
        .not.toBeUndefined();
      expect(hit!.w * hit!.h,
        `context-menu item "${expected}" is in the DOM but measured zero area — a menu that shows ` +
        `nothing lists nothing`).toBeGreaterThan(0);
    }
    await page.keyboard.press('Escape');
  });

  await softStep('Scenario 1 Step 4-5: Copy as SMILES — balloon + valid round-tripping SMILES on clipboard', async () => {
    await clearClipboard(page);
    const t0 = await page.evaluate(() => Date.now());
    await clickCellMenuItem(page, 0, 'Copy as SMILES');
    const trace = await readActionTrace(page, t0);
    console.log(`[cell-actions] Copy as SMILES trace = ${JSON.stringify(trace)}`);
    const balloons = trace.balloons.map((b: any) => b.text);
    expect(trace.errors.length,
      `Copy as SMILES raised console errors — grok.shell.info runs unconditionally after ` +
      `writeText, so an error before it means the handler threw: ${JSON.stringify(trace.errors)}`)
      .toBe(0);
    expect(balloons.some((t: string) => /smiles copied to clipboard/i.test(t)),
      `expected "Smiles copied to clipboard" balloon, recorder saw ${JSON.stringify(balloons)}`).toBe(true);
    const smiles = await readClipboardText(page);
    clip.smiles = smiles;
    expect(smiles.length, 'clipboard SMILES must be non-empty').toBeGreaterThan(0);
    expect(/\n/.test(smiles), 'a SMILES string carries no embedded newline').toBe(false);
    const cellSmiles = (page as any).__cell0Smiles as string;
    const [canonClip, canonCell] = await Promise.all([canonicalize(page, smiles), canonicalize(page, cellSmiles)]);
    console.log(`[cell-actions] canonical clipboard = ${JSON.stringify(canonClip)}, ` +
      `canonical cell = ${JSON.stringify(canonCell)}`);
    // technical: _convertMolNotation swallows a parse failure and RETURNS a sentinel instead of
    // throwing (convert-notation-utils.ts#L60-62 picks 'MALFORMED_INPUT_VALUE' for a SMILES target
    // and MALFORMED_MOL_V2000 / MALFORMED_MOL_V3000 for a molblock target; #L97-100 returns it from
    // the finally block) — with RDKit degraded both sides collapse to the same sentinel and the
    // equality below would hold no matter what Copy as SMILES put on the clipboard
    for (const [what, canon] of [['clipboard', canonClip], ['cell', canonCell]] as [string, string][]) {
      expect(canon.length,
        `canonicalizing the ${what} SMILES yielded an empty string — the round-trip never ran`)
        .toBeGreaterThan(0);
      expect(canon === 'MALFORMED_INPUT_VALUE' || /\bMalformed\b/.test(canon),
        `canonicalizing the ${what} SMILES returned an RDKit failure sentinel ` +
        `(${JSON.stringify(canon)}), so the equality below would compare sentinel to sentinel`)
        .toBe(false);
    }
    expect(canonClip, 'clipboard SMILES must round-trip to the same molecule as the cell').toBe(canonCell);
  });

  await softStep('Scenario 1 Step 6: Copy as MOLFILE V2000 — blank-header V2000 molblock ending M  END', async () => {
    await clearClipboard(page);
    const t0 = await page.evaluate(() => Date.now());
    await clickCellMenuItem(page, 0, 'Copy as MOLFILE V2000');
    console.log(`[cell-actions] Copy as MOLFILE V2000 trace = ${JSON.stringify(await readActionTrace(page, t0))}`);
    const mol = await readClipboardText(page);
    clip.v2000 = mol;
    expect(mol.length, 'clipboard V2000 molfile must be non-empty').toBeGreaterThan(0);
    expect(/^\s*\r?\n/.test(mol) || mol.startsWith('\n'), 'a V2000 molfile begins with a blank header line').toBe(true);
    expect(mol.includes('V2000'), 'V2000 molfile must carry the V2000 version tag').toBe(true);
    expect(mol.includes('M  END'), 'a molfile ends with the M  END terminator').toBe(true);
    expect(mol.includes('V3000'), 'a V2000 molfile must not carry the V3000 tag').toBe(false);
  });

  await softStep('Scenario 1 Step 7: Copy as MOLFILE V3000 — V3000 molblock distinct from the V2000 text', async () => {
    await clearClipboard(page);
    const t0 = await page.evaluate(() => Date.now());
    await clickCellMenuItem(page, 0, 'Copy as MOLFILE V3000');
    console.log(`[cell-actions] Copy as MOLFILE V3000 trace = ${JSON.stringify(await readActionTrace(page, t0))}`);
    const mol = await readClipboardText(page);
    clip.v3000 = mol;
    expect(mol.length, 'clipboard V3000 molfile must be non-empty').toBeGreaterThan(0);
    expect(mol.includes('V3000'), 'V3000 molfile must carry the V3000 version tag').toBe(true);
    expect(mol.includes('M  END'), 'a molfile ends with the M  END terminator').toBe(true);
    expect(mol, 'the V3000 molfile text must differ from the V2000 molfile text').not.toBe(clip.v2000);
  });

  await softStep('Scenario 1 Step 8: Copy as SMARTS — SMARTS token present, distinct from the copied SMILES', async () => {
    await clearClipboard(page);
    const t0 = await page.evaluate(() => Date.now());
    await clickCellMenuItem(page, 0, 'Copy as SMARTS');
    console.log(`[cell-actions] Copy as SMARTS trace = ${JSON.stringify(await readActionTrace(page, t0))}`);
    const smarts = await readClipboardText(page);
    clip.smarts = smarts;
    expect(smarts.length, 'clipboard SMARTS must be non-empty').toBeGreaterThan(0);
    // technical: get_smarts() always emits bracket atoms ([#6], [#7], ...); a plain SMILES and a
    // molfile both carry '-'/'=' /'#', so only the bracket-atom form discriminates
    expect(/\[#\d/.test(smarts),
      `SMARTS text must contain a bracket-atom token like "[#6]" — got ${JSON.stringify(smarts.slice(0, 120))}`)
      .toBe(true);
    expect(smarts, 'the SMARTS text must differ character-for-character from the copied SMILES').not.toBe(clip.smiles);
  });

  await softStep('Scenario 1 Step 9: Copy as Image — balloon + non-trivial PNG blob on clipboard', async () => {
    await clearClipboard(page);
    const t0 = await page.evaluate(() => Date.now());
    await clickCellMenuItem(page, 0, 'Copy as Image');
    await page.waitForFunction((since) => ((window as any).__balloonLog ?? [])
      .some((b: any) => b.t >= since && /image copied to clipboard/i.test(b.text)),
    t0, {timeout: 10_000}).catch(() => {});
    const balloons = (await readActionTrace(page, t0)).balloons.map((b: any) => b.text);
    console.log(`[cell-actions] Copy as Image balloons = ${JSON.stringify(balloons)}`);
    expect(balloons.some((t: string) => /image copied to clipboard/i.test(t)),
      `expected "Image copied to clipboard" balloon, saw ${JSON.stringify(balloons)}`).toBe(true);
    const bytes = await readClipboardImageBytes(page);
    expect(bytes, 'no image/png blob found on the clipboard after Copy as Image').toBeGreaterThan(0);
    expect(bytes, 'the copied PNG must be a non-trivial render (> 500 bytes)').toBeGreaterThan(500);
  });

  await softStep('Scenario 1 Step 10: Export as SVG — molecule.svg downloads with drawing content, no failure warning', async () => {
    // technical: kept handled — if the action below throws, softStep catches it and this promise
    // would otherwise reject unhandled 15 s later and tear down the browser context mid-run
    const downloadPromise = page.waitForEvent('download', {timeout: 15_000});
    downloadPromise.catch(() => {});
    const t0 = await page.evaluate(() => Date.now());
    await clickCellMenuItem(page, 0, 'Export as SVG');
    const svgBalloons = (await readActionTrace(page, t0)).balloons.map((b: any) => b.text);
    console.log(`[cell-actions] Export as SVG balloons = ${JSON.stringify(svgBalloons)}`);
    expect(svgBalloons.filter((t: string) => /failed to export structure as svg/i.test(t)),
      'exportAsSvg raises the "Failed to export structure as SVG" warning balloon and returns without ' +
      `downloading when the renderer yields nothing; the recorder saw ${JSON.stringify(svgBalloons)}`).toEqual([]);
    const download = await downloadPromise;
    expect(download.suggestedFilename(), 'the exported file must be named molecule.svg').toBe('molecule.svg');
    const stream = await download.createReadStream();
    const chunks: Buffer[] = [];
    for await (const c of stream) chunks.push(c as Buffer);
    const svg = Buffer.concat(chunks).toString('utf8');
    console.log(`[cell-actions] Export as SVG head = ${JSON.stringify(svg.slice(0, 120))}, bytes = ${svg.length}`);
    expect(/^\s*(<\?xml[^>]*\?>\s*)?(<!DOCTYPE[^>]*>\s*)?<svg[\s>]/.test(svg),
      `the downloaded file must be an SVG document; head = ${JSON.stringify(svg.slice(0, 120))}`).toBe(true);
    expect(/<(path|line|circle)\b/.test(svg), 'the SVG must contain at least one drawing primitive (path/line/circle)').toBe(true);
  });

  await softStep('Scenario 2 Step 3 (partial): the four clipboard formats are pairwise distinct', async () => {
    const texts = [clip.smiles, clip.v2000, clip.v3000, clip.smarts];
    console.log(`[cell-actions] clipboard lengths = ${JSON.stringify(texts.map((t) => t?.length ?? 0))}`);
    for (const t of texts) expect(t?.length ?? 0, 'each format must have produced clipboard text').toBeGreaterThan(0);
    const uniq = new Set(texts);
    expect(uniq.size, 'SMILES / V2000 / V3000 / SMARTS must be four pairwise-distinct strings').toBe(4);
  });

  await softStep('Scenario 2 Step 1-2: Sort by similarity places the query molecule first and pins it', async () => {
    await openChemTable(page, DATASET);
    const query = await page.evaluate(({r, col}) => (window as any).__df.get(col, r), {r: queryRow, col: MOL_COL});
    const orderBefore = await page.evaluate(() => grok.shell.tv.grid.getRowOrder().slice(0, 6));
    expect(orderBefore[0], 'before sorting the grid is in natural order (row 0 = table row 0)').toBe(0);

    await clickCellMenuItem(page, queryRow, 'Sort by similarity');

    // technical: barrier on the last effect sortBySimilarity writes (package.ts#L2129-2130)
    await page.waitForFunction(() => {
      const pinned = grok.shell.tv.grid.props.pinnedRowValues;
      return pinned != null && pinned.length > 0;
    }, null, {timeout: 60_000});

    const after = await page.evaluate((col) => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const order = grid.getRowOrder();
      const fpCol = df.col('~' + col + '.Morgan');
      return {
        firstTableIdx: order[0],
        firstMol: df.get(col, order[0]),
        pinnedValues: grid.props.pinnedRowValues,
        pinnedCols: grid.props.pinnedRowColumnNames,
        columns: df.columns.names(),
        fpExport: fpCol == null ? null :
          {csv: fpCol.meta.includeInCsvExport, binary: fpCol.meta.includeInBinaryExport},
      };
    }, MOL_COL);
    console.log(`[cell-actions] sort-by-similarity probe = ${JSON.stringify(after)}`);
    expect(after.firstTableIdx, 'the right-clicked (query) table row must lead the sorted grid').toBe(queryRow);
    expect(after.firstMol, 'row 0 of the sorted grid must show the query molecule').toBe(query);
    expect(after.pinnedValues, 'the query molecule must be pinned as the sort anchor').toContain(query);
    expect(after.pinnedCols, 'the pinned column must be the molecule column').toContain(MOL_COL);
    expect(after.columns,
      `Sort by similarity fingerprints the molecule column, and saveColumn (chem-searches.ts#L176-187) ` +
      `caches the Morgan fingerprints on the table as the working column "~${MOL_COL}.Morgan"`)
      .toContain(`~${MOL_COL}.Morgan`);
    expect(after.fpExport,
      `the "~${MOL_COL}.Morgan" cache column must be readable back to check its export flags`)
      .not.toBeNull();
    expect(after.fpExport!.csv,
      `saveColumn sets includeInCsvExport = false on the fingerprint cache column ` +
      `(chem-searches.ts#L185) — without it the cache leaks into every CSV export`).toBe(false);
    expect(after.fpExport!.binary,
      `saveColumn sets includeInBinaryExport = false on the fingerprint cache column ` +
      `(chem-searches.ts#L186) — without it the cache leaks into every binary export`).toBe(false);
    expect(after.columns,
      'the similarity SCORE is not attached to the table — similarityResultDf builds smiles/score/indexes in a ' +
      'separate DataFrame (chem-similarity-viewer.ts#L311-334) and sortBySimilarity consumes only "indexes"')
      .not.toContain('score');
  });

  await softStep('Scenario 2 Step 3 (scores): similarity does not increase down the sorted grid, and the query scores 1.0', async () => {
    const probe = await page.evaluate(async ({col, qr}) => {
      const df = grok.shell.tv.dataFrame;
      // technical: getRowOrder() is an Int32Array — mapping it in place would return another
      // Int32Array and truncate every similarity to 0, making the monotonicity check unfalsifiable
      const order: number[] = Array.from(grok.shell.tv.grid.getRowOrder().slice(0, 5));
      const scoreCol = await grok.chem.getSimilarities(df.columns.byName(col), df.get(col, qr));
      return {order, scores: order.map((i: number) => scoreCol.get(i))};
    }, {col: MOL_COL, qr: queryRow});
    console.log(`[cell-actions] sorted rows 0-4 = ${JSON.stringify(probe.order)}, ` +
      `Tanimoto/Morgan similarity to the query = ${JSON.stringify(probe.scores)}`);
    expect(probe.scores.length, 'the sorted grid must expose five leading rows to score').toBe(5);
    for (const s of probe.scores)
      expect(typeof s === 'number' && Number.isFinite(s),
        `every scored row needs a real similarity, got ${JSON.stringify(probe.scores)}`).toBe(true);
    expect(probe.scores[0], 'the query molecule scores 1.0 against itself at row 0').toBeCloseTo(1.0, 6);
    for (let i = 1; i < probe.scores.length; i++)
      expect(probe.scores[i],
        `similarity must not rise from sorted row ${i - 1} to row ${i} — ${JSON.stringify(probe.scores)}`)
        .toBeLessThanOrEqual(probe.scores[i - 1]);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
