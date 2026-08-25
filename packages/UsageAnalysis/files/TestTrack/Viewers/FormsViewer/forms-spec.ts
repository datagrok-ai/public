/* ---
realizes: [formsviewer.cp.color-and-renderer-presentation, formsviewer.edge.molecule-column-render, formsviewer.edge.fit-semtype-promotes-renderer-size]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {HOST, ORDINARY, CURRENT, withConsoleErrorCount} from '../../helpers/forms';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const curvesPath = 'System:DemoFiles/curves.csv';

function toRgb(s: string | null): [number, number, number] | null {
  if (!s) return null;
  let m = s.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)/i);
  if (m) return [+m[1], +m[2], +m[3]];
  m = s.match(/^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i);
  return m ? [parseInt(m[1], 16), parseInt(m[2], 16), parseInt(m[3], 16)] : null;
}

async function colorCheck(page: Page, cardSel: string, col: string):
  Promise<{bg: string; ref: string | null; colorInt: number} | null> {
  return page.evaluate(({sel, c}) => {
    const df = grok.shell.t;
    const card = document.querySelector(sel);
    const el = card ? card.querySelector(`[column="${c}"]`) as HTMLElement | null : null;
    if (!el) return null;
    const row = df.currentRowIdx;
    const colorInt = df.col(c).meta.colors.getColor(row);
    const ref = (colorInt && colorInt !== 4294967295) ? DG.Color.toHtml(colorInt) : null;
    return {bg: getComputedStyle(el).backgroundColor, ref, colorInt};
  }, {sel: cardSel, c: col});
}

async function styleCheck(page: Page, cardSel: string, col: string):
  Promise<{fieldAlign: string; refAlign: string | null; fieldFont: string; expectFont: string} | null> {
  return page.evaluate(({sel, c}) => {
    const grid = grok.shell.tv.grid;
    const gc = grid.col(c);
    const card = document.querySelector(sel);
    const el = card ? card.querySelector(`[column="${c}"]`) as HTMLElement | null : null;
    if (!el || !gc) return null;
    const cs = gc.contentCellStyle;
    const probe = document.createElement('div');
    probe.style.font = cs && cs.font ? cs.font : '';
    return {
      fieldAlign: getComputedStyle(el).textAlign,
      refAlign: cs ? (cs.horzAlign ?? null) : null,
      fieldFont: el.style.font,
      expectFont: probe.style.font,
    };
  }, {sel: cardSel, c: col});
}

async function sizeCheck(page: Page, col: string):
  Promise<{tag: string; aw: number | null; ah: number | null; ew: number; eh: number; size: string} | null> {
  return page.evaluate(({sel, c}) => {
    const grid = grok.shell.tv.grid;
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
    const card = document.querySelector(sel);
    const el = card ? card.querySelector(`[column="${c}"]`) as any : null;
    if (!el || !vw) return null;
    const dpr = window.devicePixelRatio;
    const size = vw.props.rendererSize as string;
    const renderer = grid.col(c) ? grid.col(c).renderer : null;
    const w = renderer ? renderer.defaultWidth : null;
    const h = renderer ? renderer.defaultHeight : null;
    let cw: number; let ch: number;
    if (!w || !h) {
      if (size === 'normal') { cw = 200; ch = 100; }
      else if (size === 'large') { cw = 300; ch = 150; }
      else { cw = 120; ch = 60; }
    } else {
      if (size === 'normal') { cw = w; ch = h; }
      else if (size === 'large') { cw = Math.floor(w * 1.5); ch = Math.floor(h * 1.5); }
      else { cw = Math.floor(w * 0.66); ch = Math.floor(h * 0.66); }
    }
    return {tag: el.tagName, aw: el.width ?? null, ah: el.height ?? null,
      ew: Math.floor(cw * dpr), eh: Math.floor(ch * dpr), size};
  }, {sel: CURRENT, c: col});
}

const rendererSizeProp = (page: Page): Promise<string | null> => page.evaluate(() => {
  const vw = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
  return (vw?.props?.rendererSize as string) ?? null;
});

async function setRendererSizeViaPanel(page: Page, value: 'small' | 'normal' | 'large'): Promise<void> {
  // re-open the gear for the CURRENT Forms viewer: a new openTable/addViewer leaves the
  // property grid bound to the previous view's viewer, so edits never reach this one
  await v.clickViewerTitlebarIcon(page, 'Forms', 'icon-font-icon-settings').catch(() => {});
  await page.waitForTimeout(600);
  await v.ensurePropertyCategory(page, 'Forms', 'misc', 'renderer-size');
  // the row can sit in a COLLAPSED category: display stays 'table-row' but the box is 0x0,
  // so isVisible()-style gating passes while clicks and selectOption reach nothing. Expand
  // every header until the row has a real box before touching the editor.
  const row = page.locator('.property-grid tr[name="prop-renderer-size"]').first();
  const sized = () => page.evaluate(() => {
    const r = document.querySelector('.property-grid tr[name="prop-renderer-size"]') as HTMLElement | null;
    if (!r) return false;
    const b = r.getBoundingClientRect();
    return b.width > 0 && b.height > 0;
  });
  if (!await sized()) {
    for (const h of await page.locator('[name^="prop-category-"]').all())
      if (await h.isVisible().catch(() => false)) await h.click().catch(() => {});
    await v.pollValue(sized, (ok) => ok, 3000, 100);
  }

  await row.locator('td').last().click();
  const sel = page.locator('[name="prop-renderer-size"] select').first();
  await sel.waitFor({state: 'visible', timeout: 5000});
  await sel.selectOption(value);
  await v.pollValue(() => rendererSizeProp(page), (s) => s === value, 3000, 100);
  expect(await rendererSizeProp(page), `rendererSize did not commit to "${value}"`).toBe(value);
}

async function setColorCodeViaPanel(page: Page, on: boolean): Promise<void> {
  await v.ensurePropertyCategory(page, 'Forms', 'misc', 'color-code');
  await v.setPropertyGridCheckbox(page, 'color-code', on, 'misc');
}

test('Forms viewer — colour coding and renderer presentation (p2)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  let colorRow = 0;
  await softStep('Scenario 1 setup — colour-code AGE by background, centre it, give it a distinct font', async () => {
    await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});

    const info = await page.evaluate(() => {
      const df = grok.shell.t;
      const gc = grok.shell.tv.grid.col('AGE');
      gc.contentCellStyle.horzAlign = 'center';
      gc.contentCellStyle.font = 'italic bold 14px "Times New Roman"';
      df.col('AGE').meta.colors.setLinear();

      const col = df.col('AGE');
      let row = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const c = col.meta.colors.getColor(i);
        if (c && c !== 4294967295) { row = i; break; }
      }
      df.currentRowIdx = row;
      return {row, textCoded: !!gc.isTextColorCoded};
    });
    colorRow = info.row;
    expect(info.textCoded).toBe(false);
  });

  await softStep('Step 1a — the AGE field background equals the AGE colour scheme, both normalized', async () => {

    await expect.poll(async () => {
      const r = await colorCheck(page, CURRENT, 'AGE');
      if (!r || !r.ref) return false;
      const a = toRgb(r.bg); const b = toRgb(r.ref);
      return !!a && !!b && a[0] === b[0] && a[1] === b[1] && a[2] === b[2];
    }, {timeout: 20_000}).toBe(true);
    const r = await colorCheck(page, CURRENT, 'AGE');
    expect(r).not.toBeNull();
    expect(r!.colorInt).not.toBe(4294967295);
    expect(r!.ref).not.toBeNull();
  });

  await softStep('Step 1b — the AGE field horizontal alignment and font equal the grid column style', async () => {
    const s = await styleCheck(page, CURRENT, 'AGE');
    expect(s).not.toBeNull();

    expect(s!.refAlign).toBe('center');
    expect(s!.fieldAlign).toBe(s!.refAlign);

    expect(s!.expectFont.length).toBeGreaterThan(0);
    expect(s!.fieldFont).toBe(s!.expectFont);
  });

  await softStep('Step 1c — Color Code OFF drops the background to an uncoloured field; ON restores it', async () => {
    await setColorCodeViaPanel(page, false);

    await expect.poll(async () => {
      const r = await page.evaluate((sel) => {
        const card = document.querySelector(sel);
        const a = card ? card.querySelector('[column="AGE"]') as HTMLElement | null : null;
        const b = card ? card.querySelector('[column="USUBJID"]') as HTMLElement | null : null;
        if (!a || !b) return null;
        return {a: getComputedStyle(a).backgroundColor, b: getComputedStyle(b).backgroundColor};
      }, CURRENT);
      return !!r && r.a === r.b;
    }, {timeout: 20_000}).toBe(true);

    await setColorCodeViaPanel(page, true);
    await expect.poll(async () => {
      const r = await colorCheck(page, CURRENT, 'AGE');
      if (!r || !r.ref) return false;
      const a = toRgb(r.bg); const b = toRgb(r.ref);
      return !!a && !!b && a[0] === b[0] && a[1] === b[1] && a[2] === b[2];
    }, {timeout: 20_000}).toBe(true);
  });

  let smallW = 0; let smallH = 0; let normalW = 0; let normalH = 0;
  await softStep('Scenario 2 setup — open spgi-100, add Forms, set a current row', async () => {
    await v.openTable(page, {path: spgiPath, semTypeTimeoutMs: 3000, settleMs: 2000});
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});
    await page.evaluate(() => { grok.shell.t.currentRowIdx = 0; });
    await page.locator(CURRENT).first().waitFor({timeout: 30_000});
  });

  await softStep('Step 2a — at the default small, the Structure canvas is base × 0.66 (floored) × dpr', async () => {

    const def = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').props.rendererSize);
    expect(def).toBe('small');

    await expect.poll(async () => {
      const r = await sizeCheck(page, 'Structure');
      return !!r && r.tag === 'CANVAS' && r.aw === r.ew && r.ah === r.eh;
    }, {timeout: 20_000}).toBe(true);
    const r = await sizeCheck(page, 'Structure');
    smallW = r!.aw!; smallH = r!.ah!;
  });

  await softStep('Step 2b — normal enlarges the canvas to base × 1 × dpr, above the small values', async () => {
    await setRendererSizeViaPanel(page, 'normal');
    await expect.poll(async () => {
      const r = await sizeCheck(page, 'Structure');
      return !!r && r.size === 'normal' && r.aw === r.ew && r.ah === r.eh;
    }, {timeout: 20_000}).toBe(true);
    const r = await sizeCheck(page, 'Structure');
    normalW = r!.aw!; normalH = r!.ah!;
    expect(normalW).toBeGreaterThan(smallW);
    expect(normalH).toBeGreaterThan(smallH);
  });

  await softStep('Step 2c — large enlarges the canvas to base × 1.5 (floored) × dpr, above the normal values', async () => {
    await setRendererSizeViaPanel(page, 'large');
    await expect.poll(async () => {
      const r = await sizeCheck(page, 'Structure');
      return !!r && r.size === 'large' && r.aw === r.ew && r.ah === r.eh;
    }, {timeout: 20_000}).toBe(true);
    const r = await sizeCheck(page, 'Structure');
    expect(r!.aw!).toBeGreaterThan(normalW);
    expect(r!.ah!).toBeGreaterThan(normalH);
  });

  await softStep('Step 3a — the Structure field is a CANVAS carrying its column name; no paint error', async () => {
    const errCount = await withConsoleErrorCount(page, async () => {
      await setRendererSizeViaPanel(page, 'normal');
      await page.evaluate(() => { grok.shell.t.currentRowIdx = 1; grok.shell.t.currentRowIdx = 0; });
      await expect.poll(async () => {
        const r = await sizeCheck(page, 'Structure');
        return !!r && r.tag === 'CANVAS';
      }, {timeout: 20_000}).toBe(true);
    });
    const kind = await page.evaluate((sel) => {
      const el = document.querySelector(`${sel} [column="Structure"]`);
      return el ? el.tagName : null;
    }, CURRENT);
    expect(kind).toBe('CANVAS');
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  await softStep('Step 3b — picking Structure, Core, Primary Series Name yields two molecular canvases per card', async () => {

    await page.evaluate(() => {
      const df = grok.shell.t;
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({fieldsColumnNames: ['Structure', 'Core', 'Primary Series Name']});

      df.mouseOverRowIdx = -1;
      df.currentRowIdx = 0;
      df.selection.setAll(false);
      [3, 7].forEach((r) => df.selection.set(r, true));
    });

    await expect.poll(async () => page.evaluate((sel) => {
      const df = grok.shell.t;
      const picked = ['Structure', 'Core', 'Primary Series Name'];
      const expectedMol = picked.filter((n) => df.col(n) && df.col(n).semType === 'Molecule').length;
      const nonEmpty = Array.from(document.querySelectorAll(sel)).filter((c) => c.querySelector('[column]'));
      const results = nonEmpty.map((card) => {
        const canvases = picked.filter((n) => {
          const el = card.querySelector(`[column="${n}"]`);
          return el && el.tagName === 'CANVAS';
        });
        const psn = card.querySelector('[column="Primary Series Name"]');
        return {
          molCanvases: canvases.length,
          namesOk: JSON.stringify(canvases) === JSON.stringify(['Structure', 'Core']),
          psnIsInput: !!psn && psn.tagName === 'INPUT',
        };
      });
      return JSON.stringify({
        ready: true,
        hasNonEmpty: results.length > 0,
        expectedMol,
        allTwoMol: results.every((r) => r.molCanvases === expectedMol),
        allNamesOk: results.every((r) => r.namesOk),
        allPsnInput: results.every((r) => r.psnIsInput),
      });
    }, ORDINARY), {timeout: 20_000}).toBe(JSON.stringify({
      ready: true, hasNonEmpty: true, expectedMol: 2,
      allTwoMol: true, allNamesOk: true, allPsnInput: true,
    }));
  });

  await softStep('Scenario 4 setup — open curves, add Forms without touching Renderer Size, pick smiles + multiple prefit', async () => {

    await v.openTable(page, {path: curvesPath, semTypeTimeoutMs: 3000, settleMs: 3000});
    const hasFit = await page.evaluate(() =>
      grok.shell.t.columns.names().some((n: string) => grok.shell.t.col(n).semType === 'fit'));
    expect(hasFit).toBe(true);

    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});
    await page.evaluate(() => {
      grok.shell.t.currentRowIdx = 0;
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({fieldsColumnNames: ['smiles', 'multiple prefit']});
    });
    await page.locator(CURRENT).first().waitFor({timeout: 30_000});
  });

  await softStep('Step 4a — the curve canvas is at the normal step of the ladder without touching Renderer Size', async () => {

    await expect.poll(async () => page.evaluate((sel) => {
      const grid = grok.shell.tv.grid;
      const card = document.querySelector(sel);
      const el = card ? card.querySelector('[column="multiple prefit"]') as any : null;
      if (!el) return JSON.stringify({ready: false});
      const dpr = window.devicePixelRatio;
      const renderer = grid.col('multiple prefit') ? grid.col('multiple prefit').renderer : null;
      const w = renderer ? renderer.defaultWidth : null;
      const h = renderer ? renderer.defaultHeight : null;
      const css = (size: string): [number, number] => {
        if (!w || !h) return size === 'normal' ? [200, 100] : size === 'large' ? [300, 150] : [120, 60];
        return size === 'normal' ? [w, h] : size === 'large'
          ? [Math.floor(w * 1.5), Math.floor(h * 1.5)] : [Math.floor(w * 0.66), Math.floor(h * 0.66)];
      };
      const [nw, nh] = css('normal'); const [sw, sh] = css('small');
      return JSON.stringify({
        ready: true,
        isCanvas: el.tagName === 'CANVAS',
        atNormal: el.width === Math.floor(nw * dpr) && el.height === Math.floor(nh * dpr),
        differsFromSmall: el.width !== Math.floor(sw * dpr) || el.height !== Math.floor(sh * dpr),
      });
    }, CURRENT), {timeout: 20_000})
      .toBe(JSON.stringify({ready: true, isCanvas: true, atNormal: true, differsFromSmall: true}));
  });

  await softStep('Step 4b — smiles and multiple prefit both render as canvases; no raw-JSON input, no error', async () => {
    const errCount = await withConsoleErrorCount(page, async () => {
      await page.evaluate(() => { grok.shell.t.currentRowIdx = 2; grok.shell.t.currentRowIdx = 0; });
      await page.waitForTimeout(800);
    });
    const kinds = await page.evaluate((sel) => {
      const card = document.querySelector(sel);
      const s = card ? card.querySelector('[column="smiles"]') : null;
      const p = card ? card.querySelector('[column="multiple prefit"]') : null;
      return {smiles: s ? s.tagName : null, prefit: p ? p.tagName : null};
    }, CURRENT);
    expect(kinds.smiles).toBe('CANVAS');
    expect(kinds.prefit).toBe('CANVAS');
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  await softStep('Step 4c — with three rows selected, every drawn card holds both a molecule and a curve canvas', async () => {
    await page.evaluate(() => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;

      df.currentRowIdx = 0;
      df.selection.setAll(false);
      [2, 5, 8].forEach((r) => df.selection.set(r, true));
    });

    await expect.poll(async () => page.evaluate((sel) => {
      const df = grok.shell.t;
      let selCount = 0;
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selCount++;
      const expected = (df.currentRowIdx >= 0 ? 1 : 0) + selCount;
      const cards = Array.from(document.querySelectorAll(sel));
      let both = 0; let smilesOnly = 0; let prefitOnly = 0;
      for (const c of cards) {
        const s = c.querySelector('[column="smiles"]');
        const p = c.querySelector('[column="multiple prefit"]');
        const sc = !!s && s.tagName === 'CANVAS';
        const pc = !!p && p.tagName === 'CANVAS';
        if (sc && pc) both++;
        else if (sc) smilesOnly++;
        else if (pc) prefitOnly++;
      }

      return JSON.stringify({
        bothMatchesExpected: both === expected, expectedPositive: expected > 0, smilesOnly, prefitOnly,
      });
    }, ORDINARY), {timeout: 20_000})
      .toBe(JSON.stringify({
        bothMatchesExpected: true, expectedPositive: true, smilesOnly: 0, prefitOnly: 0,
      }));
  });

  const patterns = ['c1ccncc1', 'c1ccccc1', 'C(=O)O', 'C(=O)N', 'Cl', 'F', '[nH]', 'c1ccc2ccccc2c1'];
  let idCol = '';
  let pattern = '';
  let survivorId = '';

  let filterRepaintErrors = 0;
  const filterErrHandler = (msg: {type(): string}) => { if (msg.type() === 'error') filterRepaintErrors++; };

  await softStep('Scenario 5 setup — open spgi-100, add Forms, select four rows split by a substructure query', async () => {
    await v.openTable(page, {path: spgiPath, semTypeTimeoutMs: 3000, settleMs: 3000});
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});

    const defaults = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return {sel: vw.props.showSelectedRows, cur: vw.props.showCurrentRow, mo: vw.props.showMouseOverRow};
    });
    expect(defaults).toEqual({sel: true, cur: true, mo: true});

    const setup = await page.evaluate(async ({pats}) => {
      const df = grok.shell.t;
      const names = df.columns.names();
      let id: string | null = null;
      for (const n of names) {
        const col = df.col(n);
        if (col.semType !== 'Molecule' && col.type === 'string') { id = n; break; }
      }
      if (!id)
        for (const n of names) {
          const col = df.col(n);
          if (col.semType !== 'Molecule' && (col.type === 'int' || col.type === 'bigint')) { id = n; break; }
        }
      if (!id) return {error: 'no non-molecule identifier column on spgi-100'};

      let chosen: string | null = null; let matched: number[] = []; let unmatched: number[] = [];
      for (const p of pats) {
        const bs = await grok.chem.searchSubstructure(df.col('Structure'), p);
        const tc = bs.trueCount;
        if (tc > 0 && tc < df.rowCount) {
          chosen = p; matched = []; unmatched = [];
          for (let i = 0; i < df.rowCount; i++) (bs.get(i) ? matched : unmatched).push(i);
          break;
        }
      }
      if (!chosen || matched.length < 1 || unmatched.length < 1)
        return {error: 'no substructure pattern produced a proper subset on spgi-100'};

      const selMatched = matched.slice(0, 2);
      const selUnmatched = unmatched.slice(0, 2);
      const selected = [...selMatched, ...selUnmatched];
      df.selection.setAll(false);
      for (const r of selected) df.selection.set(r, true);
      df.mouseOverRowIdx = -1;
      let current = 0;
      while (selected.includes(current)) current++;
      df.currentRowIdx = current;

      const survivorRow = selMatched[0];
      return {
        id, pattern: chosen,
        survivorId: df.col(id).getString(survivorRow),
        selectedCount: selected.length,
      };
    }, {pats: patterns});

    if ((setup as any).error) throw new Error(`[Scenario 5 setup] ${(setup as any).error}`);
    idCol = (setup as any).id;
    pattern = (setup as any).pattern;
    survivorId = (setup as any).survivorId;

    await expect.poll(async () => page.evaluate(({sel, idc}) => {
      const cards = Array.from(document.querySelectorAll(sel)).slice(2);
      return cards.map((c) => {
        const e = c.querySelector(`[column="${idc}"]`) as HTMLInputElement | null;
        return e ? e.value : null;
      }).filter((x) => x !== null).length;
    }, {sel: ORDINARY, idc: idCol}), {timeout: 20_000}).toBe(4);

    await expect.poll(async () => page.evaluate(({sel, idc, idv}) => {
      for (const c of Array.from(document.querySelectorAll(sel))) {
        const e = c.querySelector(`[column="${idc}"]`) as HTMLInputElement | null;
        if (e && e.value === idv) {
          const cv = c.querySelector('[column="Structure"]') as HTMLCanvasElement | null;
          if (cv && cv.tagName === 'CANVAS') { cv.setAttribute('data-repaint-probe', 'forms5b'); return true; }
        }
      }
      return false;
    }, {sel: ORDINARY, idc: idCol, idv: survivorId}), {timeout: 20_000}).toBe(true);

    filterRepaintErrors = 0;
    page.on('console', filterErrHandler);
  });

  await softStep('Step 5a — after the substructure filter, the selected-row cards equal selection ∩ filter', async () => {
    await page.evaluate(async ({pat}) => {
      const df = grok.shell.t;
      const bs = await grok.chem.searchSubstructure(df.col('Structure'), pat);
      df.filter.init((i: number) => bs.get(i));
    }, {pat: pattern});

    await expect.poll(async () => page.evaluate(({sel, idc}) => {
      const df = grok.shell.t;
      const inter: string[] = []; let selCount = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.selection.get(i)) selCount++;
        if (df.selection.get(i) && df.filter.get(i)) inter.push(df.col(idc).getString(i));
      }
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => (c.querySelector(`[column="${idc}"]`) as HTMLInputElement | null)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({excluded: inter.length < selCount, match: JSON.stringify(inter) === JSON.stringify(tail)});
    }, {sel: ORDINARY, idc: idCol}), {timeout: 25_000})
      .toBe(JSON.stringify({excluded: true, match: true}));
  });

  await softStep('Step 5b — the surviving card\'s Structure canvas is rebuilt and drawn (not empty) after the filter', async () => {

    page.off('console', filterErrHandler);
    expect(filterRepaintErrors).toBe(0);

    await expect.poll(async () => page.evaluate(({sel, idc, idv}) => {
      const markerGone = !document.querySelector('[data-repaint-probe="forms5b"]');
      let present = false; let freshCanvas = false; let drawn = 0;
      for (const c of Array.from(document.querySelectorAll(sel))) {
        const e = c.querySelector(`[column="${idc}"]`) as HTMLInputElement | null;
        if (e && e.value === idv) {
          present = true;
          const cv = c.querySelector('[column="Structure"]') as HTMLCanvasElement | null;
          if (cv && cv.tagName === 'CANVAS' && !cv.hasAttribute('data-repaint-probe')) {
            freshCanvas = true;
            try {
              const d = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
              for (let i = 0; i < d.length; i += 4) {
                const a = d[i + 3]; const r = d[i]; const g = d[i + 1]; const b = d[i + 2];
                if (a > 10 && !(r >= 245 && g >= 245 && b >= 245)) drawn++;
              }
            } catch (_) { drawn = -1; }
          }
        }
      }
      return JSON.stringify({markerGone, present, freshCanvas, hasContent: drawn > 20});
    }, {sel: ORDINARY, idc: idCol, idv: survivorId}), {timeout: 20_000})
      .toBe(JSON.stringify({markerGone: true, present: true, freshCanvas: true, hasContent: true}));
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
