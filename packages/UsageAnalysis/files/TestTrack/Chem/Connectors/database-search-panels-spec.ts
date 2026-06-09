import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Chem | Context Panel — External Database Search Panels', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const g = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
    (window as any).__df = df;
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 4000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
  await page.locator('[name="viewer-Grid"] canvas').first().waitFor({timeout: 30_000});

  await page.evaluate(() => {
    const w = window as any; const DG = w.DG; const g = w.grok;
    w.__selectMol = async (smiles: string) => {
      g.shell.o = DG.SemanticValue.fromValueType(smiles, 'Molecule');
      await new Promise((r) => setTimeout(r, 4000));
    };
    w.__expand = async (root: Element, title: string, waitMs = 2000) => {
      const pane = root.querySelector(`.d4-accordion-pane[d4-title="${title}"]`);
      if (!pane) return null;
      const header = pane.querySelector('.d4-accordion-pane-header') as HTMLElement;
      if (!header.classList.contains('expanded')) header.click();
      await new Promise((r) => setTimeout(r, waitMs));
      return pane;
    };
    w.__dbContent = () => document
      .querySelector('.grok-prop-panel .d4-accordion-pane[d4-title="Databases"] .d4-accordion-pane-content');
    w.__dbInit = async () => {
      for (let i = 0; i < 10; i++) {
        try { await g.functions.call('DrugBank:initDrugBank'); return; }
        catch (e) { await new Promise((r) => setTimeout(r, 2000)); }
      }
    };
  });

  await softStep('Scenario 1 — Databases group appears on molecule cell', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any; const df = w.__df;
      await w.__selectMol(df.col('canonical_smiles').get(0));
      const groups = Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane'))
        .map((p) => p.getAttribute('d4-title'));
      const top = ['Actions', 'Chemistry', 'Databases', 'Biology', 'Structure'].filter((t) => groups.includes(t));
      const dbPane = await w.__expand(document.querySelector('.grok-prop-panel'), 'Databases', 3000);
      const providersAll = Array.from(dbPane.querySelectorAll('.d4-accordion-pane'))
        .map((p: any) => p.getAttribute('d4-title'));
      return {top, hasFour: ['ChEMBL', 'Chemspace', 'PubChem', 'DrugBank'].every((p) => providersAll.includes(p))};
    });
    expect(res.top).toEqual(['Actions', 'Chemistry', 'Databases', 'Biology', 'Structure']);
    expect(res.hasFour).toBe(true);
  });

  await softStep('Scenario 2 — ChEMBL API search panels', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const chembl = await w.__expand(w.__dbContent(), 'ChEMBL', 2500);
      const subTitles = Array.from(chembl.querySelectorAll('.d4-accordion-pane'))
        .map((c: any) => c.getAttribute('d4-title')).filter(Boolean);
      const sim = await w.__expand(chembl, 'Similarity Search API', 9000);
      const simC = sim.querySelector('.d4-accordion-pane-content');
      const simCards = simC.querySelectorAll('canvas').length;
      const simHasScore = /Score:\s*1\.00/.test(simC.innerText);
      const openIcon = !!simC.querySelector('i.fa-arrow-square-down');
      const ss = await w.__expand(chembl, 'Substructure Search API', 10000);
      const ssC = ss.querySelector('.d4-accordion-pane-content');
      const ssOk = ssC.querySelectorAll('canvas').length > 0 || /No matches/.test(ssC.innerText) || ssC.innerText.trim() === '';
      const viewsBefore = Array.from(w.grok.shell.views).length;
      simC.querySelector('i.fa-arrow-square-down').click();
      await new Promise((r) => setTimeout(r, 3000));
      const viewNames = Array.from(w.grok.shell.views).map((v: any) => v.name);
      return {subTitles, simCards, simHasScore, openIcon, ssOk, viewsBefore, viewsAfter: viewNames.length,
        opened: viewNames.some((n: string) => /ChEMBL/i.test(n))};
    });
    expect(res.subTitles.some((t: string) => t.includes('Similarity Search API'))).toBe(true);
    expect(res.subTitles.some((t: string) => t.includes('Substructure Search API'))).toBe(true);
    expect(res.simCards).toBeGreaterThan(0);
    expect(res.simHasScore).toBe(true);
    expect(res.openIcon).toBe(true);
    expect(res.ssOk).toBe(true);
    expect(res.viewsAfter).toBeGreaterThan(res.viewsBefore);
    expect(res.opened).toBe(true);
  });

  await softStep('Scenario 3 — Chemspace form + Similar/Substructure', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const views = Array.from(w.grok.shell.views);
      w.grok.shell.v = views.find((v: any) => v.name === 'Table');
      await new Promise((r) => setTimeout(r, 1000));
      await w.__selectMol(w.__df.col('canonical_smiles').get(0));
      await w.__expand(document.querySelector('.grok-prop-panel'), 'Databases', 1500);
      const chembl = w.__dbContent().querySelector('.d4-accordion-pane[d4-title="ChEMBL"] .d4-accordion-pane-header');
      if (chembl && chembl.classList.contains('expanded')) chembl.click();
      const cs = await w.__expand(w.__dbContent(), 'Chemspace', 10000);
      const csC = cs.querySelector('.d4-accordion-pane-content');
      const inputs = Array.from(csC.querySelectorAll('.ui-input-root'))
        .map((r: any) => r.querySelector('.ui-input-label,label')?.textContent).filter(Boolean);
      const nested = Array.from(csC.querySelectorAll('.d4-accordion-pane')).map((p: any) => p.getAttribute('d4-title'));
      const authError = /unauthor|forbidden|invalid token|api key/i.test(csC.innerText);
      const similar = await w.__expand(cs, 'Similar', 2000);
      const simC = similar.querySelector('.d4-accordion-pane-content');
      const resolved = () => simC.querySelectorAll('canvas').length > 0 || /No matches|Score/.test((simC as HTMLElement).innerText);
      let simResolved = false;
      for (let i = 0; i < 90; i++) {
        if (resolved()) { simResolved = true; break; }
        await new Promise((r) => setTimeout(r, 500));
      }
      return {inputs, nested, authError, simResolved};
    });
    expect(res.inputs).toContain('Ship to country');
    expect(res.inputs).toContain('Category');
    expect(res.nested).toContain('Similar');
    expect(res.nested).toContain('Substructure');
    expect(res.authError).toBe(false);
    expect(res.simResolved).toBe(true);
  });

  await softStep('Scenario 4 — PubChem search panels resolve without crash', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const cs = w.__dbContent().querySelector('.d4-accordion-pane[d4-title="Chemspace"] .d4-accordion-pane-header');
      if (cs && cs.classList.contains('expanded')) cs.click();
      const pc = await w.__expand(w.__dbContent(), 'PubChem', 2500);
      const pcC = pc.querySelector('.d4-accordion-pane-content');
      const nested = Array.from(pcC.querySelectorAll('.d4-accordion-pane'))
        .map((c: any) => c.getAttribute('d4-title')).filter(Boolean);
      for (const p of Array.from(pcC.querySelectorAll('.d4-accordion-pane'))) {
        const h = (p as Element).querySelector('.d4-accordion-pane-header') as HTMLElement;
        if (!h.classList.contains('expanded')) h.click();
      }
      await new Promise((r) => setTimeout(r, 30000));
      const report = Array.from(pcC.querySelectorAll('.d4-accordion-pane')).map((p: any) => {
        const c = p.querySelector('.d4-accordion-pane-content');
        const resolved = c.querySelectorAll('canvas').length > 0 || /No matches|Score|Compound/i.test(c.innerText);
        return {title: p.getAttribute('d4-title'), resolved};
      });
      return {nested, report};
    });
    expect(res.nested.some((t: string) => t.includes('Substructure Search'))).toBe(true);
    expect(res.nested.some((t: string) => t.includes('Similarity Search'))).toBe(true);
    expect(res.nested.some((t: string) => t.includes('Identity Search'))).toBe(true);
    for (const p of res.report) expect(p.resolved).toBe(true);
  });

  await softStep('Scenario 5 — DrugBank substructure + similarity on a known drug', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      await w.__dbInit();
      const pc = w.__dbContent().querySelector('.d4-accordion-pane[d4-title="PubChem"] .d4-accordion-pane-header');
      if (pc && pc.classList.contains('expanded')) pc.click();
      await w.__selectMol('CC(=O)Oc1ccccc1C(=O)O');
      await w.__expand(document.querySelector('.grok-prop-panel'), 'Databases', 1500);
      const drb = await w.__expand(w.__dbContent(), 'DrugBank', 2500);
      const subTitles = Array.from(drb.querySelectorAll('.d4-accordion-pane'))
        .map((c: any) => c.getAttribute('d4-title')).filter(Boolean);
      const sim = await w.__expand(drb, 'Similarity Search', 2000);
      const ss = await w.__expand(drb, 'Substructure Search', 2000);
      const simC = sim.querySelector('.d4-accordion-pane-content');
      const ssC = ss.querySelector('.d4-accordion-pane-content');
      const done = (c: Element) => c.querySelectorAll('canvas').length > 0 || /No matches/.test((c as HTMLElement).innerText);
      for (let i = 0; i < 80; i++) {
        if (done(simC) && done(ssC)) break;
        await new Promise((r) => setTimeout(r, 500));
      }
      return {subTitles,
        simCards: simC.querySelectorAll('canvas').length, simScore: /Score:\s*1\.00/.test(simC.innerText),
        ssCards: ssC.querySelectorAll('canvas').length};
    });
    expect(res.subTitles.some((t: string) => t.includes('Substructure Search'))).toBe(true);
    expect(res.subTitles.some((t: string) => t.includes('Similarity Search'))).toBe(true);
    expect(res.simCards).toBeGreaterThan(0);
    expect(res.simScore).toBe(true);
    expect(res.ssCards).toBeGreaterThan(0);
  });

  await softStep('Scenario 6 — empty molecule handled gracefully', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      await w.__selectMol('');
      await w.__expand(document.querySelector('.grok-prop-panel'), 'Databases', 1500);
      const drb = await w.__expand(w.__dbContent(), 'DrugBank', 2500);
      const txt = drb ? drb.querySelector('.d4-accordion-pane-content').innerText : '';
      return {emptyHandled: /SMILES is empty/.test(txt)};
    });
    expect(res.emptyHandled).toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
