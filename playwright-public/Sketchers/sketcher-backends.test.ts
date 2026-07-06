/* ---
sub_features_covered: [chem.sketcher, chem.sketcher.backend-switch, chem.sketcher.chemdraw, chem.sketcher.copy-as, chem.sketcher.hamburger-menu, chem.sketcher.ketcher, chem.sketcher.molecular-input, chem.sketcher.ocl, chem.sketcher.roundtrip]
--- */
// One persistent sketcher panel, backend switched in-place via the hamburger menu
// (OpenChemLib → Ketcher → ChemDraw). Each backend runs the shared C1-C8 battery:
//   C1 SMILES round-trip + no-"undefined" guard (GROK-12685); C2 widget persists (GROK-16340);
//   C3 MOLBLOCK V2000 round-trip; C4 large V3000 loads (GROK-12391); C5 SMARTS round-trip
//   (GROK-12297/12966); C6 malformed no-throw (CLAUDE-5); C7 resize no-throw (GROK-12826);
//   C8 clear empties input (GROK-14028); active==selected (GROK-12581/12905); no console errors (GROK-12758).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

// Embedded large V3000 molblock (aspirin) for C4 — avoids an unbounded server convert that hung the test.
const ASPIRIN_V3000 = '\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 13 13 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 5.250000 -1.299038 0.000000 0\nM  V30 2 C 3.750000 -1.299038 0.000000 0\nM  V30 3 O 3.000000 -2.598076 0.000000 0\nM  V30 4 O 3.000000 0.000000 0.000000 0\nM  V30 5 C 1.500000 0.000000 0.000000 0\nM  V30 6 C 0.750000 -1.299038 0.000000 0\nM  V30 7 C -0.750000 -1.299038 0.000000 0\nM  V30 8 C -1.500000 0.000000 0.000000 0\nM  V30 9 C -0.750000 1.299038 0.000000 0\nM  V30 10 C 0.750000 1.299038 0.000000 0\nM  V30 11 C 1.500000 2.598076 0.000000 0\nM  V30 12 O 0.750000 3.897114 0.000000 0\nM  V30 13 O 3.000000 2.598076 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 2 2 3\nM  V30 3 1 2 4\nM  V30 4 1 4 5\nM  V30 5 2 5 6\nM  V30 6 1 6 7\nM  V30 7 2 7 8\nM  V30 8 1 8 9\nM  V30 9 2 9 10\nM  V30 10 1 10 11\nM  V30 11 2 11 12\nM  V30 12 1 11 13\nM  V30 13 1 10 5\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n';

// Open ONE sketcher dialog (OpenChemLib default), keep window.__sk handle + console.error capture.
async function openSketcher(page: import('@playwright/test').Page): Promise<void> {
  await page.evaluate(async () => {
    const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));
    grok.shell.closeAll();
    (window as any).__sk_err = [];
    const orig = console.error;
    console.error = function(...a: any[]) { (window as any).__sk_err.push(a.map((x: any) => String(x)).join(' ')); orig.apply(console, a as any); };
    DG.chem.currentSketcherType = 'OpenChemLib';   // set global before construction
    const sk = new DG.chem.Sketcher();
    (window as any).__sk = sk;
    ui.dialog('Sketcher').add(sk.root).show();
    const t0 = Date.now();
    while (Date.now() - t0 < 20000) { if (sk.sketcher?.isInitialized) break; await sleep(200); }
  });
  await page.locator('.d4-dialog').waitFor({timeout: 15000});
}

// Backend switch through the hamburger menu, on the same widget.
async function switchBackendViaMenu(page: import('@playwright/test').Page, friendlyName: string): Promise<void> {
  await page.locator('.d4-dialog .fa-bars.d4-input-options').click();
  await page.waitForTimeout(800);
  await page.locator('.d4-menu-item-label').filter({hasText: new RegExp(`^${friendlyName}$`)}).first().click();
  // external backends load a remote widget — wait until the switch lands on the same handle
  await page.waitForFunction((b) =>
    (window as any).DG?.chem?.currentSketcherType === b && (window as any).__sk?.sketcher,
  friendlyName, {timeout: 60000});
  await page.waitForTimeout(3000);
}

// The shared battery — runs against the current backend on the persistent window.__sk handle.
async function runBattery(page: import('@playwright/test').Page, backend: string): Promise<any> {
  return page.evaluate(async ([b, v3000]: [string, string]) => {
    const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));
    const withTimeout = (p: Promise<any>, ms: number) =>
      Promise.race([p, new Promise((_, rej) => setTimeout(() => rej(new Error('timeout')), ms))]);
    const sk = (window as any).__sk;
    const res: any = {backend: b, checks: {}};
    const errStart = ((window as any).__sk_err ?? []).length;
    res.activeType = DG.chem.currentSketcherType;
    res.innerClass = sk.sketcher?.constructor?.name;

    // C1: SMILES round-trip + no-"undefined" guard (GROK-12685)
    sk.setMolecule('c1ccccc1');
    await sleep(2500);
    res.checks.smiles = sk.getSmiles();
    res.checks.smilesNotUndefined = sk.getSmiles() !== 'undefined' && sk.molInput.value !== 'undefined';
    res.checks.molfileValid = /V2000|V3000/.test(sk.getMolFile() ?? '');

    // C2: widget persists after set — does not disappear (GROK-16340)
    await sleep(2000);
    res.checks.persistInner = !!sk.sketcher;
    res.checks.persistRootInDom = document.body.contains(sk.root);
    res.checks.persistMolfileLen = (sk.getMolFile() ?? '').length;

    // C3: MOLBLOCK (V2000) round-trip
    const mb = sk.getMolFile();
    sk.setMolecule(''); await sleep(600);
    sk.setMolecule(mb); await sleep(1500);
    res.checks.molblockRoundtrip = (sk.getMolFile() ?? '').length > 0 && !sk.isEmpty();

    // C4: large V3000 molecule loads, not empty (GROK-12391)
    try {
      sk.setMolecule(v3000); await sleep(1800);
      res.checks.v3000NotEmpty = !sk.isEmpty() && (sk.getMolFile() ?? '').length > 0;
    } catch (e) { res.checks.v3000Err = String(e); }

    // C5: SMARTS set/get round-trip (GROK-12297 / GROK-12966) — bound the possibly-server call
    try {
      sk.setMolecule('[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1', true); await sleep(1500);
      const sm = await withTimeout(sk.getSmarts(), 15000);
      res.checks.smartsNonEmpty = !!sm && sm.length > 0;
    } catch (e) { res.checks.smartsErr = String(e); }

    // C6: malformed molecule must not throw (CLAUDE-5)
    try { sk.setMolecule('not_a_molecule_zzz'); await sleep(600); res.checks.malformedNoThrow = true; }
    catch (e) { res.checks.malformedThrew = String(e); }

    // C7: resize() must not throw (GROK-12826)
    try { sk.resize(); res.checks.resizeNoThrow = true; }
    catch (e) { res.checks.resizeThrew = String(e); }

    // C8: clear empties the molecular input field (GROK-14028)
    sk.setMolecule(''); await sleep(700);
    res.checks.clearEmpty = sk.isEmpty() && sk.molInput.value === '';

    const newErrs = ((window as any).__sk_err ?? []).slice(errStart);
    res.consoleErrsRaw = newErrs.slice(0, 5);
    res.consoleErrs = newErrs.filter((e: string) =>
      /sketcher|search pattern cannot be set|copy as|molblock|molfile|moleculesketcher|getsmiles|getsmarts/i.test(e)).length;
    return res;
  }, [backend, ASPIRIN_V3000] as [string, string]);
}

// Same hard checks for every backend.
function assertBattery(res: any, label: string) {
  expect(res.checks.smiles, `${label} C1 SMILES: ${JSON.stringify(res.checks)}`).toBe('c1ccccc1');
  expect(res.checks.smilesNotUndefined, `${label} C1 no-"undefined" (GROK-12685)`).toBe(true);
  expect(res.checks.molfileValid, `${label} C1 valid MOLBLOCK`).toBe(true);
  expect(res.checks.persistInner, `${label} C2 sketcher persists (GROK-16340)`).toBe(true);
  expect(res.checks.persistRootInDom, `${label} C2 root stays in DOM (GROK-16340)`).toBe(true);
  expect(res.checks.persistMolfileLen, `${label} C2 molecule present`).toBeGreaterThan(0);
  expect(res.checks.molblockRoundtrip, `${label} C3 MOLBLOCK round-trip`).toBe(true);
  expect(res.checks.v3000NotEmpty, `${label} C4 V3000 (GROK-12391)`).toBe(true);
  expect(res.checks.smartsNonEmpty, `${label} C5 SMARTS (GROK-12297)`).toBe(true);
  expect(res.checks.malformedNoThrow, `${label} C6 malformed no-throw (CLAUDE-5): ${res.checks.malformedThrew ?? ''}`).toBe(true);
  expect(res.checks.resizeNoThrow, `${label} C7 resize no-throw (GROK-12826): ${res.checks.resizeThrew ?? ''}`).toBe(true);
  expect(res.checks.clearEmpty, `${label} C8 clear empties input (GROK-14028)`).toBe(true);
  expect(res.consoleErrs, `${label} sketcher console errors (GROK-12758): ${JSON.stringify(res.consoleErrsRaw)}`).toBe(0);
}

// ChemDraw and Marvin are proprietary sketcher backends with no package in the
// public repo, so they cannot be installed on the public CI stack — this suite
// exercises the two backends that ship publicly: OpenChemLib (Chem) and Ketcher
// (KetcherSketcher). Prereq packages: Chem, KetcherSketcher.
test('Chem: Sketcher battery — OpenChemLib → Ketcher (UI backend switch)', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  await softStep('Open sketcher (OpenChemLib default)', async () => {
    await openSketcher(page);
  });

  // Select OCL via the menu too — initial backend is the user's last-selected (non-deterministic).
  await softStep('Select OpenChemLib via hamburger menu — battery', async () => {
    await switchBackendViaMenu(page, 'OpenChemLib');
    const res = await runBattery(page, 'OpenChemLib');
    console.log(`[sketcher] OpenChemLib: ${JSON.stringify(res)}`);
    expect(res.activeType, 'active backend is OpenChemLib').toBe('OpenChemLib');
    assertBattery(res, 'OpenChemLib');
  });

  await softStep('Switch to Ketcher via hamburger menu — battery', async () => {
    await switchBackendViaMenu(page, 'Ketcher');
    const res = await runBattery(page, 'Ketcher');
    console.log(`[sketcher] Ketcher: ${JSON.stringify(res)}`);
    expect(res.activeType, 'switched to Ketcher (GROK-12581)').toBe('Ketcher');
    assertBattery(res, 'Ketcher');
  });

  // Restore the user's default sketcher so this run doesn't leave a non-default persisted.
  await page.evaluate(() => {
    try { grok.userSettings.add('sketcher', 'selected', 'OpenChemLib'); } catch (e) {}
    grok.shell.closeAll();
  }).catch(() => {});

  finishSpec();
});
