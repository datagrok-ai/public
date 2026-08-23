/* ---
realizes: [chem.calculate.mpo-score]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

declare const grok: any;
declare const DG: any;

// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser reference):
//   .chem-mpo-action-button (text "Create profile") — MPO Profiles grok-view (Chem:mpoProfilesApp),
//     first of two action buttons (Create profile / Upload); observed live 2026-08-19 via chrome-devtools MCP. Not in chem.md.
//   .chem-mpo-profiles-table a.ui-link — profile-name links in the MPO Profiles list table; observed 2026-08-19. Not in chem.md.
//   .chem-mpo-actions-button (text "⋮", [name="button-⋮"]) — per-row actions button in the profiles table,
//     opens a .d4-menu-popup with Edit/Clone/Download/Delete .d4-menu-item-label items; observed 2026-08-19. Not in chem.md.
//   .chem-profile-header (contentEditable h1) — profile-name field in the create/edit view; observed 2026-08-19. Not in chem.md.
//   .chem-profile-editor-container .statistics-mpo-property-cell input.ui-input-editor — property-name input in the
//     desirability editor; add-property affordance is a.ui-link "+ Add Property" inside .statistics-mpo-empty-state;
//     observed 2026-08-19. Not in chem.md.
//   .d4-dialog "Delete profile" + [name="button-OK"] — delete-confirmation dialog; observed 2026-08-19. Not in chem.md.
//   .d4-ribbon-name .ui-breadcrumbs-text-element — breadcrumb segments in the profile editor's ribbon
//     (Home / MPO Profiles / <name>). Source-verified, not MCP-observed: ui.breadcrumbs builds each segment as
//     ui.link(..., 'ui-breadcrumbs-text-element') (js-api/src/widgets/specialized.ts:44) and setupMpoBreadcrumbs
//     injects the root into .d4-ribbon-name (Chem/src/mpo/utils.ts:323); clicking the "MPO Profiles" segment sets
//     grok.shell.v to the already-open list view (utils.ts:310). Not in chem.md.

const MPO_DIR = 'System:AppData/Chem/mpo';
const SYNC_PROFILE = 'SyncTest-Profile';
const SCORE_PROFILE = 'ScoreTest-Profile';
const SCORE_COLUMN = 'HeavyAtomCount';

async function openMpoProfilesApp(page: Page): Promise<void> {
  await page.evaluate(async () => {
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    const v = await grok.functions.call('Chem:mpoProfilesApp', {});
    grok.shell.addView(v);
  });
  await page.locator('.chem-mpo-action-button').first().waitFor({timeout: 30_000, state: 'attached'});
  await page.locator('.chem-mpo-profiles-table').waitFor({timeout: 30_000, state: 'attached'});
}

// A profile's file name is DERIVED from its name, not equal to it: generateMpoFileName
// (libraries/statistics/src/mpo/utils.ts:21) slugifies whitespace and slashes and appends
// -2, -3 ... on collision. Asserting on the file name therefore tests the slug rule, and
// misses the case that matters — a right-named file holding the wrong profile name. Every
// identity check below goes through the JSON's internal `name`.
// `properties` carries the profile's rule keys because the name alone cannot tell a real save
// from an inert one: the property editor opens pre-filled with the registered default
// `NewProperty1` (chem.md:1460), so a profile whose name committed but whose property field
// never did persists as valid JSON under the right name with the wrong single rule.
interface DiskProfile {
  fileName: string;
  name: string | null;
  properties: string[] | null;
  error: string | null;
}

// files.list cannot answer this question: it is served through wrapCached
// (datlas/lib/src/routers/connectors.dart:260-305) and stamped
// `Cache-Control: max-age=<up to a year>, immutable` (:290-292), so the browser re-fetches only
// when the `flag` query parameter changes — re-reading the same URL returns the same listing
// however long you wait. readFilesAsString goes to the /folder/ route (connectors.dart:32 ->
// getConnectionFolder :745), which reads storage directly and answers with Content-Type only,
// so neither layer can cache it. It also hands back the file contents, which is what identity
// needs here (see the note on DiskProfile above).
//
// The JSON is parsed here rather than by readFilesAsJson because that method swallows a parse
// failure and omits the file (js-api/src/dapi.ts:2150-2152) — a corrupt profile would read back
// as one that is not there, and the absence would look clean. Parsing in the spec keeps it an
// explicit unreadable for assertAllProfilesReadable to trip on.
//
// An unreachable or non-204 error response throws out of readFilesAsBlobs
// (js-api/src/dapi.ts:2116-2136) rather than degrading to an empty object, so [] here means an
// HTTP 204 — a directory the server really did report as empty — never a read that failed.
async function readDiskProfiles(page: Page): Promise<DiskProfile[]> {
  return await page.evaluate(async (dir) => {
    const byFile: {[fileName: string]: string} = await grok.dapi.files.readFilesAsString(dir, false, 'json');
    return Object.entries(byFile).map(([fileName, text]: [string, string]) => {
      try {
        const parsed = JSON.parse(text);
        const props = parsed?.properties;
        return {
          fileName,
          name: typeof parsed?.name === 'string' ? parsed.name : null,
          properties: props !== null && typeof props === 'object' ? Object.keys(props) : null,
          error: null,
        };
      } catch (e) {
        return {fileName, name: null, properties: null, error: String(e)};
      }
    });
  }, MPO_DIR);
}

// A profile that could not be read is neither present nor absent, and an absence
// assertion would pass on it silently. Both callers of a negative check run this first.
function assertAllProfilesReadable(profiles: DiskProfile[]): void {
  const broken = profiles.filter((p) => p.error !== null || p.name === null);
  expect(broken.length,
    `unreadable MPO profile files — presence/absence cannot be decided: ${JSON.stringify(broken)}`).toBe(0);
}

// The candidate read goes through readDiskProfiles, so it is the uncached /folder/ route and
// cannot report a deleted profile as still present. files.exists stays the authority on each
// candidate anyway — it is a HEAD with no `flag` (files_client.dart:56-59,
// connectors.dart:499-527), so the two independent uncached reads must agree before this
// reports a file alive. Costs one call per candidate, normally none.
// This read is its own, not the callers' — it must carry the callers' guards too. An HTTP-204
// (empty) or corrupt-name read yields no candidates, files.exists is never called, and the
// caller's toEqual([]) reports a clean delete having observed nothing at all.
async function profileFilesStillOnDisk(page: Page, profileName: string): Promise<string[]> {
  const profiles = await readDiskProfiles(page);
  assertAllProfilesReadable(profiles);
  expect(profiles.length,
    `the profiles directory read back empty while checking for ${profileName} — an absence ` +
    'check against nothing proves nothing').toBeGreaterThan(0);
  const candidates = profiles
    .filter((p) => p.name === profileName)
    .map((p) => p.fileName);
  return await page.evaluate(async ({dir, files}) => {
    const alive: string[] = [];
    for (const f of files)
      if (await grok.dapi.files.exists(`${dir}/${f}`)) alive.push(f);
    return alive;
  }, {dir: MPO_DIR, files: candidates});
}

async function deleteDiskProfilesNamed(page: Page, profileName: string): Promise<void> {
  const victims = (await readDiskProfiles(page))
    .filter((p) => p.name === profileName)
    .map((p) => p.fileName);
  await page.evaluate(async ({dir, files}) => {
    for (const f of files) {
      try { await grok.dapi.files.delete(`${dir}/${f}`); } catch (e) {}
    }
  }, {dir: MPO_DIR, files: victims});
}

async function deleteDiskProfile(page: Page, fileName: string): Promise<void> {
  await page.evaluate(async ({dir, fileName}) => {
    try { await grok.dapi.files.delete(`${dir}/${fileName}`); } catch (e) {}
  }, {dir: MPO_DIR, fileName});
}

// Expanding is a Setup concern and happens exactly once, here. It is separated from the
// read because assigning `expanded` can re-invoke Chem's tree-browser handler, whose
// refresh() drops every item and repopulates from MpoProfileManager.ensureLoaded() — the
// in-memory cache (Chem/src/package.ts:3170-3175). The UI delete empties that cache one
// line BEFORE it fires its event (mpo-profile-manager.ts:69-72), so a read that re-expands
// would rebuild the list from an already-emptied cache and satisfy the Step-10 guard even
// if the event never fired. Priming once and reading read-only removes that branch
// entirely: after this call no assertion path ever touches `expanded` again.
async function primeBrowseTree(page: Page): Promise<string> {
  return await page.evaluate(async () => {
    const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
    const mainTree = grok.shell.browsePanel?.mainTree;
    if (!mainTree) return '__NO_TREE__';
    const apps = (mainTree.children || []).find((c: any) => c.text === 'Apps');
    if (!apps) return '__NO_APPS__';
    apps.expanded = true; await wait(1500);
    const chem = (apps.children || []).find((c: any) => c.text === 'Chem');
    if (!chem) return '__NO_CHEM__';
    chem.expanded = true; await wait(1800);
    const mpo = (chem.children || []).find((c: any) => /MPO/i.test(c.text));
    if (!mpo) return '__NO_MPO__';
    mpo.expanded = true; await wait(2500);
    return 'OK';
  });
}

// Strictly read-only: it locates nodes and reads `children`, and assigns nothing. That is
// what makes the Step-10 guard unarguable — the tree can only have lost the deleted entry
// because the product's own event drove refresh(), since nothing here could have rebuilt
// the item list. Verified end-to-end on dev 2026-08-22: after a real UI create the walk
// went 15 -> 16 and contained the new name, and after a real UI delete it went back to 15
// and did not, with the platform's changed/deleted events observed firing. Do not
// "simplify" this into a fresh Browse-panel build, and do not re-add an `expanded`
// assignment: either one refreshes the surface under test and makes the guard a tautology.
async function readBrowseTreeMpoProfiles(page: Page): Promise<string[]> {
  return await page.evaluate(() => {
    const mainTree = grok.shell.browsePanel?.mainTree;
    if (!mainTree) return ['__NO_TREE__'];
    const apps = (mainTree.children || []).find((c: any) => c.text === 'Apps');
    if (!apps) return ['__NO_APPS__'];
    const chem = (apps.children || []).find((c: any) => c.text === 'Chem');
    if (!chem) return ['__NO_CHEM__'];
    const mpo = (chem.children || []).find((c: any) => /MPO/i.test(c.text));
    if (!mpo) return ['__NO_MPO__'];
    return (mpo.children || []).map((c: any) => c.text);
  });
}

// The sentinels above are failed reads, not tree states, and a one-element ['__NO_MPO__']
// satisfies both `not.toContain(profile)` and a length comparison against a baseline that
// was itself a sentinel — the GROK-19624 guard would report satisfied on a tree it never
// saw. Every read is decided here before it is used for presence, absence, or a count.
function assertTreeRead(tree: string[]): void {
  expect(tree.every((t) => !t.startsWith('__')),
    `Browse-tree read failed — presence, absence and count are all undecidable: ${JSON.stringify(tree)}`).toBe(true);
  // `every` is vacuously true on [], and an empty tree passes `not.toContain(profile)` and a
  // 0-vs-0 count while having observed nothing at all. The baseline is asserted non-empty in
  // Setup, so an empty read here is a broken read, not a legitimately empty stand.
  expect(tree.length,
    'Browse-tree read came back empty — an absence check against nothing proves nothing').toBeGreaterThan(0);
}

// The 6 s is a deadline, not a grace. MpoProfileManager fires the event only on a SUCCEEDED
// save/delete (mpo-profile-manager.ts:146 and :72; both catch paths return without firing),
// so resolving on the timer means the event-driven propagation that IS the sync contract did
// not happen — and a settle that swallows which branch won cannot tell the two apart. The
// listener is installed in its own awaited round-trip so it is provably armed before the
// gesture, and the winning branch is returned for the caller to assert on.
async function actAndAwaitProfileEvent(page: Page, eventName: string, act: () => Promise<void>): Promise<string> {
  await page.evaluate((name) => {
    (window as any).__mpoSettled = new Promise((res) => {
      const sub = grok.events.onCustomEvent(name).subscribe(() => { sub.unsubscribe(); res('event'); });
      setTimeout(() => res('timeout'), 6000);
    });
  }, eventName);
  await act();
  return await page.evaluate(() => (window as any).__mpoSettled as Promise<string>);
}

async function listProfileLinkTexts(page: Page): Promise<string[]> {
  return await page.locator('.chem-mpo-profiles-table a.ui-link').allInnerTexts();
}

test.use(specTestOptions);

test('Chem: MPO profile CRUD and Browse tree sync (GROK-19624)', async ({page}) => {
  test.setTimeout(360_000);

  await loginToDatagrok(page);

  await softStep('Setup: remove any leftover test profiles, open the MPO Profiles app, record the baseline Browse-tree count', async () => {
    await deleteDiskProfilesNamed(page, SYNC_PROFILE);
    await deleteDiskProfilesNamed(page, SCORE_PROFILE);
    // Scenario 3 writes its profile to this exact path, so that one file is purged by path too.
    await deleteDiskProfile(page, `${SCORE_PROFILE}.json`);
    await openMpoProfilesApp(page);

    // Both baselines must be non-empty before anything downstream compares against them.
    // A stand with no profiles collapses every later absence check into a statement about
    // nothing: `expect([]).not.toContain(x)` and `expect(0).toBe(0)` both pass having
    // observed no tree and no directory. Dev carries 15 stored profiles (chem.md), so an
    // empty baseline is itself the signal that the read, not the stand, is broken.
    const baselineDisk = await readDiskProfiles(page);
    assertAllProfilesReadable(baselineDisk);
    expect(baselineDisk.length,
      'no MPO profiles on disk at baseline — every later disk absence check would be vacuous').toBeGreaterThan(0);

    const primed = await primeBrowseTree(page);
    expect(primed, 'Browse tree could not be primed to Apps > Chem > MPO profiles').toBe('OK');
    const baselineTree = await readBrowseTreeMpoProfiles(page);
    assertTreeRead(baselineTree);
    expect(baselineTree).not.toContain(SYNC_PROFILE);
    await page.evaluate((n) => { (window as any).__baseTreeCount = n; }, baselineTree.length);
  });

  await softStep('Scenario 1 Step 2-3: click "Create profile", add a property, name it, and Save', async () => {
    await page.locator('.chem-mpo-action-button', {hasText: 'Create profile'}).first().click();
    await page.locator('.chem-profile-header').waitFor({timeout: 30_000, state: 'attached'});

    await page.locator('.chem-profile-editor-container a.ui-link', {hasText: 'Add Property'}).click();
    const propInput = page.locator('.chem-profile-editor-container .statistics-mpo-property-cell input.ui-input-editor').first();
    await propInput.waitFor({timeout: 10_000, state: 'attached'});
    // Same actuation class as the header below: this editor is a Dart-backed ui.input.string
    // (mpo-profile-editor.ts:304 -> InputBase.forInputType, js-api/ui.ts:1100/:1051), and a
    // synthesized value write does not reach the model (chem.md § Actuating Chem inputs).
    // Click, select the pre-filled default, then type with real keystrokes.
    await propInput.click();
    await propInput.press('Control+a');
    await page.keyboard.type(SCORE_COLUMN);
    // chem.md:1344 — assert the box holds the intended name BEFORE the commit; a dropped
    // keystroke otherwise commits the registered default NewProperty1 (chem.md:1460), and
    // every downstream assertion in this scenario still passes on that profile.
    await expect(propInput).toHaveValue(SCORE_COLUMN, {timeout: 10_000});
    await propInput.press('Enter');

    const header = page.locator('.chem-profile-header');
    await header.click();
    // Control+A does not select the contenteditable's contents here, so trusted typing
    // appends onto "Untitled Profile" (Gate B failure: saved as "Untitled ProfileSyncTest-Profile").
    // Establish a full-content selection via the Selection API (what a triple-click does), then
    // trusted type replaces it — verified live: header commits exactly "SyncTest-Profile".
    await header.evaluate((h: HTMLElement) => {
      h.focus();
      const range = document.createRange();
      range.selectNodeContents(h);
      const sel = window.getSelection()!;
      sel.removeAllRanges();
      sel.addRange(range);
    });
    await page.keyboard.type(SYNC_PROFILE);
    await expect(header).toHaveText(SYNC_PROFILE, {timeout: 10_000});
    // Re-read at the commit point: naming the header re-renders the editor rows, so the value
    // asserted above has to still be there when Save reads the model.
    await expect(propInput).toHaveValue(SCORE_COLUMN, {timeout: 10_000});
    const saveButton = page.locator('[name="button-Save"]').first();
    await expect(saveButton).not.toHaveClass(/d4-disabled/, {timeout: 15_000});
    const settledBy = await actAndAwaitProfileEvent(page, 'chem-mpo-profile-changed',
      () => saveButton.click());
    expect(settledBy,
      'Save must announce itself through chem-mpo-profile-changed within 6 s — that event is what the ' +
      'MPO Profiles list and the Browse tree both subscribe to, so no event means no cross-surface sync')
      .toBe('event');
  });

  await softStep('Scenario 1 Step 4: the new profile is persisted and appears in the MPO Profiles list', async () => {
    const onDisk = await readDiskProfiles(page);
    console.log(`Step 4: ${MPO_DIR} returned ${onDisk.length} file(s); profile names: ` +
      JSON.stringify(onDisk.map((p) => p.name)));
    assertAllProfilesReadable(onDisk);
    const saved = onDisk.filter((p) => p.name === SYNC_PROFILE);
    expect(saved.length,
      `exactly one persisted profile must carry the name typed into the header; on disk: ${JSON.stringify(onDisk)}`).toBe(1);
    console.log(`Step 4: persisted properties of ${SYNC_PROFILE}: ${JSON.stringify(saved[0].properties)}`);
    // The name proves the header committed; only the rule keys prove the property field did.
    // A profile saved with the untouched default persists as {"NewProperty1": {...}} and would
    // satisfy every other assertion in this scenario.
    expect(saved[0].properties,
      `persisted profile has no properties object: ${JSON.stringify(saved[0])}`).not.toBeNull();
    expect(saved[0].properties!,
      `the persisted profile must carry the property typed into the editor, not the registered ` +
      `default NewProperty1; properties: ${JSON.stringify(saved[0].properties)}`).toContain(SCORE_COLUMN);
    // Returning to the list must NOT rebuild the view: closeAll() + a fresh
    // Chem:mpoProfilesApp call re-reads from disk, which is a stronger reset than the page
    // reload this anchor excludes, and it would hide the very failure the step is for — a
    // list that never picks up the saved profile. The editor's own breadcrumb switches to
    // the already-open list view instead (setupMpoBreadcrumbs, Chem/src/mpo/utils.ts:302).
    await page.locator('.d4-ribbon-name .ui-breadcrumbs-text-element', {hasText: 'MPO Profiles'})
      .first().click();
    await page.locator('.chem-mpo-profiles-table').waitFor({timeout: 30_000, state: 'attached'});
    const links = await listProfileLinkTexts(page);
    expect(links,
      `the open MPO Profiles list must show the saved profile without being rebuilt; links: ${JSON.stringify(links)}`)
      .toContain(SYNC_PROFILE);
  });

  await softStep('Scenario 1 Step 5: the Browse tree gains the new profile — count is baseline + 1, no manual refresh', async () => {
    const tree = await readBrowseTreeMpoProfiles(page);
    assertTreeRead(tree);
    const base = await page.evaluate(() => (window as any).__baseTreeCount as number);
    expect(tree).toContain(SYNC_PROFILE);
    expect(tree.length).toBe(base + 1);
  });

  await softStep('Scenario 2 Step 6-8: delete SyncTest-Profile via the row actions menu and confirm', async () => {
    const row = page.locator('.chem-mpo-profiles-table tr', {hasText: SYNC_PROFILE});
    await row.locator('.chem-mpo-actions-button').click();
    await page.locator('.d4-menu-popup, .d4-menu-item-label').first().waitFor({timeout: 10_000, state: 'attached'});
    await page.locator('.d4-menu-item-label', {hasText: /^Delete$/}).click();

    const deleteDialog = page.locator('.d4-dialog', {hasText: 'Delete profile'});
    await deleteDialog.waitFor({timeout: 10_000, state: 'attached'});
    const settledBy = await actAndAwaitProfileEvent(page, 'chem-mpo-profile-deleted',
      () => deleteDialog.locator('[name="button-OK"]').click());
    expect(settledBy,
      'Delete must announce itself through chem-mpo-profile-deleted within 6 s — GROK-19624 was fixed by ' +
      'firing that event instead of waiting for a manual refresh, so a silent delete is the regression')
      .toBe('event');
  });

  await softStep('Scenario 2 Step 9: after delete, the profile row is absent from the MPO Profiles list (and from disk) immediately', async () => {
    const onDisk = await readDiskProfiles(page);
    assertAllProfilesReadable(onDisk);
    expect(onDisk.length, 'the profiles directory read back empty — a disk absence check on it is vacuous')
      .toBeGreaterThan(0);
    expect(await profileFilesStillOnDisk(page, SYNC_PROFILE),
      'deleting through the row menu must remove the profile file from disk').toEqual([]);
    const links = await listProfileLinkTexts(page);
    // Positive control: an empty list satisfies the absence check on its own, so a list that
    // failed to render would look like a successful delete.
    expect(links.length, 'the profiles list must still be rendering other profiles for its absence check to mean anything')
      .toBeGreaterThan(0);
    expect(links, `the deleted profile must be gone from the rendered list; links: ${JSON.stringify(links)}`)
      .not.toContain(SYNC_PROFILE);
  });

  await softStep('Scenario 2 Step 10: the Browse tree returns to the baseline count — no stale entry retained (GROK-19624 regression guard)', async () => {
    const tree = await readBrowseTreeMpoProfiles(page);
    assertTreeRead(tree);
    const base = await page.evaluate(() => (window as any).__baseTreeCount as number);
    expect(tree).not.toContain(SYNC_PROFILE);
    expect(tree.length).toBe(base);
  });

  await softStep('Scenario 3 Step 11-13: score smiles.csv against a sloped MPO profile through the scoring pipeline', async () => {
    const result = await page.evaluate(async ({dir, profileName, column}) => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((res) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(res, 8000);
      });
      const col = df.col(column);
      const profile = {
        type: 'MPO Desirability Profile', name: profileName, description: '',
        properties: {[column]: {functionType: 'numerical', weight: 1, mode: 'freeform',
          min: col.min, max: col.max, line: [[col.min, 1.0], [col.max, 0.0]], rangeUserSet: true}},
        aggregation: 'Average',
      };
      await grok.dapi.files.writeAsText(`${dir}/${profileName}.json`, JSON.stringify(profile));
      // Chem:getMpoProfileNames() calls MpoProfileManager.load() unconditionally, re-reading the
      // mpo dir from disk. Without it, mpoScoreByProfile's ensureLoaded() keeps the stale singleton
      // cache from an earlier load and throws "MPO profile ... not found" (the prior Gate B failure).
      const names = await grok.functions.call('Chem:getMpoProfileNames', {});
      if (!names.includes(profileName))
        return {outName: null, columnsAfter: [], reloadMissing: true};
      const fn = DG.Func.find({package: 'Chem', name: 'mpoScoreByProfile'})[0];
      const call = await fn.prepare({table: df, profileName, columnMapping: '',
        aggregation: 'Average', createDesirabilityColumns: false}).call();
      const outName = call.getOutputParamValue()?.name ?? null;
      return {outName, reloadMissing: false, columnsAfter: df.columns.names().filter((n: string) => /^MPO /.test(n))};
    }, {dir: MPO_DIR, profileName: SCORE_PROFILE, column: SCORE_COLUMN});
    expect(result.reloadMissing, `${SCORE_PROFILE} absent from getMpoProfileNames() after write — registry did not refresh`).toBe(false);
    expect(result.outName).toBe(`MPO ${SCORE_PROFILE}`);
    expect(result.columnsAfter).toContain(`MPO ${SCORE_PROFILE}`);
  });

  await softStep('Scenario 3 Step 14: the MPO score column holds varying values in [0, 1] — not uniformly zero or one', async () => {
    const stats = await page.evaluate((profileName) => {
      const df = grok.shell.t;
      const col = df.col(`MPO ${profileName}`);
      if (!col) return null;
      const vals: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const v = col.get(i);
        if (v != null && !isNaN(v)) vals.push(v);
      }
      const distinct = new Set(vals.map((v) => v.toFixed(4)));
      return {
        count: vals.length,
        distinctCount: distinct.size,
        anyBelowOne: vals.some((v) => v < 1),
        allInUnit: vals.every((v) => v >= 0 && v <= 1),
      };
    }, SCORE_PROFILE);
    console.log(`Step 14: MPO ${SCORE_PROFILE} stats: ${JSON.stringify(stats)}`);
    expect(stats).not.toBeNull();
    expect(stats!.count).toBeGreaterThan(0);
    expect(stats!.distinctCount).toBeGreaterThan(1);
    expect(stats!.anyBelowOne).toBe(true);
    expect(stats!.allInUnit).toBe(true);
  });

  await softStep('Scenario 3 Step 15: delete ScoreTest-Profile to restore the baseline state', async () => {
    await deleteDiskProfile(page, `${SCORE_PROFILE}.json`);
    const onDisk = await readDiskProfiles(page);
    assertAllProfilesReadable(onDisk);
    expect(onDisk.length, 'the profiles directory read back empty — a cleanup check on it is vacuous')
      .toBeGreaterThan(0);
    expect(await profileFilesStillOnDisk(page, SCORE_PROFILE),
      'cleanup left the scoring profile behind').toEqual([]);
  });

  finishSpec('Chem MPO profile sync failures');
});
