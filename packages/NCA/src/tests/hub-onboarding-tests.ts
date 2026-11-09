/**
 * Hub onboarding (topic: nca-hub-onboarding) — the demo-tables entry point
 * ("Open demo dataset" link, which replaced the ★ icon 2026-10-12) and the
 * header's "Learn more" link.
 *
 * These tests pin the assumption the feature rests on that is NOT verifiable by
 * reading TypeScript, because it is resolved inside a third-party component:
 *
 *  - `AppHeader` (`@datagrok-libraries/u2`) renders "Learn more" under the
 *     bullets and runs its `onClick` in place, without navigating away: the Hub
 *     relies on that to open its in-app help panel.
 *
 * The third test pins idempotence against the real AppData path — the property
 * the user asked for by name ("if these tables are not already in the list").
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AppHeader} from '@datagrok-libraries/u2';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {ensureDemoTables, loadSampleByName} from '../views/nca-launch';
import {
  DEMO_TABLE_FILES, DEMO_ID_TAG, PROVENANCE_TAG,
} from '../views/nca-demo-tables';
import {
  HUB_HELP_MD, INFO_ICON_TOKEN, ADD_ICON_TOKEN, SHARE_ICON_TOKEN,
  buildHubHelpElement, openHubHelp, restoreHubHelpSync,
  DEMO_BALLOON_MD,
} from '../views/nca-hub-help';

/** Demo table DISPLAY names — decoupled from the filename (user decision D2/D5:
 *  compound + route + species), so they are read from the entry, never derived
 *  from the file. */
const DEMO_NAMES = DEMO_TABLE_FILES.map((e) => e.name);

/** Drop any demo tables this category opened, so a re-run starts cold and the
 *  idempotence test is testing the code rather than a leftover from itself.
 *  Matches on the stable id tag AS WELL AS the name — the rename test leaves a
 *  demo frame whose name no longer matches anything, and a name-only sweep
 *  would leak it into the next test. */
function closeDemoTables(): void {
  const files = new Set(DEMO_TABLE_FILES.map((e) => e.file));
  const names = new Set(DEMO_NAMES);
  const doomed = grok.shell.tables.filter((t: DG.DataFrame) => {
    const id = t.getTag(DEMO_ID_TAG) as string | null | undefined;
    return (id != null && id !== '' && files.has(id)) || names.has(t.name);
  });
  for (const t of doomed)
    grok.shell.closeTable(t);
}

category('Hub onboarding: demo tables (★)', () => {
  test('loads both demo tables through the standard AppData path', async () => {
    closeDemoTables();
    try {
      const tables = await ensureDemoTables();
      expect(tables.length, DEMO_TABLE_FILES.length, 'one table per demo file');
      expect(tables.map((t) => t.name).join(','), DEMO_NAMES.join(','),
        'returned in DEMO_TABLE_FILES order — the Hub seeds its picker with [0]');
      // Non-empty is the real assertion: a wrong AppData path resolves to an
      // error, but a mis-parsed CSV would hand back a hollow frame.
      for (const t of tables) {
        expect(t.rowCount > 0, true, `${t.name} has rows`);
        expect(t.columns.length > 0, true, `${t.name} has columns`);
      }
      // Both demo files carry SMILES — that is why these two were chosen over
      // the plain fixtures, so pin it rather than trusting the filename.
      // (The stronger property — that the column is TYPED `Molecule` — is not
      // assertable in this harness; see the F0 test below for why.)
      for (const t of tables) {
        const hasSmiles = Array.from(t.columns.names() as string[])
          .some((n: string) => /smiles/i.test(n));
        expect(hasSmiles, true, `${t.name} carries a SMILES column`);
      }
      for (const name of DEMO_NAMES) {
        expect(grok.shell.tables.some((t: DG.DataFrame) => t.name === name), true,
          `${name} is in the platform table list`);
      }
    } finally {
      closeDemoTables();
    }
  });

  test('addTable does NOT auto-detect semantic types — the loader must do it (F0)', async () => {
    // WHY THIS TEST EXISTS. `ensureDemoTables` calls `grok.data.detectSemanticTypes`
    // so the demo's SMILES column is typed `Molecule`; without typing,
    // nca-add-to-study.ts:206 reads the LIVE semType once and persists
    // `isStructure: false` permanently, so the demo's structures would never
    // render in a study grid. That call is only justified if the platform does
    // not already do it. Measured here rather than assumed — the API surface
    // documents neither behaviour, and inferring platform behaviour from types
    // is exactly the failure mode this test was written to avoid.
    closeDemoTables();
    const df = await loadSampleByName(DEMO_TABLE_FILES[0].file);
    try {
      const smiles = Array.from(df.columns.names() as string[])
        .find((n: string) => /smiles/i.test(n))!;
      expect(smiles !== undefined, true, 'the fixture has a SMILES column to type');

      // AC-F1b — `loadSampleByName`'s docstring promises "the filename stem
      // (e.g. 'theoph')", and the demo entries are now NESTED paths. A plain
      // extension-strip would yield 'demo/indomethacin_iv' and quietly break
      // that contract for any future caller. Asserted here rather than in its
      // own test because the frame is already loaded — no extra cost.
      // NOTE: this cannot be observed through `ensureDemoTables`, which
      // overwrites `df.name` with the display name immediately after load; only
      // a direct `loadSampleByName` call can see the stem.
      expect(df.name, 'indomethacin_iv',
        'nested path -> directory stripped as well as extension');
      const flat = await loadSampleByName('theoph.csv');
      expect(flat.name, 'theoph', 'flat path unchanged — the fix is behaviour-neutral');

      // Observed 2026-07-28 against a live local server: unset (`undefined`).
      expect(df.col(smiles)!.semType === DG.SEMTYPE.MOLECULE, false,
        'a freshly parsed CSV carries no semantic type');

      grok.shell.addTable(df);
      // Observed 2026-07-28: STILL unset. This is the load-bearing assertion —
      // it says the platform will not type the column for us, which is what
      // makes the explicit detectSemanticTypes call in `ensureDemoTables`
      // necessary rather than redundant. If this ever fails, the platform
      // started auto-detecting and that call can be removed.
      expect(df.col(smiles)!.semType === DG.SEMTYPE.MOLECULE, false,
        'grok.shell.addTable performs NO implicit semantic-type detection — ' +
        'so ensureDemoTables must call detectSemanticTypes explicitly');

      // THE HARNESS BOUNDARY, asserted so it cannot be mistaken for a result.
      // The complementary half of F0 — "does an explicit detectSemanticTypes
      // call actually yield Molecule for OUR column?" — is NOT answerable here:
      // `grok test` loads only the package under test, so no package-contributed
      // detector exists in this context (Chem's detectMolecule is unreachable:
      // "Unable to get project asset"). Detection therefore runs over an EMPTY
      // detector set and returns null for every column, including a canonical
      // SMILES frame. That null is a property of the harness, not of the data,
      // and it is why AC-F1a is discharged by live visual QA instead.
      // This assertion is a tripwire: if the harness ever gains detectors it
      // fails, and the fix is to upgrade the test to assert Molecule directly.
      expect(DG.Func.find({tags: ['semTypeDetector']}).length, 0,
        'no semantic-type detector is loaded in the package-test harness — if ' +
        'this fails the harness gained detectors: assert semType === MOLECULE here');
    } finally {
      closeDemoTables();
    }
  });

  test('is idempotent — a second click reuses, never duplicates', async () => {
    closeDemoTables();
    try {
      const first = await ensureDemoTables();
      const countAfterFirst = grok.shell.tables.length;
      const second = await ensureDemoTables();
      expect(grok.shell.tables.length, countAfterFirst,
        'a repeat ★ click adds no table');
      // Reused, not reloaded: the same DataFrame instances come back, which is
      // what keeps a picker already pointed at a demo table valid.
      for (let i = 0; i < first.length; i++)
        expect(second[i] === first[i], true, `${DEMO_NAMES[i]} instance reused`);
      for (const name of DEMO_NAMES) {
        expect(grok.shell.tables.filter((t: DG.DataFrame) => t.name === name).length, 1,
          `exactly one ${name} open`);
      }
    } finally {
      closeDemoTables();
    }
  });

  test('idempotence survives a user rename (AC-F3a)', async () => {
    // The old key was `df.name`, which the user can change at will — renaming a
    // demo table then clicking ★ again added a duplicate. The key is now the
    // stable file id carried on a frame tag.
    closeDemoTables();
    try {
      const first = await ensureDemoTables();
      const countAfterFirst = grok.shell.tables.length;
      first[0].name = 'My renamed demo table';

      const second = await ensureDemoTables();
      expect(grok.shell.tables.length, countAfterFirst,
        'a ★ click after a rename adds no table');
      expect(second[0] === first[0], true,
        'the RENAMED table is the one reused — resolution is by id, not by name');
      expect(second[0].name, 'My renamed demo table',
        'reuse does not silently rename the user table back');
    } finally {
      closeDemoTables();
    }
  });

  test('a pre-tag frame is adopted by name, and the tag is backfilled (AC-F3b/AC-F3d)', async () => {
    // A demo table opened BEFORE the id tag existed carries the right name but
    // no tag. It must be adopted rather than duplicated — and the tag must be
    // stamped ON ADOPTION, or the frame re-enters this same branch on every
    // later click and re-runs the whole-frame detection scan. A table-count
    // assertion cannot see that: the count stays right while the work repeats.
    closeDemoTables();
    try {
      const entry = DEMO_TABLE_FILES[0];
      const legacy = await loadSampleByName(entry.file);
      legacy.name = entry.name; // named as the demo, deliberately untagged
      grok.shell.addTable(legacy);
      const before = grok.shell.tables.length;

      const first = await ensureDemoTables();
      expect(first[0] === legacy, true, 'the untagged look-alike is adopted, not duplicated');
      expect(grok.shell.tables.filter((t: DG.DataFrame) => t.name === entry.name).length, 1,
        'exactly one table under the demo name');

      const tag = legacy.getTag(DEMO_ID_TAG) as string | null | undefined;
      expect(tag === null || tag === undefined || tag === '', false,
        'adoption BACKFILLS the id tag, so the next click takes the tag branch');
      expect(tag, entry.file, 'the backfilled tag is this entry file id');

      const second = await ensureDemoTables();
      expect(second[0] === legacy, true, 'the second click reuses the same frame');
      expect(grok.shell.tables.length, before + DEMO_TABLE_FILES.length - 1,
        'no table was added by either call beyond the other demo entry');
    } finally {
      closeDemoTables();
    }
  });

  test('both frames carry a provenance tag, and no column is added (AC-F4a/AC-F4b)', async () => {
    closeDemoTables();
    try {
      const tables = await ensureDemoTables();
      for (let i = 0; i < tables.length; i++) {
        const entry = DEMO_TABLE_FILES[i];
        const tag = tables[i].getTag(PROVENANCE_TAG) as string | null | undefined;
        expect(tag === null || tag === undefined || tag === '', false,
          `${entry.name} carries a provenance tag`);
        expect(tag, entry.provenance, `${entry.name} provenance text is the entry's`);
      }
      // Positive provenance on BOTH frames, not only the synthetic one: a tag
      // present only on the simulated table would make its ABSENCE the signal.
      const real = tables[0].getTag(PROVENANCE_TAG) ?? '';
      const sim = tables[1].getTag(PROVENANCE_TAG) ?? '';
      expect(/real data/i.test(real), true, 'the indomethacin frame states its real lineage');
      expect(/simulated/i.test(sim), true, 'the multiroute frame states the numbers are simulated');
      expect(/SEED=\d+/.test(sim), true, 'the simulated frame names its generator seed');

      // AC-F4b — tagging is frame-level metadata and must not touch the column
      // set, which is what the mapping pane's pickers read. Compared against a
      // clean re-load of the same fixture rather than a remembered list.
      for (let i = 0; i < tables.length; i++) {
        const clean = await loadSampleByName(DEMO_TABLE_FILES[i].file);
        expect(Array.from(tables[i].columns.names() as string[]).join(','),
          Array.from(clean.columns.names() as string[]).join(','),
          `${DEMO_TABLE_FILES[i].name} column set is byte-identical to the fixture`);
      }
    } finally {
      closeDemoTables();
    }
  });
});

/**
 * These pin the third-party SEAM the Hub's onboarding affordance stands on,
 * rather than the Hub view itself. (The `addOptions` seam test left with the ★
 * icon — demo data is now an "Open demo dataset" link, 2026-10-12.)
 *
 * **Why not mount `NcaAppView`?** It cannot run here. `buildHeader` reads
 * `_package.webRoot`, and inside the test bundle `src/package.ts`'s `_package`
 * is a bare `new DG.Package()` that the platform never binds — `package-test.ts`
 * declares its OWN `_package` (line 43), and that is the one the platform
 * initializes. `grok_Package_Get_WebRoot` throws on the unbound handle. No
 * existing test mounts a `_package`-touching view for exactly this reason; this
 * is a harness boundary, not something to work around by contorting the view.
 *
 * So the tests target what actually carries risk. The Hub's own wiring is one
 * readable call into this seam; the seam is a third-party contract that can
 * change under us with no compile error. If it breaks, Learn more silently stops
 * opening the help panel — a failure no type check catches.
 */
category('Hub onboarding: seams the Hub stands on', () => {
  test('AppHeader renders Learn more under the bullets and runs its onClick in place', async () => {
    let clicks = 0;
    const header = new AppHeader({icon: 'flask', description: '- bullet one\n- bullet two\n',
      learnMore: {onClick: () => clicks++}});
    try {
      const bullets = header.root.querySelectorAll('.u2-app-header-item');
      expect(bullets.length, 2, 'one item per markdown bullet');
      const more = header.root.querySelector('.u2-app-header-learn-more') as HTMLAnchorElement | null;
      expect(more !== null, true, 'Learn more is rendered');
      expect(more?.textContent, 'Learn more', 'the default link text');
      expect(more !== null && (bullets[1].compareDocumentPosition(more) & Node.DOCUMENT_POSITION_FOLLOWING) !== 0,
        true, 'Learn more follows the bullets');
      expect(more?.getAttribute('href') ?? null, null, 'an onClick link does not navigate away');
      more?.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      expect(clicks, 1, 'the click runs the onClick (the Hub opens its help panel there)');
    } finally {
      header.dispose();
    }
  });
});

category('Hub onboarding: demo balloon (★)', () => {
  test('renders as a markdown list, one item per demo table (AC-F5a)', async () => {
    const balloonMd = DEMO_BALLOON_MD;
    const el = ui.markdown(balloonMd);
    const items = el.querySelectorAll('li');
    expect(items.length, DEMO_TABLE_FILES.length,
      'one <li> per demo table — a plain string balloon would render no list at all');
    expect(el.querySelector('ul') !== null, true, 'the bullets are a real <ul>');
  });

  test('every table name is derived from DEMO_TABLE_FILES, never typed twice (AC-F5b)', async () => {
    const balloonMd = DEMO_BALLOON_MD;
    for (const e of DEMO_TABLE_FILES) {
      expect(balloonMd.includes(e.name), true,
        `the balloon names ${e.name} — and gets it from the entry, so a rename cannot drift`);
      expect(balloonMd.includes(e.unlocks), true,
        `the balloon says what ${e.name} unlocks`);
    }
    // The species qualifier is the point of D5 — assert it survived into the copy.
    expect(/\(human\)/.test(balloonMd) && /\(rat\)/.test(balloonMd), true,
      'both display names carry their species');
  });

  test('the simulated-data caveat is present, italic, and within its 4-word budget', async () => {
    // The ★ balloon is the ONLY user-visible disclosure that the multiroute
    // numbers are simulated: AC-F4c measured that the frame provenance tag
    // renders in no platform surface whatsoever. So this test guards the
    // caveat's EXISTENCE, not just its wording.
    //
    // It also guards the user's 4-word budget (2026-07-28): the first version
    // was a full sentence plus a meta-clause about its own persistence, which
    // the user reviewed on screen and asked to be cut to one short italic line.
    // Both bounds matter — losing the caveat breaks disclosure, and letting it
    // grow back re-creates the copy the user rejected.
    const balloonMd = DEMO_BALLOON_MD;
    const withCaveat = DEMO_TABLE_FILES.filter((e) => e.qualifier !== undefined);
    expect(withCaveat.length > 0, true,
      'at least one demo table is simulated and must say so');

    for (const e of withCaveat) {
      expect(balloonMd.includes(`*${e.qualifier}*`), true,
        `${e.name} carries its caveat, and carries it ITALIC (single asterisks)`);
      const words = e.qualifier!.trim().split(/\s+/).length;
      expect(words <= 4, true,
        `caveat is ${words} words — the budget is 4: "${e.qualifier}"`);
      // Spend the budget on the two facts that matter, not on hedging.
      expect(/simulat/i.test(e.qualifier!), true, 'the caveat names the data as simulated');
      expect(/illustrativ/i.test(e.qualifier!), true, 'the caveat qualifies the structure');
    }
    // A table with nothing to disclose must NOT grow a caveat — absence has to
    // keep meaning "nothing to disclose", or the signal is worthless.
    for (const e of DEMO_TABLE_FILES.filter((x) => x.qualifier === undefined))
      expect(balloonMd.includes(`${e.unlocks} *`), false, `${e.name} has no caveat`);

    // The machine-readable half must survive alongside it — the balloon is
    // dismissible, the tag is not.
    closeDemoTables();
    try {
      const tables = await ensureDemoTables();
      const sim = tables[1].getTag(PROVENANCE_TAG) ?? '';
      expect(/simulated/i.test(sim), true,
        'the frame tag still states the numbers are simulated');
      expect(/illustrative/i.test(sim), true, 'the frame tag still qualifies the structure');
    } finally {
      closeDemoTables();
    }
  });

  test('uses no glyph the balloon font cannot render (AC-F5d)', async () => {
    const balloonMd = DEMO_BALLOON_MD;
    // Same ban as the help panel: a fullwidth ＋ shipped invisibly once
    // (2026-07-17) and rendered as a blank gap — no error, no tofu, just a
    // missing word. Nothing but a screenshot catches it, so ban the family.
    // Escape sequences, not literal characters: spelling these ranges out
    // literally would put the very codepoints under test into this source file.
    const fullwidth = balloonMd.match(/[\uFF00-\uFFEF]/g);
    expect(fullwidth, null,
      `fullwidth glyphs do not render: ${JSON.stringify(fullwidth)}`);
    const allowed = new Set(['★', '▶', '∞', '→', '≥']);
    const exotic = (balloonMd.match(/[\u2190-\u2BFF\u{1F000}-\u{1FAFF}]/gu) ?? [])
      .filter((c) => !allowed.has(c));
    expect(exotic.length, 0,
      `unverified pictograph(s) ${JSON.stringify(exotic)} — confirm they render, then allowlist`);
  });
});

category('Hub onboarding: help content', () => {
  test('carries the three briefed sections and quotes real control labels', async () => {
    expect(HUB_HELP_MD.includes('## Quick start'), true, 'Quick start section');
    expect(HUB_HELP_MD.includes('## Computations'), true, 'Computations section (placeholder)');
    // The entity model is carried by the steps now (save into a study, compare
    // on the Studies view), not a bold-term intro — so assert the nouns are
    // present, not a specific markup.
    expect(/non-compartmental analysis/i.test(HUB_HELP_MD), true, 'intro states what the app is');
    expect(/analys/i.test(HUB_HELP_MD) && /stud/i.test(HUB_HELP_MD), true,
      'names both the analysis unit and the study collection');
    // Copy-drift guard: these map to live tooltips in nca-mapping-pane.ts. If a
    // tooltip is reworded there and not here, the help starts naming a control
    // the user cannot find.
    expect(HUB_HELP_MD.includes('Confirm the setup and compute'), true,
      'compute step quotes the Confirm tooltip (the ribbon ▶ icon was removed, D17)');
    expect(HUB_HELP_MD.includes('Run interactive guide'), true,
      'points at the existing in-view tour');
    // All three control icons are spliced-in elements, not glyphs — their
    // anchor tokens must be present or the icons have nowhere to land.
    expect(HUB_HELP_MD.includes(INFO_ICON_TOKEN), true,
      'the %INFO_ICON% splice anchor is present in the copy');
    expect(HUB_HELP_MD.includes(ADD_ICON_TOKEN), true,
      'the %ADD_ICON% splice anchor is present in the copy');
    expect(HUB_HELP_MD.includes(SHARE_ICON_TOKEN), true,
      'the %SHARE_ICON% splice anchor is present in the copy');
  });

  test('renders ALL three real app icons and leaves no token behind', async () => {
    // The user asked to show the actual info, add-to-study, and share icons, not
    // text or emoji (2026-07-17). This proves the splice puts the real
    // FontAwesome nodes — the same glyphs the app buttons use — into the
    // rendered help, with none of the %…% anchors left behind. A font-dependent
    // emoji could not be asserted this way; a real element can.
    const el = buildHubHelpElement();
    expect(el.querySelector('i.fa-info-circle') !== null, true,
      'the real info-circle icon node was spliced in');
    expect(el.querySelector('i.fa-layer-plus') !== null, true,
      'the real layer-plus icon node was spliced in');
    expect(el.querySelector('i.fa-share') !== null, true,
      'the real share icon node was spliced in');
    const text = el.textContent ?? '';
    const leaked = [INFO_ICON_TOKEN, ADD_ICON_TOKEN, SHARE_ICON_TOKEN]
      .some((tok) => text.includes(tok));
    expect(leaked, false, 'no raw %…% token leaked into the rendered text');
  });

  test('uses no glyph the help panel font cannot render', async () => {
    // A fullwidth ＋ (U+FF0B) shipped in the add-to-study step and rendered as
    // a blank gap — "Press the  layers icon" (user-reported 2026-07-17). The
    // failure is silent: no error, no tofu box, just a missing word, so nothing
    // but a live screenshot catches it. Ban the family that caused it.
    const fullwidth = HUB_HELP_MD.match(/[\uFF00-\uFFEF]/g);
    expect(fullwidth, null,
      `fullwidth glyphs do not render in the help panel: ${JSON.stringify(fullwidth)} — ` +
      'name the control by its tooltip instead');
    // The ranges below are escape sequences on purpose: spelling them with
    // literal characters puts invisible codepoints into the source — the very
    // bug under test, and something no-irregular-whitespace rejects outright.
    // ★ (U+2605), ▶ (U+25B6), and the math symbols ∞ → ≥ are text-presentation
    // characters that render reliably (unlike emoji). Everything else in the
    // swept ranges (symbols / arrows / dingbats / emoji) is a gamble against
    // font coverage — confirm on screen first, then allowlist deliberately.
    // Typographic punctuation and sub/superscripts sit below U+2190 and are not
    // swept up.
    const allowed = new Set(['★', '▶', '∞', '→', '≥']);
    const exotic = (HUB_HELP_MD.match(/[\u2190-\u2BFF\u{1F000}-\u{1FAFF}]/gu) ?? [])
      .filter((c) => !allowed.has(c));
    expect(exotic.length, 0,
      `unverified pictograph(s) ${JSON.stringify(exotic)} — confirm they render, then allowlist`);
  });

  test('openHubHelp disables current-object sync; restoreHubHelpSync scopes it to the Hub', async () => {
    // The first click on a CLOSED help panel showed the default package help,
    // not the NCA help (user-reported 2026-07-17). Root cause: the panel opens
    // with syncCurrentObject = true and the open-time sync to grok.shell.o
    // clobbers showHelp()'s content. The fix switches sync OFF before showHelp.
    // But syncCurrentObject is a GLOBAL platform singleton — leaving it off
    // forever would break contextual help for every other package (decision-
    // audit finding). So the disable is SCOPED: restoreHubHelpSync puts the
    // prior value back when the Hub stops being the current view. Pin both
    // halves so a refactor cannot drop either (bring back the two-click bug, or
    // leak the global disable session-wide). Simulate the pre-first-click state
    // (sync ON) so the assertions prove openHubHelp is what turns it off.
    const panel = grok.shell.windows.help;
    const prevVisible = panel.visible;
    const prevSync = panel.syncCurrentObject;
    try {
      panel.syncCurrentObject = true;
      openHubHelp();
      expect(panel.syncCurrentObject, false,
        'openHubHelp disables syncCurrentObject before handing the panel its content');
      expect(panel.visible, true, 'the help panel is shown');
      restoreHubHelpSync();
      expect(panel.syncCurrentObject, true,
        'restoreHubHelpSync restores the prior value — the disable is scoped to the Hub, not the session');
    } finally {
      panel.syncCurrentObject = prevSync;
      panel.visible = prevVisible;
    }
  });
});
