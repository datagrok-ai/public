/**
 * NCA Studies — URL router category.
 *
 * @see src/views/nca-studies-router.ts (under test)
 *
 * Coverage:
 * - `parseStudyPath` maps ''/`study`/`study,analysis` to the right state; tolerates
 *   trailing slashes + empty segments; decodes `encodeURIComponent`-escaped label
 *   segments; the query-param fallback (AC-A3) promotes a depth-1 path + `analysisId`
 *   to the analysis state via the SAME state machine.
 * - `updateNcaStudyPath` writes the app-RELATIVE tail `/<study label>/<analysis label>`
 *   to `view.path` (the platform prepends `view.basePath`), encoding labels, degrading
 *   correctly on missing labels, and giving a view without THIS app's root `/apps/Nca`.
 * - encode(write) → decode(parse) round-trips a label verbatim, including spaces and
 *   a literal `/` (which must stay inside a single segment, not split the path).
 * - `resolveDeepLink` (SD-26): two visible studies of one label land on the Hub with a toast; the
 *   stale-label and moved-standalone branches resolve as the router documents.
 * - live: the real address bar carries ONE app prefix (nca-url-double-app-prefix).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, delay} from '@datagrok-libraries/test/src/test';
import {uniqueLabel} from '../adapter';
import type {StudyHeader} from '../adapter';
import {asNcaNodeMeta} from '../views/nca-tree-browser';

import {
  parseStudyPath, updateNcaStudyPath, updateNcaAnalysisPath, ncaDeepLinkPath, ncaStandaloneLinkPath, NCA_APP_ROOT,
  resolveDeepLink,
} from '../views/nca-studies-router';

category('NCA Studies: router', () => {
  test('parseStudyPath: empty → landing', async () => {
    expect(parseStudyPath('').kind, 'landing');
    expect(parseStudyPath(null).kind, 'landing');
    expect(parseStudyPath('/').kind, 'landing');
  });

  test('parseStudyPath: single segment → overview', async () => {
    const s = parseStudyPath('proj-1');
    expect(s.kind, 'overview');
    expect(s.kind === 'overview' ? s.studyLabel : '', 'proj-1');
  });

  test('parseStudyPath: two segments → analysis', async () => {
    const s = parseStudyPath('proj-1/an-2');
    expect(s.kind, 'analysis');
    expect(s.kind === 'analysis' ? s.studyLabel : '', 'proj-1');
    expect(s.kind === 'analysis' ? s.analysisLabel : '', 'an-2');
  });

  test('parseStudyPath: trailing slash + empty segments tolerated', async () => {
    const s = parseStudyPath('/proj-1/an-2/');
    expect(s.kind, 'analysis');
    expect(s.kind === 'analysis' ? s.analysisLabel : '', 'an-2');
  });

  test('parseStudyPath: decodes percent-encoded label segments', async () => {
    const s = parseStudyPath('Study%20A/Analysis%201');
    expect(s.kind, 'analysis');
    expect(s.kind === 'analysis' ? s.studyLabel : '', 'Study A');
    expect(s.kind === 'analysis' ? s.analysisLabel : '', 'Analysis 1');
  });

  test('parseStudyPath: query-param fallback promotes depth-1 to analysis (AC-A3)', async () => {
    const s = parseStudyPath('proj-1', 'an-9');
    expect(s.kind, 'analysis');
    expect(s.kind === 'analysis' ? s.studyLabel : '', 'proj-1');
    expect(s.kind === 'analysis' ? s.analysisLabel : '', 'an-9');
  });

  test('parseStudyPath: empty query param leaves depth-1 as overview', async () => {
    expect(parseStudyPath('proj-1', '').kind, 'overview');
    expect(parseStudyPath('proj-1', null).kind, 'overview');
  });

  test('updateNcaStudyPath writes the app-relative /<study>/<analysis> tail and degrades', async () => {
    const stub = {path: '', basePath: NCA_APP_ROOT} as unknown as DG.ViewBase;
    updateNcaStudyPath(stub, 'p1', 'a1');
    expect(stub.path, '/p1/a1');
    updateNcaStudyPath(stub, 'p1');
    expect(stub.path, '/p1');
    updateNcaStudyPath(stub);
    expect(stub.path, '');
    // An analysis label with no study label is meaningless → drops to study-less.
    updateNcaStudyPath(stub, null, 'a1');
    expect(stub.path, '');
    updateNcaAnalysisPath(stub, 'My NCA');
    expect(stub.path, '/analyses/My%20NCA');
  });

  test('path writers give a view without THIS app root /apps/Nca, keep an NCA root', async () => {
    const fresh = {path: '', basePath: null} as unknown as DG.ViewBase;
    updateNcaAnalysisPath(fresh, 'X');
    expect(fresh.basePath, NCA_APP_ROOT, 'null root → the app root');
    const browse = {path: '', basePath: `/browse${NCA_APP_ROOT}`} as unknown as DG.ViewBase;
    updateNcaStudyPath(browse, 'S');
    expect(browse.basePath, `/browse${NCA_APP_ROOT}`, 'Browse-mode root kept');
    const foreign = {path: '', basePath: '/apps/Chem'} as unknown as DG.ViewBase;
    updateNcaStudyPath(foreign, 'S');
    expect(foreign.basePath, NCA_APP_ROOT, 'another app root → the NCA root');
  });

  test('a table view: the root rides in path when the write drops it, else the tail alone', async () => {
    // Two measured platform behaviours of a TABLE view's `path` setter: BUG-20 (reset basePath to
    // '') and an app-opened view (basePath kept as the app root, prepended — 2026-11-02).
    const tableView = (resets: boolean): DG.ViewBase => {
      const v = {type: DG.VIEW_TYPE.TABLE_VIEW, basePath: resets ? NCA_APP_ROOT : '', raw: ''};
      return Object.defineProperty(v, 'path', {
        get: () => `${v.basePath}${v.raw}`,
        set: (p: string) => {v.raw = p; v.basePath = resets ? '' : NCA_APP_ROOT;},
      }) as unknown as DG.ViewBase;
    };
    const reset = tableView(true);
    updateNcaAnalysisPath(reset, 'X');
    expect(reset.path, `${NCA_APP_ROOT}/analyses/X`, 'BUG-20: the root rides in path');
    const appOpened = tableView(false);
    updateNcaAnalysisPath(appOpened, 'X');
    expect(appOpened.path, `${NCA_APP_ROOT}/analyses/X`, 'an app-opened view: one prefix, not two');
  });

  test('updateNcaStudyPath encodes label segments (spaces / reserved chars)', async () => {
    const stub = {path: '', basePath: NCA_APP_ROOT} as unknown as DG.ViewBase;
    updateNcaStudyPath(stub, 'Study A', 'Analysis 1');
    expect(stub.path, '/Study%20A/Analysis%201');
  });

  test('encode(write) → decode(parse) round-trips a label verbatim, incl. a slash', async () => {
    const stub = {path: '', basePath: NCA_APP_ROOT} as unknown as DG.ViewBase;
    const studyLabel = 'PK/TK Study A';
    const analysisLabel = 'IV bolus (0.5 mg/kg)';
    updateNcaStudyPath(stub, studyLabel, analysisLabel);
    // view.path IS the metaUrl tail shape (leading '/' tolerated by parseStudyPath).
    const tail = stub.path;
    const s = parseStudyPath(tail);
    expect(s.kind, 'analysis');
    expect(s.kind === 'analysis' ? s.studyLabel : '', studyLabel);
    expect(s.kind === 'analysis' ? s.analysisLabel : '', analysisLabel);
  });

  test('ncaDeepLinkPath: label-keyed, encoded, round-trips through the router (copy-link contract)', async () => {
    // The exact shape copy-link buttons build — must match what routeNcaApp resolves.
    expect(ncaDeepLinkPath('New study-1', 'fourth'), '/apps/Nca/New%20study-1/fourth');
    expect(ncaDeepLinkPath('PK/TK study'), '/apps/Nca/PK%2FTK%20study'); // slash stays ONE segment
    expect(ncaDeepLinkPath(), '/apps/Nca');
    // A pasted copy-link tail must parse back to the SAME labels — guards the
    // uuid→label + unencoded-slash contract-drift that shipped a broken analysis link.
    const tail = ncaDeepLinkPath('PK/TK study', 'IV bolus').slice(NCA_APP_ROOT.length);
    const s = parseStudyPath(tail);
    expect(s.kind, 'analysis');
    expect(s.kind === 'analysis' ? s.studyLabel : '', 'PK/TK study');
    expect(s.kind === 'analysis' ? s.analysisLabel : '', 'IV bolus');
  });
});

/** A visible study / standalone container as `listStudies` returns it (only what the resolve reads). */
function header(id: string, label: string, kind: 'study' | 'standalone',
  analyses: Array<{id: string; label: string; movedFromLabel?: string}>): StudyHeader {
  return {id, label, kind, updatedAt: '', analysisCount: analyses.length,
    access: {canEdit: true, canDelete: true, canShare: true},
    analyses: analyses.map((a) => ({id: a.id, studyId: id, label: a.label, engineVersion: null, resultsBasis: null,
      hasCache: false,
      stageCache: null, pipelineFlagCount: 0, movedFromLabel: a.movedFromLabel ?? null,
      updatedAt: ''}))};
}

// SD-26: the resolve step, apart from the Hub it falls back to (NcaAppView cannot render under grok test).
category('NCA Studies: deep-link resolve', () => {
  const twins = [header('s1', 'Rat PK', 'study', [{id: 'a1', label: 'IV'}]),
    header('s2', 'Rat PK', 'study', [{id: 'a2', label: 'IV'}])];

  test('SD-26: two visible studies share the label → the Hub with a toast, never a guess', async () => {
    for (const path of ['Rat%20PK', 'Rat%20PK/IV']) {
      const t = resolveDeepLink(parseStudyPath(path), twins);
      expect(t.kind, 'landing', path);
      expect(t.kind === 'landing' ? t.notice : '', 'More than one study is named \'Rat PK\'. Open it from the Hub.');
    }
  });

  test('one study: the overview, the analysis, a stale analysis → overview, a stale study → Hub', async () => {
    const one = [twins[0]];
    expect(JSON.stringify(resolveDeepLink(parseStudyPath('Rat%20PK'), one)), '{"kind":"overview","studyId":"s1"}');
    expect(JSON.stringify(resolveDeepLink(parseStudyPath('Rat%20PK/IV'), one)),
      '{"kind":"analysis","analysisId":"a1"}');
    expect(JSON.stringify(resolveDeepLink(parseStudyPath('Rat%20PK/PO'), one)), '{"kind":"overview","studyId":"s1"}');
    expect(JSON.stringify(resolveDeepLink(parseStudyPath('Dog%20PK'), one)), '{"kind":"landing"}');
  });

  test('analyses/<label>: own → open; moved once → there; moved twice → Hub with a toast; gone → Hub', async () => {
    const own = [header('u1', 'Scratch', 'standalone', [{id: 'x1', label: 'Scratch'}])];
    expect(JSON.stringify(resolveDeepLink(parseStudyPath('analyses/Scratch'), own)),
      '{"kind":"analysis","analysisId":"x1"}');
    const movedOnce = [header('s1', 'Rat PK', 'study', [{id: 'm1', label: 'IV', movedFromLabel: 'Scratch'}])];
    expect(JSON.stringify(resolveDeepLink(parseStudyPath('analyses/Scratch'), movedOnce)),
      '{"kind":"analysis","analysisId":"m1"}');
    const movedTwice = [...movedOnce,
      header('s2', 'Dog PK', 'study', [{id: 'm2', label: 'IV', movedFromLabel: 'Scratch'}])];
    const t = resolveDeepLink(parseStudyPath('analyses/Scratch'), movedTwice);
    expect(t.kind === 'landing' ? t.notice : '', '\'Scratch\' was moved into 2 studies — open it from the hub.');
    expect(JSON.stringify(resolveDeepLink(parseStudyPath('analyses/Gone'), own)), '{"kind":"landing"}');
  });
});

// nca-url-double-app-prefix: the platform builds the address bar as `view.basePath + view.path`
// (plain concatenation), so the stub tests above cannot see a doubled prefix — these read the REAL
// URL. A view opened inside the app carries basePath `/apps/Nca`; writing the absolute
// `/apps/NCA/…` into `view.path` produced `/apps/Nca/apps/NCA/…` (user QA 2026-10-06).
category('NCA Studies: router (live URL)', () => {
  /** Open a throw-away view with `basePath` (null = leave the platform default), apply `write`,
   *  make it current and return the address-bar pathname. */
  async function urlAfter(basePath: string | null, write: (v: DG.View) => void): Promise<string> {
    const v = DG.View.create();
    grok.shell.addView(v);
    try {
      if (basePath !== null) v.basePath = basePath;
      write(v);
      grok.shell.v = v;
      await delay(300);
      return window.location.pathname;
    } finally {
      v.close();
    }
  }

  test('NCA_APP_ROOT is /apps/ + the app package name (case-exact)', async () => {
    const app = DG.Func.find({tags: ['app']}).find((f) => f.package?.name?.toLowerCase() === 'nca');
    expect(app !== undefined, true, 'the NCA app function is registered');
    expect(NCA_APP_ROOT, `/apps/${app!.package.name}`);
  });

  test('app-launched standalone analysis: ONE app prefix, equal to its copy-link', async () => {
    const label = 'Indomethacin IV (human)';
    const url = await urlAfter(NCA_APP_ROOT, (v) => updateNcaAnalysisPath(v, label));
    expect(url, `${NCA_APP_ROOT}/analyses/Indomethacin%20IV%20(human)`);
    expect(url, ncaStandaloneLinkPath(label), 'the copy-link is the address-bar URL');
  });

  test('app-launched study analysis: ONE app prefix, equal to its copy-link', async () => {
    const url = await urlAfter(NCA_APP_ROOT,
      (v) => updateNcaStudyPath(v, 'Indomethacin IV (human)', 'NCA · PO 1 mg'));
    expect(url, `${NCA_APP_ROOT}/Indomethacin%20IV%20(human)/NCA%20%C2%B7%20PO%201%20mg`);
    expect(url, ncaDeepLinkPath('Indomethacin IV (human)', 'NCA · PO 1 mg'), 'the copy-link is the address-bar URL');
  });

  test('a view without the app root (fresh View / preview) still gets /apps/Nca', async () => {
    expect(await urlAfter(null, (v) => updateNcaAnalysisPath(v, 'X')), `${NCA_APP_ROOT}/analyses/X`);
    expect(await urlAfter(null, (v) => updateNcaStudyPath(v, 'S')), `${NCA_APP_ROOT}/S`);
  });

  // BUG-20: an analysis is a TABLE view, whose `path` setter resets `basePath` (core
  // table_view.dart) — the stub-free View tests above could not see the lost prefix.
  test('a table view (every analysis) keeps /apps/Nca, also on a second write', async () => {
    const v = grok.shell.addTableView(grok.data.demo.demog(5));
    try {
      updateNcaAnalysisPath(v, 'Indomethacin IV (human)');
      grok.shell.v = v;
      await delay(300);
      expect(window.location.pathname, ncaStandaloneLinkPath('Indomethacin IV (human)'));
      updateNcaStudyPath(v, 'S', 'A'); // a later bind (save / add to study) rewrites the path
      await delay(300);
      expect(window.location.pathname, ncaDeepLinkPath('S', 'A'));
    } finally {
      v.close();
    }
  });

  test('a Browse-mode root is kept (/browse/apps/Nca/…)', async () => {
    expect(await urlAfter(`/browse${NCA_APP_ROOT}`, (v) => updateNcaAnalysisPath(v, 'X')),
      `/browse${NCA_APP_ROOT}/analyses/X`);
  });
});

// The source-of-truth dedup underpinning label-keyed routing + tree captions. One
// pure function backs every gate — study create/rename, add-to-study, analysis rename.
category('NCA Studies: label uniqueness', () => {
  test('uniqueLabel: free name unchanged; collision → label-<n> at the first free slot', async () => {
    expect(uniqueLabel('Study A', new Set<string>()), 'Study A'); // free → unchanged
    expect(uniqueLabel('Study A', new Set(['Study A'])), 'Study A-2'); // taken → -2
    expect(uniqueLabel('Study A', new Set(['Study A', 'Study A-2'])), 'Study A-3'); // next taken → -3
    expect(uniqueLabel('X', new Set(['X', 'X-2', 'X-4'])), 'X-3'); // fills the -3 gap
  });
});

category('NCA Studies: Browse-tree node meta', () => {
  test('only our study / analysis values preview; the app node and markers do not', async () => {
    expect(asNcaNodeMeta(null) === null, true, 'a ⚠ marker (null value)');
    // The app node carries the platform's own value — it once opened `undefined` (QA 2026-10-21).
    const app = DG.Func.find({package: 'Nca', name: 'ncaApp'})[0];
    expect(asNcaNodeMeta(app) === null, true, 'the app node (a Func value)');
    expect(asNcaNodeMeta({kind: 'analysis'}) === null, true, 'an analysis value without an id');
    expect(asNcaNodeMeta({kind: 'standalone', analysisId: 'a'}) === null, true, 'an unknown kind');
    expect(asNcaNodeMeta({kind: 'study', studyId: 's'})?.kind, 'study');
    expect(asNcaNodeMeta({kind: 'analysis', analysisId: 'a'})?.kind, 'analysis');
  });
});
