/**
 * NCA Hub app view — the persistent entry point registered as `ncaApp()` in
 * `package.ts`. Mirrors the DiffStudio Hub pattern: header + "New analysis"
 * form + "Projects" section (the sample cards grid, currently hidden). Replaces
 * the UC-01 Phase 6 launcher card.
 *
 * @see packages/DiffStudio/src/hub.ts (class
 *      skeleton + render lifecycle + buildModelCard + cards grid container)
 *
 * **Header.** Built with `u2.appHeader` from `@datagrok-libraries/utils` —
 * same idiom DiffStudio uses (`hub.ts:62-74`). Spec §5.1 calls for "3-4 lines
 * about NCA": parameter scope, target audience, single-compound serial-
 * sampling scope, PKNCA validation. Description rendered as markdown bullet
 * list (one bullet per point) so the header reads as a scannable feature list,
 * not prose.
 *
 * **Library cards.** Two static entries (theoph + indometh) with PNG icons at
 * `files/icons/<name>.png` (user-added 2026-05-18). Double-click
 * triggers `loadSampleByName(file) → launchAnalysis(df, viewName)`.
 *
 * **CSS.** Classes are in `css/app-styles.css`; loaded via the
 * `//@ts-ignore` import shim at the top of this file (webpack picks up
 * through `style-loader` + `css-loader` — see webpack.config.js).
 */

// @ts-ignore  — CSS module not declared in tsconfig; webpack resolves at build time.
import '../../css/app-styles.css';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {u2} from '@datagrok-libraries/utils/src/u2';

import {_package} from '../package';
import {loadSampleByName, launchAnalysis, ensureDemoTables} from './nca-launch';
import {openHubHelp, restoreHubHelpSync, buildDemoBalloonElement} from './nca-hub-help';
import {
  renderStudiesListInto, createStudyFlow, blueIcon, UNGROUPED_SECTION, UNGROUPED_TOOLTIP,
} from './nca-studies-list';
import {renderAnalysesListInto} from './nca-analyses-list';
import {NCA_STUDIES_CHANGED_EVENT} from '../adapter';

interface SampleEntry {
  readonly file: string;
  readonly viewName: string;
  readonly icon: string;
  readonly description: string;
}

/** Built-in Library samples — file lives under `files/`, icon under
 *  `files/icons/`. Same set the UC-01 Phase 6 launcher offered. */
const SAMPLE_FILES: ReadonlyArray<SampleEntry> = [
  {
    file: 'theoph.csv',
    viewName: 'Theophylline',
    icon: 'theoph.png',
    description: 'Oral PK study with 12 subjects following a single PO dose.',
  },
  {
    file: 'indometh.csv',
    viewName: 'Indomethacin',
    icon: 'indometh.png',
    description: 'IV bolus PK study with 6 subjects following a single intravenous dose.',
  },
];

/** NCA Hub app view — see header doc-comment. */
export class NcaAppView {
  /** `view.temp` key marking the NCA Hub view so `ncaApp()` can find it in
   *  `grok.shell.views` without relying on the view name (which may collide
   *  with views from other apps) or JS-reference equality (wrappers may be
   *  re-instantiated by the platform). `temp` is a dart-side map, so the
   *  marker survives JS-wrapper re-creation. Mirrors DiffStudio.HUB_TAG. */
  static readonly HUB_TAG = 'nca-hub';

  /** Find the open NCA Hub view, if any. */
  static findHubView(): DG.View | null {
    for (const v of grok.shell.views)
      if (v.temp?.[NcaAppView.HUB_TAG] === true) return v;
    return null;
  }

  /** Mark a freshly created view as the NCA Hub. */
  static stampHubView(view: DG.View): void {
    view.temp[NcaAppView.HUB_TAG] = true;
  }

  readonly view: DG.View;
  private readonly root: HTMLDivElement;
  /** Subscription that re-renders the Projects list when the Hub re-activates. */
  private viewActivatedSub?: {unsubscribe(): void};
  /** Subscription that re-renders the studies + analyses lists whenever a study or
   *  analysis is created / renamed / deleted / added-to / moved, from any surface
   *  (the coordinator fires `NCA_STUDIES_CHANGED_EVENT` on save + delete). Fires
   *  UNCONDITIONALLY — the Hub stays current with no manual refresh (2026-10-05). */
  private studiesChangedSub?: {unsubscribe(): void};

  constructor() {
    this.root = ui.divV([]);
    this.view = DG.View.fromRoot(this.root);
    this.view.name = 'NCA';
    NcaAppView.stampHubView(this.view);
  }

  /** Render header + New analysis + Library. Idempotent — safe to re-call. */
  async render(): Promise<void> {
    ui.empty(this.root);
    this.root.appendChild(this.buildHeader());
    this.root.appendChild(this.buildNewAnalysis());
    this.root.appendChild(this.buildLibrary());
  }

  private buildHeader(): HTMLElement {
    const header = u2.appHeader({
      iconPath: `${_package.webRoot}package.png`,
      description:
        '- Perform non-compartmental analysis of preclinical PK data\n' +
        '- Compute PK parameters from concentration–time profiles\n' +
        '- Run a quick analysis on any table\n' +
        '- Group multiple analyses into a shareable Study\n',
    });
    header.classList.add('nca-hub-header');
    // "Learn more", directly under the bullets. u2.appHeader HAS a Learn-more
    // slot (`IAppInfo.learnMoreUrl`) but it only accepts a URL and renders an
    // external link that navigates away — we need the in-app help panel. So
    // build the very element u2 builds for its own Learn more (`ui.link` + the
    // same 15px indent, u2.js appHeader) and append it to u2's description
    // panel — the container that holds the markdown bullets — so the link lands
    // exactly where u2's own would. `ui.link` takes a Function target
    // (ui.js:540), which is the seam `learnMoreUrl: string` closes off.
    // Fallback to the header root if u2 ever renames its description class.
    const learnMore = ui.link('Learn more', () => openHubHelp(),
      'What an analysis is, and how to get from a table to a Study — opens the help panel',
      {style: {marginLeft: '15px'}});
    (header.querySelector('.u2-app-header-description') ?? header).appendChild(learnMore);
    return header;
  }

  private buildNewAnalysis(): HTMLElement {
    // ui.input.table — initial value omitted (DG types accept DataFrame |
    // undefined, not null); the picker renders empty / "Select table..." by
    // default. Run button gating below handles the null/undefined case.
    const tableInput = ui.input.table('Table', {
      tooltipText: 'Set the table containing subject IDs, time points, and concentration values.',
    });

    // No "Analysis name" input (user-directed 2026-10-06, hub-ungrouped-polish): the view is titled
    // with the table name, and the analysis is named where it matters — the first SAVE pre-fills the
    // table name, Add to study proposes a settings-spec name. Run is the Hub's ONE primary action, so
    // it is DG's raised `bigButton` (user-directed 2026-10-06).
    const runBtn = ui.bigButton('Run', () => {
      const df = tableInput.value;
      if (df === null)
        return; // gated by disabled state; defensive guard.
      // Provenance for "Save to project": an in-session table id. The share dialog
      // correctly flags this as owner-only (recipients can't resolve it).
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const tableId = (df as any).id as string | undefined;
      launchAnalysis(df, df.name,
        {kind: 'tableId', ref: (tableId != null && tableId.length > 0) ? tableId : df.name});
    });
    // Single STANDARD DG tooltip (lazy → reflects the disabled state). Do NOT also
    // set `runBtn.title` — the raw HTML `title` attribute renders a SECOND, native
    // tooltip on top of the DG one (the duplicate the user saw).
    ui.tooltip.bind(runBtn, () => tableInput.value === null ?
      'Pick a table first' :
      'Run NCA on the selected table');

    const refreshRunEnabled = (): void => {
      runBtn.disabled = tableInput.value === null;
    };
    tableInput.onChanged.subscribe(refreshRunEnabled);
    refreshRunEnabled();

    // Demo tables: an "Open demo dataset" link under the Table input, on the Run row (user-directed
    // 2026-10-12, replacing the ★ icon in the input's icon strip — a named link reads as the way in
    // to demo data, an unlabelled star did not). Same action and tooltip as the icon it replaces.
    const demoLink = ui.link('Open demo dataset', () => {
      void (async () => {
        try {
          const tables = await ensureDemoTables();
          // Put the user one click from Run — the whole point of the icon. Seed
          // the picker ONLY while it is empty: a table the user already chose is
          // theirs, and silently swapping it would be a worse surprise than an
          // extra click. Re-gate Run explicitly rather than trusting the
          // programmatic set to raise onChanged.
          if (tableInput.value === null && tables.length > 0) {
            tableInput.value = tables[0];
            refreshRunEnabled();
          }
          // Markdown, not a joined name list: the balloon is the demo's first
          // touch, so it says what each table is FOR — and it is the surface
          // that carries the simulated-data qualifier for the multiroute set.
          grok.shell.info(buildDemoBalloonElement());
        } catch (err) {
          grok.shell.error(
            `Could not open demo tables: ${err instanceof Error ? err.message : String(err)}`);
        }
      })();
    }, 'Open demo tables');

    // Header tooltip explains the entity hierarchy: an analysis is the unit.
    const analysisTitle = ui.h1('New analysis');
    ui.tooltip.bind(analysisTitle, () =>
      'An analysis is the basic unit — one table mapped to columns and computed. ' +
      'Group several analyses into a Study.');

    // Title, Table input and Run all live INSIDE the ui.form (user-directed 2026-10-06), so the
    // form's width is the one edge: the Table input stretches to it (app-styles.css) and Run is
    // pushed to it by DG's own `.ui-form .ui-btn {margin-left: auto}`. `ui.form` takes inputs only,
    // so the title and the Run row are added to its root. The Run row is a plain `ui.div`, NOT
    // `ui.divH` — DG pads a `.ui-form > .ui-div.d4-flex-row` 8px on the right, which would pull Run
    // off the input's right edge. The demo link opens the row, Run closes it.
    const form = ui.form([tableInput]);
    form.classList.add('nca-hub-section', 'nca-hub-new-analysis');
    form.prepend(analysisTitle);
    const runRow = ui.div([demoLink, runBtn], 'nca-hub-run-row');
    form.append(runRow);
    // The link starts under the select, not under the label. The label column's width is written
    // by DG's form-resize script, so it is measured, not restated in CSS.
    ui.onSizeChanged(form).subscribe(() => {
      const dx = tableInput.input.getBoundingClientRect().left - runRow.getBoundingClientRect().left;
      demoLink.style.marginLeft = `${Math.max(0, dx)}px`;
    });
    return form;
  }

  private buildLibrary(): HTMLElement {
    // Section "Studies" (NCA Studies workspace, AC-C1). The live studies table
    // renders into `studiesHost` (async); the sample-data cards
    // (Theophylline / Indomethacin) stay built-but-hidden so their launch
    // handlers remain live for an easy restore.
    const cards = SAMPLE_FILES.map((s) => this.buildSampleCard(s));
    const grid = ui.div(cards, 'nca-hub-grid');
    grid.style.display = 'none';
    const studiesHost = ui.div([], 'nca-projects-host');
    // nca-workflow-stages F-S6.1 (D4): standalone analyses, under Studies, same table style. Titled
    // "Ungrouped" — the analyses NOT in a study; "Analyses" read as "all analyses" (user-directed
    // 2026-10-06, hub-ungrouped-polish). The name is shared with the Browse-tree group.
    const analysesHost = ui.div([], 'nca-projects-host');
    const analysesTitle = ui.h1(UNGROUPED_SECTION);
    analysesTitle.style.margin = '24px 0 0 0';
    ui.tooltip.bind(analysesTitle, () => UNGROUPED_TOOLTIP);
    analysesTitle.style.display = 'none'; // shown once the list has entries
    // Re-read the lists and re-render, so the tables reflect studies / analyses created or changed
    // from OTHER views (e.g. "Add to study" from an analysis pane). The two lists share one read
    // (the repository coalesces concurrent list reads).
    const refresh = async (): Promise<void> => {
      try {
        // No empty Ungrouped block (user-directed 2026-09-29): title + list appear with the first
        // entry and go again when the last one leaves (the list re-renders itself — 2026-09-30).
        await Promise.all([
          renderStudiesListInto(studiesHost),
          renderAnalysesListInto(analysesHost, (n) => {
            analysesTitle.style.display = n > 0 ? '' : 'none';
            analysesHost.style.display = n > 0 ? '' : 'none';
          }),
        ]);
      } catch (err) {
        const msg = `Could not load studies: ${err instanceof Error ? err.message : String(err)}`;
        // The list now refreshes on every study/analysis change, even when the Hub is not the
        // visible view. A balloon for a transient failure of a HIDDEN Hub would appear with no
        // visible cause, so only surface it when the Hub is in front; else log quietly.
        if (grok.shell.v === this.view) grok.shell.error(msg);
        else console.warn(`[NCA Hub] ${msg}`);
      }
    };
    void refresh();
    // Re-render whenever the Hub becomes the current view again (cheap — a few small
    // file reads), so a study created/changed in another view shows without a manual
    // reload. Stored so a future teardown could unsubscribe; the Hub is session-long.
    this.viewActivatedSub?.unsubscribe();
    this.viewActivatedSub = grok.events.onCurrentViewChanged.subscribe(() => {
      if (grok.shell.v === this.view)
        void refresh();
      else
        // Left the Hub → restore the help panel's syncCurrentObject (openHubHelp
        // disabled it globally; the disable is scoped to the Hub's lifetime, not
        // the session — decision-audit fix, model-catalog.ts:69-76 pattern).
        restoreHubHelpSync();
    });
    // Re-render on EVERY study/analysis change (create / rename / delete / add-to-study / move),
    // from any surface — the coordinator fires NCA_STUDIES_CHANGED_EVENT on save + delete. The Hub
    // refreshes UNCONDITIONALLY (no `grok.shell.v === this.view` guard), exactly like the Browse-tree
    // (nca-tree-browser.ts), so adding an entity updates the Hub automatically — no manual refresh
    // (user-directed 2026-10-05, `study-add-multiple-analyses`). Re-rendering a non-visible Hub is a
    // few small file reads into its own host; the user sees it already current on return.
    this.studiesChangedSub?.unsubscribe();
    this.studiesChangedSub = grok.events.onCustomEvent(NCA_STUDIES_CHANGED_EVENT).subscribe(() => {
      void refresh();
    });

    // Header: "Studies" + New (＋), blue icon inline (AC-C1, user-directed 2026-06-29). The manual
    // Refresh icon was removed 2026-10-05 — the list is now reactive (see the event subscription
    // above). ＋ creates a new study (same as the former "New study" button).
    const newIcon = blueIcon('plus', () => {void createStudyFlow(studiesHost);}, 'New study');
    // Zero the h1's default block margins so `align-items:center` lines the icons
    // up with the title text (not its margin box) — per user-directed 2026-06-29.
    const title = ui.h1('Studies');
    title.style.margin = '0';
    // Header tooltip mirrors the New-analysis one — a Study is the collection.
    ui.tooltip.bind(title, () =>
      'A Study is a named, shareable collection of analyses you can re-open later.');
    // The action-icon group (now just ＋). Kept as a divH so the title→icons gap and
    // future icons stay consistent with the former two-icon layout.
    const icons = ui.divH([newIcon]);
    icons.style.alignItems = 'center';
    icons.style.gap = '12px';
    const header = ui.divH([title, icons], 'nca-studies-header');
    header.style.alignItems = 'center';
    header.style.gap = '8px';

    return ui.divV([
      header,
      studiesHost,
      analysesTitle,
      analysesHost,
      grid,
    ], 'nca-hub-section');
  }

  private buildSampleCard(s: SampleEntry): HTMLElement {
    const iconUrl = `${_package.webRoot}files/icons/${s.icon}`;
    const icon = ui.iconImage('nca', iconUrl);
    icon.classList.add('nca-hub-card-icon');
    const label = ui.label(s.viewName);
    const header = ui.divH([icon, label], 'nca-hub-card-header');
    const card = ui.divV([header], 'nca-hub-card');

    const onLaunch = (): void => {
      void (async () => {
        try {
          const df = await loadSampleByName(s.file);
          // Bundled sample lives in shared App storage → a shareable platformPath.
          launchAnalysis(df, s.viewName,
            {kind: 'platformPath', ref: `System:AppData/Nca/${s.file}`});
        } catch (err) {
          grok.shell.error(
            `Could not load sample '${s.file}': ` +
            (err instanceof Error ? err.message : String(err)),
          );
        }
      })();
    };

    card.addEventListener('dblclick', onLaunch);

    ui.tooltip.bind(card, () => {
      const tipIcon = ui.iconImage('nca', iconUrl);
      tipIcon.classList.add('nca-hub-tooltip-icon');
      const tipName = ui.label(s.viewName);
      tipName.classList.add('nca-hub-tooltip-name');
      const tipHeader = ui.divH([tipIcon, tipName], 'nca-hub-tooltip-header');
      const tipDesc = ui.divText(s.description, 'nca-hub-tooltip-desc');
      return ui.divV([tipHeader, tipDesc]);
    });

    return card;
  }
}
