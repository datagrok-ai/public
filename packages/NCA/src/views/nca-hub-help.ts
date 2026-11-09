/**
 * Hub help panel — the content behind the Hub header's "Learn more" link
 * (topic: nca-hub-onboarding). Answers the cold-start question a first-time
 * user actually has: *what is the unit of work here, and what do I press?*
 *
 * **Why the markdown is an inline constant and not a shipped `.md` file.**
 * `showHelp` takes either a URL or an HTMLElement (shell.d.ts:190). The
 * element path has a working precedent (ApiSamples
 * `shell/windows/shell-windows.js`, and `ui.markdown` is already the renderer
 * used by `tour/ui-describer.ts:158`); the *package-file URL* path has none —
 * the Bio precedent (`bio01b:158`) points at the platform's own `/help/**`
 * wiki, not at package files. Shipping the prose to `System:AppData/Nca/` and
 * fetching it would add an async round-trip and a silent failure mode (a
 * publish that misses the file leaves Learn more broken) for text that never
 * varies at runtime. A template literal is as editable as a `.md` at this size,
 * and it cannot drift out of the bundle.
 *
 * **Glyphs.** `★` (U+2605) and `▶` (U+25B6) are text-presentation symbols
 * verified to render in the help panel's font. Two things NOT to do, both
 * learned the hard way:
 *  - `＋` (U+FF0B FULLWIDTH PLUS SIGN) does NOT render — it shipped in the
 *    add-to-study step and showed as a blank gap (user-reported 2026-07-17), a
 *    silent failure: no error, no tofu box, just a missing word. Fullwidth
 *    forms (U+FF00–U+FFEF) are banned here and a test enforces it.
 *  - An emoji (e.g. `ℹ️` U+2139+U+FE0F) for the info/guide icon: emoji
 *    presentation is font-dependent, so unlike ★/▶ it is not guaranteed.
 * Three controls the user asked to show as their REAL icons (2026-07-17) — the
 * info/guide `info-circle`, the add-to-study `layer-plus`, and the share
 * `share` — are the actual `ui.iconFA(...)` nodes the app uses (nca-mapping-
 * pane.ts:2828 / :2860, nca-studies-list.ts:201), NOT look-alike glyphs. An
 * icon cannot live inside a markdown string, so {@link buildHubHelpElement}
 * renders the markdown and splices each icon node in at its token
 * (`%INFO_ICON%`, `%ADD_ICON%`, `%SHARE_ICON%`). That removes the font gamble
 * entirely: the panel shows the actual app icons, or a splice fails LOUDLY (a
 * visible `%…%` token a test asserts against), never a blank gap. For every
 * other control, prefer naming it by its tooltip (font-independent, and what
 * the user sees on hover).
 *
 * **Copy discipline.** Every control named below maps to a REAL label or
 * tooltip — the compute control is the `Confirm` button at the bottom of the mapping pane
 * (tooltip "Confirm the setup and compute …"; the ribbon ▶ icon was removed, nca-workflow-stages D17).
 * The spliced icons carry their tooltips on hover via `ui.iconFA`, so hovering
 * the help icon shows the same text as hovering the real control (the share one
 * is generic — see {@link HELP_ICONS}).
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {DEMO_TABLE_FILES} from './nca-demo-tables';

/**
 * Help-panel body (copy user-directed 2026-07-17). Three parts:
 *  - a one-line intro naming what the app is,
 *  - Quick start (the cold-start algorithm — table → compute → save into a
 *    study → compare on the Studies view; the analysis→Study relationship is
 *    carried by steps 5–6 rather than a separate intro paragraph),
 *  - Computations (method-transparency write-up, user-authored 2026-07-18).
 *    Its claims are COUPLED to the code they describe: thresholds and behavior
 *    were audited against `adapter/defaults.ts` + sci-comp `nca/lambda-z.ts` /
 *    `compute-nca.ts` (8/10 claims code-confirmed — e.g. minRSquared 0.85
 *    reject-not-flag, adjRSquaredFactor 1e-4 tie, minPoints 3, extrap 20/50%,
 *    Vss IV-only, CL/F routing, MRT−Tinf/2). Two were corrected to match code:
 *    the BLQ line was softened from a "same AUC as ICH M13A" equivalence (no
 *    fixture proves it) to a convention statement, and "four warnings" was
 *    scoped to "computation-quality" (the provenance.warnings stream also
 *    carries adapter advisories). If those defaults change, this copy must too.
 */
export const HUB_HELP_MD: string = `
# NCA

NCA is an application for non-compartmental analysis of preclinical pharmacokinetic data.

## Quick start

1. **Pick a table** under *New analysis*. To use demo data, click *Open demo
   dataset* under the picker.
2. **Click Run** to open the analysis view. The strip at the top shows the three
   stages — *Setup › Review λz › Report* — and what to do next. On *Setup*, map
   *Subject*, *Time*, and *Concentration* in the mapping pane; add *Dose* and
   *Route* if needed.
3. **Click the %INFO_ICON% (Run interactive guide) icon** on the top panel to
   walk through the mapping pane and the results in place.
4. **Click Confirm** at the bottom of the mapping pane
   (*Confirm the setup and compute*). The results grid and the
   concentration–time profiles appear, and the analysis moves to *Review λz*.
5. **Review the λz fits that need a decision.** A Flagged or Not estimated
   profile is *To review*: edit its λz points, or — with a reason — accept it as
   is, exclude λz (its t½, AUC∞, CL, Vz, MRT and Vss then stay out of every
   statistic), or acknowledge it, in the context panel. *Next to review* walks
   the remaining profiles.
6. **Save at any stage** with the SAVE button, leftmost on the top panel. The
   first save goes to *Ungrouped* on the app's main view; reopening continues
   where you stopped.
7. **Add the analysis to a study** on the *Report* stage: click the %ADD_ICON%
   icon and pick an existing study, or create a new one. The analysis moves from
   *Ungrouped* into the study.
8. **Find your studies and analyses** on the app's main view. Click a study to
   open it and compare analyses side by side. Click a row in the study grid to
   see that compound's NCA parameters and concentration–time profile on the
   context panel.
9. **Share studies or individual analyses with colleagues.** Click the
   %SHARE_ICON% icon to copy a link.

## Computations

Every parameter is validated against PKNCA (AUC and Cmax within 0.1%, λz within
0.5%, derived within 1%). Here is how each number is produced — and where it may
differ from another tool.

**Clearance & volumes**
- \`CL = Dose / AUC∞\`, from the *observed* last concentration. Extravascular
  routes are reported as \`CL/F\`, never \`CL\` — the value is confounded by
  unknown bioavailability, so calling it \`CL\` would be wrong.
- \`Vz = Dose / (λz · AUC∞)\`. \`Vss = Dose · MRT / AUC∞\`, **IV only** — blank for
  extravascular, where an "oral Vss" would really be Vss/F.

**Half-life (λz)** — *PKNCA convention*
- \`t½ = ln(2) / λz\`. λz is the best-fit terminal window by adjusted R²: start
  from the last 3 points, extend backwards, keep the highest adjusted R²
  (ties within 0.0001 → the longer window wins).
- The peak point is excluded from the fit. For IV bolus that peak is the
  back-extrapolated C₀, so the first *measured* point is kept.
- **A window below adjusted R² 0.85 is rejected, not flagged.** If none
  qualifies, t½ and everything downstream — AUC∞, CL, Vz, MRT, Vss — is left
  blank. When a value is blank, that is usually why. (Stricter than PKNCA, which
  reports the number and flags it at a separate step — by design.)

**MRT** — *industry-standard convention*
- \`MRT = AUMC∞ / AUC∞\`, minus \`Tinf/2\` for infusions. For extravascular it
  includes mean absorption time, so it is never fed into Vss.

**Computation-quality warnings** — four are raised:
- \`%AUC extrapolated\` above 20% (error above 50%)
- \`%AUMC extrapolated\` above 20% (MRT / Vss are fragile)
- λz used the minimum allowed number of points
- more than 50% of samples below the limit of quantification

Some checks a company SOP might expect — span ratio ≥ 2, sampling ≥ 3 × t½,
pre-dose > 5% of Cmax — are **not computed**; this is stated rather than implied.

**Below the limit of quantification (BLQ)**
- Handled by position (before first measurable / embedded / trailing),
  consistent with the zero-fill convention of ICH M13A.

## Demo studies

Five demo studies of one simulated discovery program (Series X, lead
DGK-0007) can be installed on the server: a rat IV/PO screen of eight analogs, a
mouse sparse study, a rat dose escalation, a 28-day TK study and an IV infusion.
They are shared with all users (View) and show finished analyses with recorded
λz decisions; the TK study leaves three profiles *To review*. All data are
simulated; the structures are real molecules with fictitious PK.

The demo studies are not created automatically. Run \`ncaSeedDemoStudies()\`
from the console once per server — any user can; the first caller owns the demo.
After every update of the NCA calculation library (sci-comp), the owner or an
administrator runs it again: out-of-date demo studies are re-created with new ids
(copied links stop working). Otherwise every demo row reads *stale (refreshes on
open)*, which View-only users cannot refresh.
`;

/**
 * The ★ balloon body: a lead-in plus one bullet per demo table saying what that
 * table UNLOCKS. Names and purpose clauses are projected from
 * `DEMO_TABLE_FILES`, so a renamed table cannot drift out of the balloon.
 *
 * **This copy is the ONLY user-visible disclosure that the multiroute data is
 * simulated** — AC-F4c measured that the frame-level `PROVENANCE_TAG` renders in
 * no platform surface at all, so nothing else in the app says it. That is why
 * the caveat is projected from `DemoTableEntry.qualifier` rather than typed
 * here, and why a test guards its presence, its italics, and its 4-word budget
 * in both directions.
 *
 * The wording went through the user twice (2026-07-28): a full sentence plus a
 * meta-clause about its own persistence was cut entirely, then restored as one
 * short italic line capped at four words. Do not lengthen it and do not drop it
 * — either way is a user decision, not a cleanup.
 *
 * **Still not covered, and not to be claimed as covered:** the balloon is
 * dismissible and does not follow the data into a re-opened analysis, an export,
 * or a screenshot. See `files/REGEN.md`.
 *
 * Safe as an eager `const` (mirroring {@link HUB_HELP_MD}) only because
 * `DEMO_TABLE_FILES` lives in a dependency-free leaf module — see
 * `nca-demo-tables.ts` for why that matters and what it cost to learn.
 */
export const DEMO_BALLOON_MD: string = (() => {
  const bullets = DEMO_TABLE_FILES
    .map((e) => {
      // Italic, and deliberately terse — see DemoTableEntry.qualifier for the
      // 4-word budget and why it exists. The balloon is the ONLY place this
      // caveat is visible to a user: the frame provenance tag was measured to
      // render in no platform surface at all.
      const caveat = e.qualifier === undefined ? '' : ` *${e.qualifier}*`;
      return `- **${e.name}** — ${e.unlocks}${caveat}`;
    })
    .join('\n');
  return `${DEMO_TABLE_FILES.length} demo tables are open.\n\n${bullets}\n`;
})();

/** Render {@link DEMO_BALLOON_MD} for `grok.shell.info`, which accepts an
 *  element (`shell.d.ts:82` types it `any`; the neighbouring `warning`/`error`
 *  are `string | HTMLElement`, and element-valued balloons have live precedents
 *  in Compute and Plates). Markdown, so the bullets read as a list rather than
 *  one run-on line. */
export function buildDemoBalloonElement(): HTMLElement {
  return ui.markdown(DEMO_BALLOON_MD);
}

/** Tokens in {@link HUB_HELP_MD} marking where a real FontAwesome icon node is
 *  spliced in. Plain ASCII on purpose: markdown leaves them untouched, and if a
 *  splice ever fails the token renders as visible text (loud), not a blank gap. */
export const INFO_ICON_TOKEN = '%INFO_ICON%';
export const ADD_ICON_TOKEN = '%ADD_ICON%';
export const SHARE_ICON_TOKEN = '%SHARE_ICON%';

/** token → the exact `ui.iconFA` the matching app control uses, so the help
 *  shows the ACTUAL app icon (and its real tooltip on hover), not a look-alike
 *  glyph. `info-circle` = the guide button (nca-mapping-pane.ts:2828);
 *  `layer-plus` = the add-to-study button (nca-mapping-pane.ts:2860); `share` =
 *  the copy-link icon on both study rows (nca-studies-list.ts:201) and analysis
 *  rows (nca-study-overview.ts:1205). The share tooltip is generic on purpose:
 *  it stands for two real controls whose tooltips differ only in the trailing
 *  noun ("...to this study" / "...to this analysis"). */
const HELP_ICONS: ReadonlyArray<{token: string; fa: string; tooltip: string}> = [
  {token: INFO_ICON_TOKEN, fa: 'info-circle', tooltip: 'Run interactive guide'},
  {token: ADD_ICON_TOKEN, fa: 'layer-plus', tooltip: 'Add this analysis to a study'},
  {token: SHARE_ICON_TOKEN, fa: 'share', tooltip: 'Copy a shareable link'},
];

/** Replace every `token` inside `root`'s text nodes with a fresh `ui.iconFA`
 *  node — the one splice mechanism both help icons go through. A markdown string
 *  cannot carry a live element, so the icon is injected after render. */
function spliceIcon(root: HTMLElement, token: string, fa: string, tooltip: string): void {
  const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
  const hits: Text[] = [];
  while (walker.nextNode()) {
    const t = walker.currentNode as Text;
    if (t.nodeValue !== null && t.nodeValue.includes(token))
      hits.push(t);
  }
  for (const node of hits) {
    const parts = (node.nodeValue ?? '').split(token);
    const frag = document.createDocumentFragment();
    parts.forEach((part, i) => {
      if (part.length > 0)
        frag.appendChild(document.createTextNode(part));
      if (i < parts.length - 1)
        frag.appendChild(ui.iconFA(fa, null, tooltip));
    });
    node.parentNode?.replaceChild(frag, node);
  }
}

/**
 * Render {@link HUB_HELP_MD} and splice each real FontAwesome icon into it at
 * its token. Exported so a test can assert every icon lands and no token leaks.
 */
export function buildHubHelpElement(): HTMLElement {
  const el = ui.markdown(HUB_HELP_MD);
  for (const {token, fa, tooltip} of HELP_ICONS)
    spliceIcon(el, token, fa, tooltip);
  return el;
}

/** The platform help panel's `syncCurrentObject` value from BEFORE the Hub first
 *  disabled it — captured once so {@link restoreHubHelpSync} can put it back.
 *  `null` = the Hub has not disabled it (nothing to restore). */
let priorSyncCurrentObject: boolean | null = null;

/**
 * Show the Hub help in the platform help panel.
 *
 * `windows.help` (ShellHelpPanel, shell.d.ts:182) is NOT the same thing as
 * `windows.showHelp` (a plain visibility boolean, shell.d.ts:234) — the panel
 * object is the one that can be handed content, so the visibility toggle goes
 * through `panel.visible`.
 *
 * **`syncCurrentObject = false` is load-bearing, not cleanup.** Without it the
 * FIRST click on a CLOSED panel showed the wrong help (user-reported 2026-07-17,
 * two screenshots): the panel defaults to `syncCurrentObject = true`, so opening
 * it makes the platform populate it from `grok.shell.o` — the default package
 * help ("Datagrok: Swiss Army Knife for Data") — and that open-time sync
 * overwrites whatever `showHelp()` set. The second click found the panel already
 * open, no open-time sync fired, and `showHelp()` won — exactly the reported
 * two-click behavior. Turning sync off before `showHelp` is the platform-
 * standard order (Tutorials/src/package.ts:548, GIS/gis-demo.ts, Charts/demos/
 * demo.ts). Keeping it off while the Hub is the active view also stops the NCA
 * help being swapped out when the user clicks a study row (the current object
 * changes) — a second latent bug closed here.
 *
 * **Scope is the Hub's lifetime, NOT the session.** `syncCurrentObject` is a
 * global platform singleton; leaving it off forever would break contextual (F1)
 * help sync for every OTHER package after one "Learn more" click. So the prior
 * value is captured here and {@link restoreHubHelpSync} puts it back when the
 * user navigates away from the Hub — the restore-on-leave pattern a normal
 * embedded feature uses (compute-utils model-catalog.ts:69-76 also restores on
 * leave, though it assumes the default rather than remembering the real prior
 * value as this does), not the never-restore of the full-shell-takeover demos
 * above. Restoring on navigate-away (not synchronously after `showHelp`) is what
 * avoids re-arming the open-time clobber.
 */
export function openHubHelp(): void {
  const panel = grok.shell.windows.help;
  panel.visible = true;
  if (priorSyncCurrentObject === null)
    priorSyncCurrentObject = panel.syncCurrentObject;
  panel.syncCurrentObject = false;
  panel.showHelp(buildHubHelpElement());
}

/**
 * Restore the platform help panel's `syncCurrentObject` to whatever it was
 * before the Hub last disabled it, and re-arm the capture. No-op if the Hub
 * never disabled it. The Hub view calls this when it stops being the current
 * view (its existing `onCurrentViewChanged` subscription), so the disable in
 * {@link openHubHelp} is scoped to the Hub and never leaks session-wide.
 */
export function restoreHubHelpSync(): void {
  if (priorSyncCurrentObject === null)
    return;
  grok.shell.windows.help.syncCurrentObject = priorSyncCurrentObject;
  priorSyncCurrentObject = null;
}
