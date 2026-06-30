---
phase: 15-read-only-publishing-foundation
plan: 05
type: execute
wave: 4
depends_on: ["15-01", "15-04"]
files_modified:
  - src/publishing/share-dialog.ts
autonomous: true
requirements: [PUB-04, PUB-05, PUB-08, PUB-10, PUB-13]
must_haves:
  truths:
    - "`showShareForReviewDialog(df)` is `async (df) => Promise<void>` (W-6 fix — it `await`s grok.dapi.groups.list() inside) and opens a `ui.dialog` with target text input, reviewer-group ChoiceInput (populated from `grok.dapi.groups.list()`), and note textarea"
    - "Dialog detects republish on open via `findPriorShare(target, group)`; prefills target/group/note + renders banner `'⚠ This will publish as v<N+1> and supersede <prior name>'`"
    - "Dialog renders a reactive confirmation summary block (W-4 fix for PUB-08) showing computed slug + project name + reviewer group + supersede status; updates as the user edits any input"
    - "On OK: calls `publishAnalysis(df, opts)`; surfaces success with mailto link (PUB-13, P2 sub-task) and project link; surfaces errors verbatim including D-03 exact string"
    - "Every reviewer-touchable string passes the Pitfall 14 jargon audit (banned: `DataFrame`, `tag`, `semType`, `ACL`, `viewer factory`)"
    - "Reviewer-group list is filtered to groups the publishing user can administer (D-02)"
    - "`buildMailtoUrl` is imported from `./publish-state` (NOT redefined here per B-1 fix — Plan 01 owns the export so wave-2 callers can import it without cycles)"
  artifacts:
    - path: "src/publishing/share-dialog.ts"
      provides: "showShareForReviewDialog dialog flow (async)"
      exports: ["showShareForReviewDialog"]
  key_links:
    - from: "src/publishing/share-dialog.ts"
      to: "src/publishing/publish-project.ts"
      via: "publishAnalysis call on OK"
      pattern: "publishAnalysis\\("
    - from: "src/publishing/share-dialog.ts"
      to: "src/publishing/publish-state.ts"
      via: "buildMailtoUrl + DE_COMPLETE_TAG + findPriorShare + slugifyTarget + PublishOptions imports"
      pattern: "import.*publish-state"
    - from: "src/publishing/share-dialog.ts"
      to: "grok.dapi.groups.list"
      via: "reviewer group ChoiceInput population"
      pattern: "groups\\.list"
---

<objective>
Build the user-facing dialog: target input + reviewer-group ChoiceInput + note + republish-detection banner + reactive confirmation summary (W-4). On OK, calls `publishAnalysis(df, opts)` and surfaces the result. Mailto link (PUB-13, P2 sub-feature) appears in the success state.

**Revision notes:**
- **B-1:** `buildMailtoUrl` no longer defined here — imported from `./publish-state` (Plan 01 owns it so Plan 06 can use it without wave-2-imports-wave-4 cycle).
- **B-3:** This plan is no longer marked `priority: P2` at plan level. The plan covers PUB-04, PUB-05, PUB-08, PUB-10 (all P1) + PUB-13 (P2). Only the mailto-related Task 2 sub-step is P2-tagged; the rest of the work is P1.
- **W-4:** Added reactive confirmation summary div (PUB-08) that displays computed slug + project name + reviewer group + supersede status. Refreshes on every input change. The republish banner alone is NOT a confirmation summary — it only fires on republish detection.
- **W-6:** `showShareForReviewDialog` is now `async`. Plan 07 menu handler must `await` it (or `void` it for fire-and-forget).
- **I-10:** Precondition check uses the `DE_COMPLETE_TAG` constant from publish-state.ts, not the inline `'proteomics.de_complete'` string.

Purpose: This is the analyst's only entry point to the publish flow. Every reviewer-touchable string flows through here (target label, group label, note placeholder, republish banner, confirmation summary, error messages, success confirmation), so the Pitfall 14 jargon audit is enforced HERE. Use the `showDEDialog` reactive UI pattern as the analog — same shape (ChoiceInput + threshold-like inputs + dialog `.onOK` callback).

Output: `src/publishing/share-dialog.ts` exporting `showShareForReviewDialog(df): Promise<void>`. Plan 07 wires the menu entry to it via `await`.
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/publishing/publish-state.ts
@packages/Proteomics/src/publishing/publish-project.ts
@packages/Proteomics/src/analysis/differential-expression.ts

<interfaces>
<!-- Pulled from showDEDialog (closest dialog analog) -->

From src/analysis/differential-expression.ts (showDEDialog at lines 282-364):
  ui.input.choice(label, opts): DG.InputBase
  ui.input.string(label, opts): DG.InputBase
  ui.input.textArea(label, opts): DG.InputBase   // available, takes {value, ...}
  ui.div(): HTMLDivElement
  ui.dialog(title): DG.Dialog
  dialog.add(node).add(node).onOK(callback).show()
  inputBase.onChanged.subscribe(callback)

From Plan 01 (publish-state.ts):
  PublishOptions, findPriorShare(target, group), slugifyTarget(raw), buildMailtoUrl(opts), DE_COMPLETE_TAG, MailtoOptions

From Plan 04 (publish-project.ts):
  publishAnalysis(df, opts): Promise<DG.Project>

From datagrok-api:
  grok.dapi.groups.list(): Promise<DG.Group[]>
  grok.shell.user.group: DG.Group
  grok.shell.user.email: string | null
  grok.shell.info(msg), grok.shell.warning(msg), grok.shell.error(msg)
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Build async dialog skeleton + group ChoiceInput + republish detection + reactive confirmation summary (W-4)</name>
  <files>src/publishing/share-dialog.ts</files>
  <read_first>
    - @packages/Proteomics/src/analysis/differential-expression.ts (lines 282-364 — showDEDialog full structure)
    - @packages/Proteomics/CLAUDE.md (function-naming prefix `showX(df)` + dialog suffix `...`)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-02 reviewer-group sourcing; D-04 republish UX; Pitfall 14 jargon audit; D-domain biologist-consumer pin; PUB-08 confirmation summary mandate)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 4 — analog + notes; "Audience pin" banned-word list)
    - @packages/Proteomics/src/publishing/publish-state.ts (PublishOptions, findPriorShare, slugifyTarget, DE_COMPLETE_TAG, buildMailtoUrl)
  </read_first>
  <action>
Create `src/publishing/share-dialog.ts`. Standard imports plus `findPriorShare`, `slugifyTarget`, `PublishOptions`, `DE_COMPLETE_TAG`, `buildMailtoUrl`, `MailtoOptions` from `./publish-state` and `publishAnalysis` from `./publish-project`.

**W-6 fix:** Export `showShareForReviewDialog(df: DG.DataFrame): Promise<void>` — note the `async` and `Promise<void>` return type. (Function-naming `showX(df)` per CLAUDE.md.)

**I-10 fix:** Precondition uses the constant: `if (df.getTag(DE_COMPLETE_TAG) !== 'true')` — NOT inline `'proteomics.de_complete'`.

**Skeleton:**
- Validate precondition: `if (df.getTag(DE_COMPLETE_TAG) !== 'true')` -> `grok.shell.warning('Run Differential Expression first.'); return;` (defensive — the menu handler in Plan 07 also gates).

- Load groups list: `const allGroups = await grok.dapi.groups.list();`. Filter to groups the current user can administer:
  - Simplest filter that still resolves Spike A2 / D-02 intent: `allGroups.filter(g => g.adminMemberships.includes(grok.shell.user.id))` — verify this attribute exists on `DG.Group`; if not, fall back to `allGroups.filter(g => g.name !== 'All users')` (excludes the world group; biologist-scoped groups appear).
  - If filter returns empty: `grok.shell.warning('No reviewer groups available. Ask an admin to create a reviewer group.'); return;`

- Reviewer group ChoiceInput:
  ```
  const groupInput = ui.input.choice('Share with team', {
    items: filteredGroups.map(g => g.friendlyName),
    nullable: false,
  });
  ```
  Use string items (group names) NOT `DG.Group` objects (ChoiceInput needs primitive values). Keep a lookup: `const groupByName = new Map(filteredGroups.map(g => [g.friendlyName, g]));`. On OK, resolve `groupByName.get(groupInput.value)`.

  Note: "Share with team" not "Reviewer Group" — biologist jargon audit (Pitfall 14 banned: `Reviewer Group` reads as IT-speak; "team" is biologist-natural).

- Target text input:
  ```
  const targetInput = ui.input.string('Target', {value: ''});
  targetInput.setTooltip('Free text — your team\\'s name for this target (e.g., MYH7-DMD, Cytokinetics atrophy panel)');
  ```

- Note text area:
  ```
  const noteInput = ui.input.textArea('Note for reviewers', {value: ''});
  ```

- Republish-detection banner (initially hidden):
  ```
  const banner = ui.div();
  banner.style.cssText = 'background:#fff3cd; padding:8px; border-radius:4px; margin-bottom:8px; display:none;';
  ```

- **W-4 — Reactive Confirmation Summary block (PUB-08):**
  ```
  const summary = ui.div();
  summary.style.cssText = 'background:#f0f0f0; padding:8px; border-radius:4px; margin:8px 0; font-family: monospace; font-size:11px;';
  ```
  The summary div is ALWAYS visible (not hidden behind banner condition). It refreshes on every input change. It displays four lines:
  - `Project: Proteomics-Review-<slug>-v<N>-<YYYY-MM-DD>`
  - `Will be visible to: <reviewer group friendly name>`
  - `Status: <New share | Will supersede '<prior name>' (v<N+1>)>`
  - `Date: <YYYY-MM-DD>` (today)
  
  This is the PUB-08 "confirmation summary before creating the Project" — the user sees exactly what's about to happen before clicking OK. The republish banner alone (which only appears for republish) is NOT a substitute.

- Hold the current `priorVersion` in a closure-scoped `let priorVersion: DG.Project | null = null;` variable.

- Reactive update function (updates BOTH banner AND summary):
  ```
  const updateBannerAndSummary = async () => {
    const target = targetInput.value?.trim() ?? '';
    const group = groupByName.get(groupInput.value);
    const slug = target ? slugifyTarget(target) : '<no-target>';
    const dateStr = new Date().toISOString().slice(0, 10);
    
    // Republish detection (drives both banner display + summary "Status" line)
    if (target && group) {
      priorVersion = await findPriorShare(target, group);
    } else {
      priorVersion = null;
    }
    
    const nextVersion = priorVersion ? <parse v from priorVersion.name regex> + 1 : 1;
    
    // Banner
    if (priorVersion) {
      banner.textContent = '⚠ This will share as v' + nextVersion + ' and supersede "' + priorVersion.name + '". The previous version stays available for reference.';
      banner.style.display = 'block';
    } else {
      banner.style.display = 'none';
    }
    
    // Confirmation summary (W-4 — always visible, refreshes on every change)
    const projectName = 'Proteomics-Review-' + slug + '-v' + nextVersion + '-' + dateStr;
    const groupLabel = group ? group.friendlyName : '<pick a team>';
    const statusLabel = priorVersion 
      ? 'Will supersede "' + priorVersion.name + '" (v' + nextVersion + ')'
      : 'New share';
    summary.innerHTML = ''; // clear
    summary.appendChild(ui.divText('Project: ' + projectName));
    summary.appendChild(ui.divText('Will be visible to: ' + groupLabel));
    summary.appendChild(ui.divText('Status: ' + statusLabel));
    summary.appendChild(ui.divText('Date: ' + dateStr));
  };
  
  targetInput.onChanged.subscribe(updateBannerAndSummary);
  groupInput.onChanged.subscribe(updateBannerAndSummary);
  
  // Initial render so the summary shows up immediately on dialog open (with placeholders).
  await updateBannerAndSummary();
  ```

Per saved memory `feedback_datagrok_dialog_dom.md`: dialog buttons use row-reverse layout. Use `ui.dialog('Share Analysis for Review').add(targetInput).add(groupInput).add(noteInput).add(banner).add(summary).onOK(handler).show();`

Add JSDoc note: "Pitfall 14: all reviewer-touchable strings audited — banned `DataFrame`, `tag`, `semType`, `ACL`, `viewer factory`. Replacements: 'table', 'label', 'who can see', 'view'. PUB-08 confirmation summary refreshes reactively (W-4)."
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/share-dialog.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export async function showShareForReviewDialog" src/publishing/share-dialog.ts | grep -q '^1$' &amp;&amp; grep -c "Promise<void>" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; grep -c "DE_COMPLETE_TAG" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; grep -c "buildMailtoUrl" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; grep -c "findPriorShare" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; grep -c "groups\\.list" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; grep -c "summary\\.appendChild\\|confirmation summary" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; { ! grep -E "DataFrame|semType|ACL|viewer factory" src/publishing/share-dialog.ts | grep -v "^//" | grep -v "import"; }</automated>
  </verify>
  <done>Dialog opens; signature is async returning Promise<void> (W-6); loads groups via `dapi.groups.list()` filtered to administrable groups; renders target string + reviewer-group ChoiceInput + note textArea + republish banner + ALWAYS-VISIBLE reactive confirmation summary (W-4); banner + summary update reactively on target/group change; precondition uses DE_COMPLETE_TAG constant (I-10); buildMailtoUrl is imported from publish-state (B-1); banned jargon words absent in user-facing strings; strict TypeScript compiles.</done>
</task>

<task type="auto" priority="P2" tdd="false">
  <name>Task 2 (P2 - PUB-13 mailto sub-feature): Implement onOK handler — call publishAnalysis + surface errors + success state with mailto</name>
  <files>src/publishing/share-dialog.ts</files>
  <read_first>
    - @packages/Proteomics/src/publishing/publish-project.ts (publishAnalysis signature + error contracts — exact D-03 string)
    - @packages/Proteomics/src/publishing/publish-state.ts (buildMailtoUrl signature + MailtoOptions)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (Claude's discretion: mailto subject/body shape; PUB-13 mailto)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Open Question 5" — sharer email source + null fallback; §"Don't Hand-Roll" — mailto via encodeURIComponent)
    - @packages/Proteomics/src/analysis/differential-expression.ts (lines 379-441 — try/catch + pi.close finally + grok.shell.error pattern)
  </read_first>
  <action>
**B-3 task-level priority:** This task is P2 because it covers PUB-13 (mailto) at the success-state UI. The publishAnalysis call itself is P1 (PUB-05); only the post-success mailto link rendering is the P2 portion. If P2 must be deferred, the task can ship without the mailto link in the success modal — but the rest of the onOK path (publishAnalysis call + success state + error surface) ships as P1.

Add the onOK handler to the dialog created in Task 1. **Note: `buildMailtoUrl` is IMPORTED from `./publish-state` — it is NOT defined here per B-1 fix.**

Implement the `onOK` callback for the dialog:
```
.onOK(async () => {
  const target = targetInput.value?.trim();
  if (!target) { grok.shell.warning('Target is required.'); return; }
  const group = groupByName.get(groupInput.value);
  if (!group) { grok.shell.warning('Pick a team to share with.'); return; }

  const opts: PublishOptions = {
    target,
    reviewerGroup: group,
    note: noteInput.value ?? '',
    priorVersion,
  };

  try {
    const project = await publishAnalysis(df, opts);
    // SUCCESS STATE
    // P2 (PUB-13): mailto link in success state
    const mailtoUrl = buildMailtoUrl({
      sharerEmail: grok.shell.user.email,
      sharerName: grok.shell.user.friendlyName,
      projectName: project.name,
      publishedDateStr: new Date().toISOString().slice(0, 10),
    });
    grok.shell.info('Shared: ' + project.name);
    const successBody = ui.div([
      ui.divText('Shared as ' + project.name + ' with team ' + group.friendlyName + '.'),
      ui.link('Open shared analysis', () => project.open()),
      ui.link('Send request to re-run', mailtoUrl, undefined, undefined), // P2 - PUB-13 mailto
    ]);
    ui.dialog('Shared Successfully').add(successBody).showModal(false);
  } catch (e: any) {
    // Surface error VERBATIM — Plan 04 step 7 throws the D-03 exact string
    grok.shell.error('Share failed: ' + (e?.message ?? String(e)));
  }
})
```

Surface the success popup distinctly from the error path. Per saved memory `feedback_datagrok_dialog_dom.md`: dialog buttons row-reverse; modal vs non-modal is `showModal(true|false)` — pick `false` so reviewer can keep working.

Per CONTEXT.md "audience pin" — every string in this handler avoids `DataFrame`, `tag`, `semType`, `ACL`, `viewer factory`. Replacements applied: "Shared" not "Published" if appropriate, "Send request to re-run" not "Submit reanalysis request", etc.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/share-dialog.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; { ! grep -E "^export function buildMailtoUrl" src/publishing/share-dialog.ts; } &amp;&amp; grep -c "buildMailtoUrl" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; grep -c "publishAnalysis" src/publishing/share-dialog.ts | grep -qv '^0$' &amp;&amp; grep -c "grok\\.shell\\.error" src/publishing/share-dialog.ts | grep -qv '^0$'</automated>
  </verify>
  <done>onOK handler validates inputs, builds `PublishOptions`, calls `publishAnalysis`, surfaces success modal with project-open link + mailto link (P2 sub-feature using imported `buildMailtoUrl`), surfaces errors verbatim via `grok.shell.error`; `buildMailtoUrl` is NOT exported from this file (B-1 fix — it lives in publish-state.ts); strict TypeScript compiles.</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| analyst user input (target, note) -> publishAnalysis | Freeform input; sanitized inside `slugifyTarget` (Plan 01) |
| analyst-selected reviewer group -> ACL grant | Group object retrieved from server `dapi.groups.list()`; cannot be spoofed client-side |
| sharer email -> mailto URL | Pure client construction; no network call; `encodeURIComponent` in buildMailtoUrl prevents header injection |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-03 | Tampering (target injection at dialog input) | targetInput | mitigate | Trim + slugifyTarget downstream in publishAnalysis; raw target preserved in tag/column only |
| T-15-SP-06 | Tampering (mailto header injection via subject/body) | buildMailtoUrl (in publish-state.ts) | mitigate | All variable parts wrapped in `encodeURIComponent`; subject/body templates use only sanitized strings (Plan 01 owns the helper) |
| T-15-SP-07 | Information disclosure (group list leak) | grok.dapi.groups.list filter | mitigate | Filter to groups user can administer; minimizes the visible group population to those legitimately shareable |
</threat_model>

<verification>
- TypeScript strict-mode compiles
- Dialog opens against the running server (smoke check)
- Reactive confirmation summary appears on dialog open with placeholders, refreshes as user edits inputs (W-4 PUB-08)
- Republish banner appears when target+group match a prior published Project
- Error path surfaces verbatim D-03 string when verify-and-rollback fires
- Mailto link opens user's mail client with subject/body prefilled (P2 sub-feature)
- `buildMailtoUrl` is NOT exported from this file (B-1 — owned by publish-state.ts)
- `showShareForReviewDialog` returns `Promise<void>` and is callable via `await` (W-6)
- Precondition uses `DE_COMPLETE_TAG` constant (I-10)
</verification>

<success_criteria>
- `src/publishing/share-dialog.ts` exports ONLY `showShareForReviewDialog(df): Promise<void>` (NOT `buildMailtoUrl` — that lives in publish-state.ts per B-1)
- Function signature is async returning Promise<void> (W-6)
- Reviewer-group ChoiceInput sources from `grok.dapi.groups.list()`, filtered to administrable groups
- Republish banner reactively detects via `findPriorShare`
- Reactive confirmation summary always visible (PUB-08, W-4)
- Success state offers Open Project + Send Request mailto links (PUB-13 P2 sub-feature)
- All user-facing strings pass the Pitfall 14 jargon audit
- Plan-level priority is NOT P2 (B-3 — plan covers P1 requirements; only Task 2 sub-feature for PUB-13 is P2)
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-05-SUMMARY.md` when done with:
- Confirmation that all reviewer-touchable strings are biologist-jargon-clean
- Notes on which `DG.Group` attribute was used to filter to administrable groups
- Mailto-link template verified against an open mail client (smoke)
- Confirmation that `buildMailtoUrl` is imported from publish-state.ts (B-1)
- Confirmation that function signature is async / Promise<void> (W-6)
- Confirmation that PUB-08 confirmation summary is reactive (W-4)
</output>
