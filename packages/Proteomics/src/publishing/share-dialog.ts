import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DE_COMPLETE_TAG, PublishOptions, findPriorShare, slugifyTarget,
} from './publish-state';
import {publishAnalysis} from './publish-project';
import {reviewNamePrefix, verifyPublishedDashboard} from './publish-settings';

/**
 * Analyst-facing share dialog. Single entry point into the publish flow —
 * every reviewer-touchable string flows through here (target label, group
 * label, note placeholder, republish banner, confirmation summary, error
 * messages, success state). The Pitfall 14 biologist-jargon audit is
 * enforced here.
 *
 * Async (`Promise<void>`) per W-6 — loads reviewer groups via
 * `dapi.groups.list()` before composing the dialog. Plan 07 menu handler
 * must `await` (or `void`-fire) it.
 *
 * - Republish-detection banner reactively driven by {@link findPriorShare}
 * - PUB-08 reactive confirmation summary (W-4): always-visible block that
 *   updates on every input change so the analyst sees the exact resulting
 *   Project name / team / supersede status BEFORE clicking OK
 * - Mailto link in the success state (PUB-13, P2) uses {@link buildMailtoUrl}
 *   imported from `./publish-state` (NOT redefined here per B-1)
 */
export async function showShareForReviewDialog(df: DG.DataFrame): Promise<void> {
  if (df.getTag(DE_COMPLETE_TAG) !== 'true') {
    grok.shell.warning('Run Differential Expression first.');
    return;
  }

  const allGroups = await grok.dapi.groups.list();
  const allUsersId = (DG.Group as any).defaultGroupsIds?.['All users'];
  const filteredGroups = (allGroups ?? []).filter((g: any) => {
    if (g == null) return false;
    if (g.hidden) return false;
    if (g.personal) return false;
    if (allUsersId && g.id === allUsersId) return false;
    return true;
  });

  if (filteredGroups.length === 0) {
    grok.shell.warning('No teams available to share with. Ask an admin to create a reviewer team.');
    return;
  }

  const groupByName = new Map<string, DG.Group>();
  for (const g of filteredGroups) {
    const name = (g as any).friendlyName ?? (g as any).name ?? '';
    if (name && !groupByName.has(name)) groupByName.set(name, g);
  }
  const groupItems = Array.from(groupByName.keys());

  const targetInput = ui.input.string('Target', {value: ''});
  targetInput.setTooltip('Free text — your team\'s name for this target (e.g., MYH7-DMD, Cytokinetics atrophy panel)');

  const groupInput = ui.input.choice('Share with team', {
    value: groupItems[0],
    items: groupItems,
    nullable: false,
  });

  const noteInput = (ui.input as any).textArea
    ? (ui.input as any).textArea('Note for reviewers', {value: ''})
    : ui.input.string('Note for reviewers', {value: ''});

  // Per-share verification toggle. Defaults to the package setting
  // (`verifyPublishedDashboard`), but lets the user skip the heavy reopen-and-check
  // round trip for a fast demo share without changing the deployment default.
  const verifyInput = ui.input.bool('Verify published dashboard', {value: verifyPublishedDashboard()});
  verifyInput.setTooltip(
    'Re-open the shared dashboard after publishing to confirm it survives a reload. ' +
    'On = safer (recommended for client deliverables); off = faster share.');

  const banner = ui.div();
  banner.style.cssText = 'background:#fff3cd; padding:8px; border-radius:4px; margin-bottom:8px; display:none;';

  const summary = ui.div();
  summary.style.cssText = 'background:#f0f0f0; padding:8px; border-radius:4px; margin:8px 0; font-family: monospace; font-size:11px;';

  let priorVersion: DG.Project | null = null;

  const parsePriorVersionN = (priorName: string | null | undefined): number => {
    if (!priorName) return 0;
    const m = /-v(\d+)-/.exec(priorName);
    if (!m) return 0;
    const n = parseInt(m[1], 10);
    return Number.isFinite(n) ? n : 0;
  };

  const updateBannerAndSummary = async (): Promise<void> => {
    const target = (targetInput.value ?? '').trim();
    const groupName = groupInput.value;
    const group = groupName ? groupByName.get(groupName) ?? null : null;
    const slug = target ? slugifyTarget(target) : '<no-target>';
    const dateStr = new Date().toISOString().slice(0, 10);

    if (target && group) {
      try { priorVersion = await findPriorShare(target, group); }
      catch { priorVersion = null; }
    } else {
      priorVersion = null;
    }

    const priorN = parsePriorVersionN(priorVersion ? (priorVersion as any).name : null);
    const nextVersion = priorN > 0 ? priorN + 1 : 1;

    if (priorVersion) {
      banner.textContent =
        `⚠ This will share as v${nextVersion} and supersede "${(priorVersion as any).name}". ` +
        `The previous version stays available for reference.`;
      banner.style.display = 'block';
    } else {
      banner.style.display = 'none';
    }

    const projectName = `${reviewNamePrefix()}-${slug}-v${nextVersion}-${dateStr}`;
    const groupLabel = group ? ((group as any).friendlyName ?? (group as any).name ?? '<unnamed team>') : '<pick a team>';
    const statusLabel = priorVersion
      ? `Will supersede "${(priorVersion as any).name}" (v${nextVersion})`
      : 'New share';

    summary.innerHTML = '';
    summary.appendChild(ui.divText(`Project: ${projectName}`));
    summary.appendChild(ui.divText(`Will be visible to: ${groupLabel}`));
    summary.appendChild(ui.divText(`Status: ${statusLabel}`));
    summary.appendChild(ui.divText(`Date: ${dateStr}`));
  };

  targetInput.onChanged.subscribe(() => { void updateBannerAndSummary(); });
  groupInput.onChanged.subscribe(() => { void updateBannerAndSummary(); });

  await updateBannerAndSummary();

  ui.dialog('Share Analysis for Review')
    .add(targetInput)
    .add(groupInput)
    .add(noteInput)
    .add(verifyInput)
    .add(banner)
    .add(summary)
    .onOK(async () => {
      const target = (targetInput.value ?? '').trim();
      if (!target) { grok.shell.warning('Target is required.'); return; }
      const group = groupInput.value ? groupByName.get(groupInput.value) : null;
      if (!group) { grok.shell.warning('Pick a team to share with.'); return; }

      const opts: PublishOptions = {
        target,
        reviewerGroup: group,
        note: (noteInput.value ?? '') as string,
        priorVersion,
        verify: verifyInput.value ?? true,
      };

      try {
        const project = await publishAnalysis(df, opts);
        const groupName = (group as any).friendlyName ?? (group as any).name ?? '(unnamed team)';
        // Single non-modal confirmation — the published project is browsable in the
        // Spaces tree, so no modal or action links are needed. publishAnalysis no
        // longer emits its own toast, so this is the one and only confirmation.
        grok.shell.info(`Shared as ${(project as any).name} with team ${groupName}.`);
      } catch (e: any) {
        grok.shell.error(`Share failed: ${e?.message ?? String(e)}`);
      }
    })
    .show();
}
