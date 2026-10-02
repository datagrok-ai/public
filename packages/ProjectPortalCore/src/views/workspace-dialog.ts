/** A portfolio's workspace in the app: the orphan report (root folders no program claims, folders
 * of programs that ended), **Preview park** and **Park** for approvers, and the parked folders with
 * **Unpark** each. Opened from the board's Workspace button. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {badge, button, divH, divV, span} from '@datagrok-libraries/u2';
import {PortfolioRow} from '../generated/db';
import {isApprover} from '../portal';
import {workspaceOrphans} from '../workspace';
import {parkedFoldersOf, parkOrphans, summarizeParkResult, unparkFolder} from '../parking';
import {installStyles} from './styles';

const slug = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

/** The panel: drawn again after each park or unpark. */
export async function workspacePanel(portfolio: PortfolioRow): Promise<HTMLElement> {
  installStyles();
  const root = named(divV([], 'pp-sync-history'), 'workspace-panel');
  const result = named(divV([], 'pp-units'), 'park-result');
  const approver = await isApprover(portfolio).catch(() => false);
  const failed = (e: any) => result.replaceChildren(span(e?.message ?? String(e), 'pp-reason'));
  const lines = (text: string) => text.split('\n').map((l) => span(l.trim(), 'pp-muted'));
  const draw = async () => {
    const report = await workspaceOrphans(portfolio.code);
    const parked = await parkedFoldersOf(portfolio);
    const preview = named(button('Preview park', () => {
      result.replaceChildren(span('Planning…', 'pp-muted'));
      parkOrphans(portfolio.code, true).then((r) => result.replaceChildren(...lines(summarizeParkResult(r)))).catch(failed);
    }), 'park-preview-button');
    preview.title = 'What parking would move, and where; nothing changes';
    preview.disabled = report.orphans.length === 0;
    const parkSome = (only?: string[]) => {
      result.replaceChildren(span('Parking…', 'pp-muted'));
      parkOrphans(portfolio.code, false, only).then(async (r) => {
        await draw();
        result.replaceChildren(...lines(summarizeParkResult(r)));
      }).catch(failed);
    };
    const park = named(button(`Park all ${report.orphans.length}`, () => parkSome(), {primary: true}), 'park-button');
    park.disabled = !approver || report.orphans.length === 0;
    park.title = approver ? `Move them into _Parked/<today> under ${report.root}; each can be unparked. Nothing is deleted` :
      'Only the portfolio\'s approvers park folders';
    const orphans = report.orphans.map((o) => {
      const key = slug(o.path.split('/').pop() ?? o.id);
      const one = named(span('Park', approver ? 'pp-link' : 'pp-muted'), `park-orphan-${key}`);
      one.title = approver ? `Move it into _Parked/<today> under ${report.root}` : 'Only the portfolio\'s approvers park folders';
      if (approver)
        one.onclick = () => parkSome([o.id]);
      return named(divH([badge(o.kind === 'park' ? 'program ended' : 'unclaimed', {variant: o.kind === 'park' ? 'warning' : 'default'}),
        span(o.path), span(o.reason, 'pp-muted'), one], 'pp-unit'), `orphan-${key}`);
    });
    const parkedRows = parked.map((p) => {
      const key = slug(p.name);
      const unpark = named(button('Unpark', () => {
        unparkFolder(p.id).then(async (r) => {
          await draw();
          result.replaceChildren(...lines([r.line, ...r.notes].join('\n')));
        }).catch(failed);
      }), `unpark-${key}`);
      unpark.disabled = !approver;
      unpark.title = approver ? `Move it back to ${p.from_path}` : 'Only the portfolio\'s approvers unpark folders';
      return named(divH([span(p.parked_path), span(p.reason ?? '', 'pp-muted'), unpark], 'pp-unit'), `parked-${key}`);
    });
    root.replaceChildren(
      named(span(`Workspace root: ${report.root}`), 'workspace-root'),
      ...report.notes.map((n) => span(n, 'pp-muted')),
      span(`Orphan folders (${report.orphans.length})`, 'pp-heading'),
      report.orphans.length ? named(divV(orphans, 'pp-units'), 'orphans') :
        named(span('No folder for prune-or-park to look at', 'pp-muted'), 'no-orphans'),
      named(divH([preview, park], 'pp-actions'), 'park-actions'),
      span('Parking moves a folder into _Parked/<date> under the root and records it. Nothing is deleted: ' +
        'removing a parked folder is a person\'s call, by hand.', 'pp-muted'),
      result,
      span(`Parked folders (${parked.length})`, 'pp-heading'),
      parked.length ? named(divV(parkedRows, 'pp-units'), 'parked-folders') :
        named(span('Nothing parked', 'pp-muted'), 'no-parked'),
    );
  };
  await draw();
  return root;
}

/** Opens the workspace panel of [portfolio] in a dialog. */
export async function openWorkspaceDialog(portfolio: PortfolioRow): Promise<DG.Dialog> {
  const panel = await workspacePanel(portfolio);
  const dlg = ui.dialog({title: `${portfolio.name}: workspace`}).add(panel);
  dlg.root.dataset.u2Name = 'workspace-dialog';
  dlg.show({width: 620});
  return dlg;
}
