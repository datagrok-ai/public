/** Searching gate decisions: by text (rationale, outcome, evidence), program, gate and portfolio,
 * across the portfolios the user can see. A result opens its program. */
import * as grok from 'datagrok-api/grok';
import {divH, divV, loader, span} from '@datagrok-libraries/u2';
import {decisionsOf} from '../decision-records';
import {DecisionQuery, DecisionView, describeDecision, searchDecisions} from '../model/decision-records';
import {feedSelect} from './feed-dialog';
import {check, linkTo, named, openDialog, PortalDialog, text as textInput} from './u2-controls';

/** Results shown at once. */
const SHOWN = 50;

export interface DecisionSearchOptions {
  portfolios: {code: string, name: string}[];
  /** The portfolio open on the board: searched first. */
  current: string | null;
  openProgram: (portfolioCode: string, programId: string) => void | Promise<void>;
}

export async function openDecisionSearchDialog(o: DecisionSearchOptions): Promise<PortalDialog> {
  const q: DecisionQuery = {portfolio: o.current ?? '', withoutSync: true};
  const text = (name: string, placeholder: string, key: 'text' | 'program' | 'gate') => {
    let timer: any = null;
    const input = textInput(name, undefined, {placeholder, inline: true});
    input.effect(() => {
      const v = input.value.value;
      clearTimeout(timer);
      timer = setTimeout(() => {
        if (q[key] === v || (!q[key] && !v))
          return;
        q[key] = v;
        draw();
      }, 150);
    });
    input.own(() => clearTimeout(timer));
    return input.root;
  };
  const box = (name: string, label: string, key: 'withoutSync' | 'withEvidence') =>
    check(name, label, !!q[key], {onChanged: (v) => {
      q[key] = v;
      draw();
    }}).root;
  const results = named(divV([loader('Reading decisions…')], 'pp-units'), 'decision-results');
  const summary = named(span('', 'pp-muted'), 'decision-summary');
  const dlg = openDialog({title: 'Decisions', name: 'decision-search-dialog', width: 760}, () => ({content: divV([
    divH([text('decision-text-input', 'Words in the rationale, outcome or evidence', 'text')], 'pp-actions pp-fill'),
    divH([text('decision-program-input', 'Program code or name', 'program'),
      text('decision-gate-input', 'Gate (stage code or name)', 'gate')], 'pp-actions pp-fill'),
    divH([feedSelect('decision-portfolio-input', 'Portfolio', [{code: '', name: 'All I can see'}, ...o.portfolios],
      q.portfolio ?? '', (v) => {
        q.portfolio = v;
        draw();
      }), box('decision-without-sync', 'Leave out registry moves', 'withoutSync'),
    box('decision-with-evidence', 'Only with evidence', 'withEvidence')], 'pp-actions'),
    summary,
    results,
  ], 'pp-feed')}));

  let all: DecisionView[] = [];
  const draw = () => {
    const found = searchDecisions(all, {...q, portfolio: q.portfolio || undefined});
    summary.textContent = `${found.length} of ${all.length} decision${all.length === 1 ? '' : 's'}` +
      (found.length > SHOWN ? `; the newest ${SHOWN} shown` : '');
    results.replaceChildren(...(found.length ? found.slice(0, SHOWN).map((d, i) => {
      const line = linkTo(describeDecision(d), `decision-hit-${i}`, () => {
        dlg.close();
        void Promise.resolve(o.openProgram(d.portfolio, d.programId)).catch((x) => grok.shell.error(x?.message ?? String(x)));
      }, {title: `Open ${d.programCode} in ${d.portfolio}`});
      const ev = d.evidence.length ? [named(span(`Evidence: ${d.evidence.map((e) => e.name).join('; ')}`, 'pp-muted'),
        `decision-hit-evidence-${i}`)] : [];
      return named(divV([line, ...ev], 'pp-unit'), `decision-${i}`);
    }) : [named(span('No decision matches', 'pp-muted'), 'decision-none')]));
  };
  try {
    all = await decisionsOf();
    draw();
  } catch (e: any) {
    results.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'decision-error'));
  }
  return dlg;
}
