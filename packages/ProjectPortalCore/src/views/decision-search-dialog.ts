/** Searching gate decisions: by text (rationale, outcome, evidence), program, gate and portfolio,
 * across the portfolios the user can see. A result opens its program. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {decisionsOf} from '../decision-records';
import {DecisionQuery, DecisionView, describeDecision, searchDecisions} from '../model/decision-records';
import {feedSelect} from './feed-dialog';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

/** Results shown at once. */
const SHOWN = 50;

export interface DecisionSearchOptions {
  portfolios: {code: string, name: string}[];
  /** The portfolio open on the board: searched first. */
  current: string | null;
  openProgram: (portfolioCode: string, programId: string) => void | Promise<void>;
}

export async function openDecisionSearchDialog(o: DecisionSearchOptions): Promise<DG.Dialog> {
  installStyles();
  const dlg = ui.dialog({title: 'Decisions'});
  dlg.root.dataset.u2Name = 'decision-search-dialog';
  const q: DecisionQuery = {portfolio: o.current ?? '', withoutSync: true};
  const text = (name: string, placeholder: string, key: 'text' | 'program' | 'gate') => {
    const el = named(document.createElement('input'), name);
    el.className = 'pp-text';
    el.placeholder = placeholder;
    let timer: any = null;
    el.oninput = () => {
      clearTimeout(timer);
      timer = setTimeout(() => {
        q[key] = el.value;
        draw();
      }, 150);
    };
    return el;
  };
  const box = (name: string, label: string, key: 'withoutSync' | 'withEvidence') => {
    const el = named(document.createElement('input'), name);
    el.type = 'checkbox';
    el.checked = !!q[key];
    el.onchange = () => {
      q[key] = el.checked;
      draw();
    };
    const l = document.createElement('label');
    l.className = 'pp-evidence-pick';
    l.append(el, span(label));
    return l;
  };
  const results = named(divV([span('Reading decisions…', 'pp-muted')], 'pp-units'), 'decision-results');
  const summary = named(span('', 'pp-muted'), 'decision-summary');
  dlg.add(divV([
    divH([text('decision-text-input', 'Words in the rationale, outcome or evidence', 'text')], 'pp-actions'),
    divH([text('decision-program-input', 'Program code or name', 'program'), text('decision-gate-input', 'Gate (stage code or name)', 'gate')],
      'pp-actions'),
    divH([feedSelect('decision-portfolio-input', 'Portfolio', [{code: '', name: 'All I can see'}, ...o.portfolios],
      q.portfolio ?? '', (v) => {
        q.portfolio = v;
        draw();
      }), box('decision-without-sync', 'Leave out registry moves', 'withoutSync'),
    box('decision-with-evidence', 'Only with evidence', 'withEvidence')], 'pp-actions'),
    summary,
    results,
  ], 'pp-feed'));
  dlg.show({width: 760});

  let all: DecisionView[] = [];
  const draw = () => {
    const found = searchDecisions(all, {...q, portfolio: q.portfolio || undefined});
    summary.textContent = `${found.length} of ${all.length} decision${all.length === 1 ? '' : 's'}` +
      (found.length > SHOWN ? `; the newest ${SHOWN} shown` : '');
    results.replaceChildren(...(found.length ? found.slice(0, SHOWN).map((d, i) => {
      const line = named(span(describeDecision(d), 'pp-link'), `decision-hit-${i}`);
      line.title = `Open ${d.programCode} in ${d.portfolio}`;
      line.onclick = () => {
        dlg.close();
        void Promise.resolve(o.openProgram(d.portfolio, d.programId)).catch((x) => grok.shell.error(x?.message ?? String(x)));
      };
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
