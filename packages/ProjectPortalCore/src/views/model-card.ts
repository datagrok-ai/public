/** The model card (F-12, FR-50): the scientific question first, then what the model is for, its
 * inputs with their validity range, its outputs, applicability, performance, lineage and approval.
 * One renderer serves the program page's saved runs, the `modelCard` function and the info panel. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {badge, divH, divV, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {CardParam} from '../model/catalogue';
import {ModelCardView, modelCardOf} from '../catalogue';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const slug = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

const approvalVariant: Record<string, BadgeVariant> = {
  draft: 'warning', in_review: 'accent', approved: 'success', retired: 'error',
};
const statusVariant: Record<string, BadgeVariant> = {active: 'success', experimental: 'warning', deprecated: 'error'};

const range = (p: CardParam): string =>
  p.min === undefined && p.max === undefined ? '' : `${p.min ?? '…'} – ${p.max ?? '…'}`;

function paramTable(params: CardParam[], kind: 'input' | 'output'): HTMLElement {
  const table = document.createElement('table');
  table.className = 'pp-card-table';
  const head = table.createTHead().insertRow();
  for (const h of kind === 'input' ? ['Input', 'Units', 'Valid for'] : ['Output', 'Units'])
    head.append(Object.assign(document.createElement('th'), {textContent: h}));
  const body = table.createTBody();
  for (const p of params) {
    const row = body.insertRow();
    row.dataset.u2Name = `model-card-${kind}-${slug(p.name)}`;
    row.insertCell().textContent = p.caption ?? p.name;
    row.insertCell().textContent = p.units ?? '';
    if (kind === 'input')
      row.insertCell().textContent = range(p);
    if (p.description)
      row.title = p.description;
  }
  return named(table, `model-card-${kind}s`);
}

/** A labelled section, left out when there is nothing to say. */
function section(name: string, title: string, text: string | null | undefined): HTMLElement[] {
  if (!text)
    return [];
  return [span(title, 'pp-heading'), named(span(text, 'pp-card-text'), `model-card-${name}`)];
}

/** The card of a catalogued model. */
export function modelCardElement(v: ModelCardView): HTMLElement {
  installStyles();
  const {model, card} = v;
  const header = divH([
    badge(model.status ?? 'active', {variant: statusVariant[model.status ?? 'active'] ?? 'default'}),
    ...(card ? [named(badge(card.approval.replace('_', ' '), {variant: approvalVariant[card.approval] ?? 'default'}),
      'model-card-approval'), badge(`v${card.model_version}`, {variant: 'accent'})] : []),
    ...(model.category ? [span(model.category, 'pp-muted')] : []),
  ], 'pp-card-meta');
  const lineageNote = [card?.lineage, card?.supersedes ? `Replaces version ${card.supersedes}.` : null]
    .filter((x) => x).join(' ');
  const facts: [string, string][] = [
    ['Owner', card?.owner ?? ''],
    ['Unit operations', v.unitOperations.map((u) => u.name).join(', ')],
    ['Default run type', v.defaultRunType?.name ?? ''],
    ['Portfolios', v.portfolios.join(', ')],
    ['Approved', card?.approved_on ? card.approved_on.format('YYYY-MM-DD') : ''],
  ];
  const others = v.versions.filter((x) => x.id !== card?.id);
  return named(divV([
    named(span(model.name, 'pp-title'), 'model-card-title'),
    span(model.code, 'pp-muted'),
    header,
    ...(!v.installed ? [named(span('The model\'s function is not installed on this instance', 'pp-reason'),
      'model-card-not-installed')] : []),
    ...(!card ? [named(span('No card yet: the catalogue declares the model without one', 'pp-muted'), 'model-card-empty')] : []),
    ...section('question', 'Question', card?.question),
    ...section('intended-use', 'Intended use', card?.intended_use ?? model.description),
    ...(v.inputs.length ? [span('Inputs', 'pp-heading'), paramTable(v.inputs, 'input')] : []),
    ...(v.outputs.length ? [span('Outputs', 'pp-heading'), paramTable(v.outputs, 'output')] : []),
    ...section('applicability', 'Applicability', card?.applicability),
    ...section('performance', 'Performance', card?.performance),
    ...section('lineage', 'Lineage', lineageNote),
    ...section('limitations', 'Limitations', card?.limitations),
    named(divV(facts.filter(([, x]) => x).map(([k, x]) => divH([span(k, 'pp-muted'), span(x)], 'pp-card-fact')),
      'pp-card-facts'), 'model-card-facts'),
    ...(others.length ? [span('Other versions', 'pp-heading'), named(divV(others.map((o) =>
      span(`${o.model_version} · ${o.approval.replace('_', ' ')}`, 'pp-muted'))), 'model-card-versions')] : []),
  ], 'pp-model-card'), 'model-card');
}

/** The card of [code] (a function's nqName), or a line saying the model is not catalogued. */
export async function modelCardWidget(code: string): Promise<DG.Widget> {
  installStyles();
  const v = await modelCardOf(code);
  return DG.Widget.fromRoot(v ? modelCardElement(v) :
    named(span(`${code} is not in the portal's model catalogue`, 'pp-muted'), 'model-card-missing'));
}

/** Opens the card of [code] in a dialog. */
export async function showModelCard(code: string): Promise<void> {
  try {
    const w = await modelCardWidget(code);
    const dlg = ui.dialog({title: 'Model card'}).add(w.root);
    dlg.root.dataset.u2Name = 'model-card-dialog';
    dlg.show({width: 520});
  } catch (e: any) {
    grok.shell.error(e?.message ?? String(e));
  }
}
