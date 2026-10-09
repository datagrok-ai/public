/** The model card (F-12, FR-50): the scientific question first, then what the model is for, its
 * inputs with their validity range, its outputs, applicability, performance, lineage and approval.
 * One renderer serves the program page's saved runs, the `modelCard` function, the Workspace on Home
 * and the platform's info panel. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {badge, divH, divV, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {CardParam, cardPortfolios, catalogueCodesOf} from '../model/catalogue';
import {ModelCardView, modelCardOf} from '../catalogue';
import {catalogueCodes} from '../save-run';
import {installStyles} from './styles';
import {inMyCatalogues, myPortfolioCodes} from '../secured-rows';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const slug = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

const approvalVariant: Record<string, BadgeVariant> = {
  draft: 'warning', in_review: 'accent', approved: 'success', rejected: 'error', retired: 'error',
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

/** The card of [code], when one of the user's portfolios catalogues the model; null otherwise. The
 * model tables are readable by all users, so this is the app's scoping, not the server's. */
async function visibleCardOf(code: string): Promise<ModelCardView | null> {
  const v = await modelCardOf(code);
  return v && await inMyCatalogues(v.portfolios) ? v : null;
}

/** The card of [code] (a function's nqName), or a line saying the model is not catalogued. */
export async function modelCardWidget(code: string): Promise<DG.Widget> {
  installStyles();
  const v = await visibleCardOf(code);
  return DG.Widget.fromRoot(v ? modelCardElement(v) :
    named(span(`${code} is not in your portfolios' model catalogues`, 'pp-muted'), 'model-card-missing'));
}

/** The card beside a model's preview in the Workspace tab on Home (PowerPack's `workspaceModelPanel` role):
 * the card of the catalogued model [func] is, by its nqName; null for a function the catalogue does not list. */
export async function workspaceModelCardWidget(func: DG.Func): Promise<DG.Widget | null> {
  const codes = [...new Set([func.nqName, func.package ? `${func.package.name}:${func.name}` : ''].filter((c) => c))];
  for (const code of codes) {
    const v = await visibleCardOf(code);
    if (v) {
      installStyles();
      return DG.Widget.fromRoot(named(divV([modelCardElement(v)]), 'workspace-model-card'));
    }
  }
  return null;
}

/** Which portfolios label a card in the info panel: one portfolio, several, or none (a model no
 * portfolio declares). */
export const cardLabel = (portfolios: string[]): string =>
  !portfolios.length ? '' : `In the catalogue of ${portfolios.join(', ')}`;

/** The "Model card" info panel of [func], wherever a function is the current object: the card of
 * each catalogue entry the function is, labelled with the user's portfolios that list it; null for
 * a function no catalogue of the user's portfolios lists. A function the catalogue does not know
 * costs no query past the session's list of catalogued codes. Never throws. */
export async function modelCardPanelWidget(func: DG.Func): Promise<DG.Widget | null> {
  try {
    const codes = catalogueCodesOf(func.nqName, func.package?.name, func.name, await catalogueCodes());
    if (!codes.length)
      return null;
    const mine = await myPortfolioCodes();
    const entries: HTMLElement[] = [];
    for (const v of await Promise.all(codes.map((c) => modelCardOf(c)))) {
      const shown = v ? cardPortfolios(v.portfolios, mine) : null;
      if (!v || !shown)
        continue;
      const label = cardLabel(shown);
      entries.push(named(divV([...(label ? [named(span(label, 'pp-muted'), 'model-card-panel-portfolios')] : []),
        modelCardElement(v)]), 'model-card-panel-entry'));
    }
    if (!entries.length)
      return null;
    installStyles();
    return DG.Widget.fromRoot(named(divV(entries), 'model-card-panel'));
  } catch (e) {
    console.warn('Model card panel:', e);
    return null;
  }
}

/** Adds the "Model card" pane to the context panel of every catalogued model function, wherever
 * it is the current object. The platform does not offer `meta.role: panel` functions for entities
 * such as functions (only for data: tables, columns, cells, files), so the pane is added when the
 * function's accordion is built, as Compute2 adds its REST pane. Nothing for other objects. */
export function addModelCardPanes(): void {
  grok.events.onAccordionConstructed.subscribe(async (acc: DG.Accordion) => {
    const func = acc?.context ? DG.toJs(acc.context) : null;
    if (!(func instanceof DG.Func) || acc.getPane(MODEL_CARD_PANE))
      return;
    const w = await modelCardPanelWidget(func);
    if (w && !acc.getPane(MODEL_CARD_PANE))
      acc.addPane(MODEL_CARD_PANE, () => w.root, true);
  });
}

export const MODEL_CARD_PANE = 'Model card';

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
