import * as ui from 'datagrok-api/ui';

const COLORS = {
  pmx: '#185FA5',
  cp: '#1D9E75',
  bs: '#534AB7',
  pm: '#BA7517',
  border: '#888780',
};

interface PipelineStage {
  num: number;
  title: string;
  description: string;
  personas: ('pmx' | 'cp' | 'bs' | 'pm')[];
  roles: Record<string, 'executor' | 'consumer'>;
}

const stages: PipelineStage[] = [
  {num: 1, title: 'Input data', description: 'Dosing records, PK samples, PD biomarkers, Covariates (weight, age, GFR)',
    personas: ['pmx', 'bs'], roles: {pmx: 'executor', bs: 'consumer'}},
  {num: 2, title: 'Data assembly and formatting', description: 'NMTRAN format, BQL handling, quality checks',
    personas: ['pmx', 'bs'], roles: {pmx: 'executor', bs: 'consumer'}},
  {num: 3, title: 'Base structural model', description: 'PK structure + random effects (IIV, IOV, residual)',
    personas: ['pmx'], roles: {pmx: 'executor'}},
  {num: 4, title: 'Parameter estimation', description: 'FOCE-I / SAEM / Bayesian → θ, Ω, Σ',
    personas: ['pmx'], roles: {pmx: 'executor'}},
  {num: 5, title: 'Base model diagnostics', description: 'GOF plots, CWRES, ETA distributions, shrinkage',
    personas: ['pmx', 'cp'], roles: {pmx: 'executor', cp: 'consumer'}},
  {num: 6, title: 'Covariate model', description: 'ETA vs covariates → forward / backward selection',
    personas: ['pmx', 'bs'], roles: {pmx: 'executor', bs: 'consumer'}},
  {num: 7, title: 'Model qualification', description: 'Bootstrap, VPC, posterior predictive check',
    personas: ['pmx', 'cp'], roles: {pmx: 'executor', cp: 'consumer'}},
  {num: 8, title: 'Simulation', description: 'Monte Carlo → dosing regimens, virtual cohorts',
    personas: ['pmx', 'cp', 'pm'], roles: {pmx: 'executor', cp: 'consumer', pm: 'consumer'}},
];

interface Deliverable { title: string; subtitle: string; personas: ('pmx' | 'cp' | 'bs' | 'pm')[] }

const deliverables: Deliverable[] = [
  {title: 'Dose selection', subtitle: 'recommendations', personas: ['pmx', 'cp', 'pm']},
  {title: 'Drug label', subtitle: 'labeling text', personas: ['cp', 'pm']},
  {title: 'Trial design', subtitle: 'next phase', personas: ['bs', 'pm']},
  {title: 'Regulatory', subtitle: 'FDA/EMA report', personas: ['pmx', 'cp']},
];

const PERSONA_LABELS: Record<string, string> = {
  pmx: 'PMx', cp: 'CP', bs: 'BS', pm: 'PM',
};

const PERSONA_FULL: Record<string, string> = {
  pmx: 'Pharmacometrician', cp: 'Clinical Pharmacologist', bs: 'Biostatistician', pm: 'Program Leadership',
};

const PERSONA_GLOSSARY_KEY: Record<string, string> = {
  pmx: 'PMx', cp: 'CP', bs: 'BS', pm: 'PM',
};

let _glossary: Map<string, string> | undefined;

function personaTooltip(key: string, suffix?: string): (() => HTMLElement) | string {
  const gKey = PERSONA_GLOSSARY_KEY[key];
  const def = _glossary?.get(gKey);
  if (def && def.includes('\n')) {
    return () => {
      const parts = def.split('\n');
      const lines = parts.map((p, idx) => {
        const line = ui.element('div');
        line.textContent = idx === 0 && suffix ? `${p} — ${suffix}` : p;
        if (idx === 0)
          line.style.cssText = 'font-weight:bold;margin-bottom:4px';
        else
          line.style.cssText = 'font-style:italic;opacity:0.9';
        return line;
      });
      const el = ui.divV(lines);
      el.style.maxWidth = '360px';
      return el;
    };
  }
  return suffix ? `${PERSONA_FULL[key]} — ${suffix}` : PERSONA_FULL[key];
}

function dot(color: string, filled: boolean): HTMLElement {
  const d = ui.element('span');
  d.style.display = 'inline-block';
  d.style.width = '10px';
  d.style.height = '10px';
  d.style.borderRadius = '50%';
  d.style.marginRight = '3px';
  d.style.verticalAlign = 'middle';
  if (filled) {
    d.style.backgroundColor = color;
  } else {
    d.style.border = `2px solid ${color}`;
    d.style.boxSizing = 'border-box';
  }
  return d;
}

function personaDots(personas: ('pmx' | 'cp' | 'bs' | 'pm')[], roles: Record<string, string>): HTMLElement {
  const container = ui.element('span');
  container.style.display = 'inline-flex';
  container.style.gap = '2px';
  container.style.alignItems = 'center';
  for (const p of personas) {
    const isExecutor = roles[p] === 'executor';
    const d = dot(COLORS[p], isExecutor);
    ui.tooltip.bind(d, personaTooltip(p, isExecutor ? 'executor' : 'consumer'));
    container.appendChild(d);
  }
  return container;
}

function stageRow(stage: PipelineStage): HTMLElement {
  const numEl = ui.element('span');
  numEl.textContent = String(stage.num);
  numEl.style.cssText = 'display:inline-flex;align-items:center;justify-content:center;width:24px;height:24px;border-radius:50%;background:var(--grey-2,#e0e0e0);font-size:12px;font-weight:bold;flex-shrink:0';

  const titleEl = ui.element('span');
  titleEl.textContent = stage.title;
  titleEl.style.fontWeight = 'bold';
  titleEl.style.fontSize = '14px';

  const descEl = ui.element('span');
  descEl.textContent = stage.description;
  descEl.style.color = 'var(--grey-5,#666)';
  descEl.style.fontSize = '12px';

  const dots = personaDots(stage.personas, stage.roles);

  const textBlock = ui.divV([titleEl, descEl]);
  textBlock.style.gap = '2px';
  textBlock.style.flex = '1';
  textBlock.style.minWidth = '0';

  const row = ui.divH([numEl, textBlock, dots]);
  row.style.cssText = 'align-items:center;gap:12px;padding:8px 12px;border:1px solid var(--grey-2,#ddd);border-radius:6px;background:white';

  return row;
}

function arrow(): HTMLElement {
  const a = ui.element('div');
  a.style.cssText = 'text-align:center;color:var(--grey-3,#aaa);font-size:18px;line-height:1;margin:2px 0';
  a.textContent = '↓';
  return a;
}


function deliverableCard(d: Deliverable): HTMLElement {
  const titleEl = ui.element('div');
  titleEl.textContent = d.title;
  titleEl.style.cssText = 'font-weight:bold;font-size:13px;text-align:center';

  const subEl = ui.element('div');
  subEl.textContent = d.subtitle;
  subEl.style.cssText = 'font-size:11px;color:var(--grey-5,#666);text-align:center';

  const dotsRow = ui.element('div');
  dotsRow.style.cssText = 'display:flex;gap:3px;justify-content:center;margin-top:4px';
  for (const p of d.personas) {
    const dd = dot(COLORS[p], true);
    ui.tooltip.bind(dd, personaTooltip(p));
    dotsRow.appendChild(dd);
  }

  const card = ui.divV([titleEl, subEl, dotsRow]);
  card.style.cssText = 'padding:10px 8px;border:1px solid var(--grey-2,#ddd);border-radius:6px;background:white;min-width:120px;flex:1';

  return card;
}

export function createPipelineDiagram(glossary?: Map<string, string>): HTMLElement {
  _glossary = glossary;
  // Title
  const title = ui.h2('Population PK/PD modeling pipeline');
  const subtitle = ui.element('div');
  subtitle.textContent = 'with persona involvement by stage';
  subtitle.style.cssText = 'color:var(--grey-5,#666);margin-bottom:12px;font-size:13px';

  // Legend
  const legendItems: HTMLElement[] = [];
  for (const [key, label] of Object.entries(PERSONA_LABELS)) {
    const color = COLORS[key as keyof typeof COLORS];
    const filled = dot(color, true);
    const outline = dot(color, false);
    const nameEl = ui.element('span');
    nameEl.textContent = ` ${label}`;
    nameEl.style.cssText = 'font-size:12px;font-weight:bold;margin-right:12px';
    const item = ui.element('span');
    item.style.cssText = 'display:inline-flex;align-items:center';
    item.appendChild(filled);
    item.appendChild(outline);
    item.appendChild(nameEl);
    ui.tooltip.bind(item, personaTooltip(key, '● executor, ○ consumer'));
    legendItems.push(item);
  }

  const legendRow = ui.divH(legendItems);
  legendRow.style.cssText = 'gap:4px;flex-wrap:wrap;margin-bottom:4px';

  const legendNote = ui.element('div');
  legendNote.style.cssText = 'font-size:11px;color:var(--grey-4,#999);margin-bottom:16px';
  legendNote.textContent = '● executor  ○ consumer';

  // Stages
  const stageElements: HTMLElement[] = [];
  for (let s = 0; s < stages.length; s++) {
    stageElements.push(stageRow(stages[s]));
    if (s < stages.length - 1)
      stageElements.push(arrow());
  }

  const stagesContainer = ui.divV(stageElements);
  stagesContainer.style.gap = '0';

  // Arrow to deliverables
  const arrowToDel = arrow();
  arrowToDel.style.margin = '6px 0';

  // Deliverables
  const delTitle = ui.h3('Deliverables');
  delTitle.style.marginTop = '0';

  const delCards = ui.divH(deliverables.map(deliverableCard));
  delCards.style.cssText = 'gap:10px;flex-wrap:wrap';

  // Footer
  const footer = ui.element('div');
  footer.style.cssText = 'font-size:11px;color:var(--grey-5,#666);margin-top:16px;line-height:1.6';
  footer.innerHTML = 'PMx leads stages 2–8 · CP joins at diagnostics, simulation, and deliverables<br>' +
    'BS contributes to data, covariates, and trial design · PM consumes simulations and decisions';

  const container = ui.divV([title, subtitle, legendRow, legendNote, stagesContainer, arrowToDel, delTitle, delCards, footer]);
  container.style.cssText = 'max-width:680px;padding:16px;margin:0 auto';

  return container;
}
