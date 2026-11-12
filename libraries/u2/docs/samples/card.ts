import {Card, Control, badge, component, divH, iconButton, link, notify, signal, tag, timestamp}
  from '@datagrok-libraries/u2';

interface Compound {
  id: string;
  name: string;
  target: string;
  stage: 'Approved' | 'In review';
  tags: string[];
  updated: Date;
}

export function compoundCard(c: Compound): Card {
  const starred = signal(false);
  return new Card({
    icon: 'flask',
    title: c.name,
    subtitle: `${c.id} · ${c.target}`,
    actions: [iconButton('star', () => {}, {tooltip: 'Star', toggle: starred})],
    children: [divH([badge(c.stage, {variant: c.stage === 'Approved' ? 'success' : 'warning'}), ...c.tags.map((t) => tag(t))],
      {gap: 's'})],
    footer: [timestamp(c.updated), link('ChEMBL', `https://www.ebi.ac.uk/chembl/compound_report_card/${c.id}/`)],
    onClick: () => notify.info(`Opening ${c.name}`),
  });
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = component(() => compoundCard({id: 'CHEMBL25', name: 'Aspirin', target: 'PTGS1',
    stage: 'Approved', tags: ['NSAID', 'oral'], updated: new Date(2026, 9, 21, 14, 5)}));
  host.append(demo.root);
  return demo;
}
