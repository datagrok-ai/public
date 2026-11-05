import {ReadonlySignal, Section, divH, span} from '@datagrok-libraries/u2';

export const CITIES = [
  'Amsterdam', 'Athens', 'Barcelona', 'Berlin', 'Boston', 'Brussels', 'Budapest', 'Chicago',
  'Copenhagen', 'Dublin', 'Geneva', 'Helsinki', 'Kyiv', 'Lisbon', 'London', 'Madrid', 'Milan',
  'Munich', 'New York', 'Oslo', 'Paris', 'Prague', 'Rome', 'Seattle', 'Stockholm', 'Tokyo',
  'Toronto', 'Vienna', 'Warsaw', 'Zurich',
];

export const DRUGS = [
  {name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O'},
  {name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'},
  {name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'},
  {name: 'Paracetamol', smiles: 'CC(=O)NC1=CC=C(C=C1)O'},
  {name: 'Naproxen', smiles: 'CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O'},
];

/** A "name = value" status line: `[role=status]` with a `<label>` and an `<output>`, the generic
 * readout the bdd library reads (`value of <name> readout`). */
export function readout(name: string, source: string | ReadonlySignal<unknown>): HTMLElement {
  const label = document.createElement('label');
  label.textContent = `${name} = `;
  const out = document.createElement('output');
  out.append(span(source));
  const root = divH([label, out], 'u2demo-status');
  root.setAttribute('role', 'status');
  return root;
}

/** An A/B grid (`.u2demo-ab`) filled row by row. Each side is a section named as the page names
 * it, the API that builds it beside the title, spanning every row and taking the grid's rows as
 * its own (css/u2demo.css): one container per side, whose cells still share a row with the other
 * side's. Every item carries its row as `--ab-row` and the CSS places it, so a narrow pane can
 * stack each u2 cell under its platform twin. */
export interface AbGrid {
  readonly root: HTMLElement;
  readonly platform: Section;
  readonly u2: Section;
  /** The first column's cell, the platform's and u2's, on the next row. */
  row(first: HTMLElement, platform: HTMLElement, u2: HTMLElement): void;
  /** A full-width row: a note or a category. */
  full(el: HTMLElement): void;
}

/** Fills `root` with an {@link AbGrid}; built inside the page's scope, which disposes the sections. */
export function abGrid(root: HTMLElement, platformApi: string, u2Api: string): AbGrid {
  let rows = 1;
  const spanRows = () => root.style.setProperty('--ab-rows', String(rows));
  const corner = span('');
  at(corner, 1);
  const platform = side('Platform', platformApi, 'u2demo-ab-platform');
  const u2 = side('u2', u2Api, 'u2demo-ab-u2');
  root.append(platform.root, u2.root, corner);
  spanRows();
  return {root, platform, u2,
    row: (first, platformCell, u2Cell) => {
      const row = ++rows;
      for (const cell of [first, platformCell, u2Cell])
        at(cell, row);
      root.append(first);
      platform.body.append(platformCell);
      u2.body.append(u2Cell);
      spanRows();
    },
    full: (el) => {
      el.classList.add('u2demo-ab-full');
      at(el, ++rows);
      root.append(el);
      spanRows();
    },
  };
}

function at(el: HTMLElement, row: number): void {
  el.style.setProperty('--ab-row', String(row));
}

function side(title: string, api: string, className: string): Section {
  const section = new Section({title, collapsible: false});
  section.root.classList.add('u2demo-ab-side', className);
  section.header.classList.add('u2demo-ab-head');
  section.body.classList.add('u2demo-ab-cells');
  section.header.append(span(api, 'u2demo-code'));
  return section;
}
