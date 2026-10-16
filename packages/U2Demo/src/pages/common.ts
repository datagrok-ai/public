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
 * stack each u2 cell under its platform twin. Built inside the page's scope, which disposes the
 * sections. */
export class AbGrid {
  readonly platform: Section;
  readonly u2: Section;
  private _row = 1;

  constructor(readonly root: HTMLElement, platformApi: string, u2Api: string) {
    const corner = span('');
    AbGrid._at(corner, 1);
    this.platform = AbGrid._side('Platform', platformApi, 'u2demo-ab-platform');
    this.u2 = AbGrid._side('u2', u2Api, 'u2demo-ab-u2');
    root.append(this.platform.root, this.u2.root, corner);
    this._span();
  }

  /** The first column's cell, the platform's and u2's, on the next row. */
  row(first: HTMLElement, platform: HTMLElement, u2: HTMLElement): void {
    const row = ++this._row;
    for (const cell of [first, platform, u2])
      AbGrid._at(cell, row);
    this.root.append(first);
    this.platform.body.append(platform);
    this.u2.body.append(u2);
    this._span();
  }

  /** A full-width row: a note or a category. */
  full(el: HTMLElement): void {
    el.classList.add('u2demo-ab-full');
    AbGrid._at(el, ++this._row);
    this.root.append(el);
    this._span();
  }

  private _span(): void {
    this.root.style.setProperty('--ab-rows', String(this._row));
  }

  private static _at(el: HTMLElement, row: number): void {
    el.style.setProperty('--ab-row', String(row));
  }

  private static _side(title: string, api: string, className: string): Section {
    const side = new Section({title, collapsible: false});
    side.root.classList.add('u2demo-ab-side', className);
    side.header.append(span(api, 'u2demo-code'));
    return side;
  }
}
