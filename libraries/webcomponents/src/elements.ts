import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class DGButton extends HTMLButtonElement {
  constructor() {
    super();
    this.classList.add('ui-btn', 'ui-btn-ok');
  }
}

export interface DGButtonT extends DGButton {};

export class DGBigButton extends HTMLButtonElement {
  constructor() {
    super();
    this.classList.add('ui-btn', 'ui-btn-ok', 'ui-btn-raised');
  }
}

export interface DGBigButtonT extends DGBigButton {};

export class DGIconFA extends HTMLElement {
  _name = 'edit';
  _cursor = 'pointer';
  _animation = null as null | string;
  _tooltip = null as null | string;

  constructor() {
    super();
  }

  private render() {
    ui.empty(this);
    const t = ui.iconFA(this._name, null, this._tooltip);
    t.style.cursor = this._cursor;
    if (this._animation) t.classList.add(`fa-${this._animation}`);
    this.appendChild(t);
  }

  set name(val: string) {
    this._name = val;
    this.render();
  }

  set cursor(val: string) {
    this._cursor = val;
    this.render();
  }

  set animation(val: string | null) {
    this._animation = val;
    this.render();
  }

  set tooltip(val: string | null) {
    this._tooltip = val;
    this.render();
  }
}

export interface DGIconFAT extends DGIconFA {};
