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
