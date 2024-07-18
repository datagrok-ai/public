import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';

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

export class DGSplitH extends HTMLElement {
  // DG API
  private resize = false;

  private defaultSlot!: HTMLSlotElement;
  private inner = null as HTMLElement | null;

  constructor() {
    super();

    if (!this.shadowRoot)
      this.attachShadow({mode: 'open'});

    this.defaultSlot = this.shadowRoot!.appendChild(ui.element('slot') as HTMLSlotElement);
  }

  static get observedAttributes() {
    return ['resize'];
  }

  connectedCallback() {
    if (!this.inner) {
      this.inner = this.appendChild(ui.splitH(
        Array.from(this.defaultSlot.assignedElements() as HTMLElement[]),
        {},
        true,
      ));
    }
  }
}

export interface DGSplitHT extends DGSplitH {};
