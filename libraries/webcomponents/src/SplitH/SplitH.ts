
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class DGSplitH extends HTMLElement {
  // DG API
  private _resize = true;

  private defaultSlot!: HTMLSlotElement;
  private inner = null as HTMLElement | null;

  private initialSlots = [] as HTMLElement[];
  private inited = false;

  constructor() {
    super();

    if (!this.shadowRoot)
      this.attachShadow({mode: 'open'});

    this.defaultSlot = this.shadowRoot!.appendChild(ui.element('slot') as HTMLSlotElement);
  }

  private render() {
    if (!this.inited) {
      this.initialSlots = this.defaultSlot.assignedElements() as HTMLElement[];
      this.inited = true;
    }

    const newInner = this.appendChild(ui.splitH(
      this.initialSlots,
      {style: {height: '100%', width: '100%'}},
      this._resize,
    ));

    Array.from(newInner.children).forEach((child) => (child as HTMLElement).style.maxWidth = 'fit-content');

    if (this.inner) {
      ui.empty(this.inner);
      this.inner.remove();
    };

    this.inner = newInner;
  }

  set resize(val: boolean) {
    this._resize = val;
    this.render();
  }
}

export interface DGSplitHT extends DGSplitH {};
