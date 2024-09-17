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

export class DGToggleInput extends HTMLElement {
  private _caption: string = 'Caption';

  private innerEditor = null as null | DG.InputBase<boolean>;

  constructor() {
    super();
  }

  private render() {
    ui.empty(this);
    this.innerEditor = ui.input.toggle(this._caption);
    this.appendChild(this.innerEditor.root);

    this.innerEditor.onChanged.subscribe(() =>
      this.dispatchEvent(new CustomEvent('value-changed', {detail: this.innerEditor!.value})),
    );
  }

  set caption(caption: string) {
    if (!this.innerEditor) this.render();

    this.innerEditor!.caption = caption;
  }

  set nullable(nullable: boolean) {
    if (!this.innerEditor) this.render();

    this.innerEditor!.nullable = nullable;
  }

  set value(value: boolean) {
    if (!this.innerEditor) this.render();

    this.innerEditor!.value = value;
  }

  set tooltip(tooltip: string) {
    if (!this.innerEditor) this.render();

    this.innerEditor!.setTooltip(tooltip);
  }
}

export interface DGToggleInputT extends DGToggleInput {};

export class DGIconFA extends HTMLElement {
  _name = 'edit';
  _cursor = 'pointer';
  _animation = null as null | string;
  _tooltip = null as null | string;
  _faStyle = 'fal' as 'fal' | 'fas' | 'far' | 'fad';

  constructor() {
    super();
  }

  private render() {
    ui.empty(this);
    const t = ui.iconFA(this._name, null, this._tooltip);
    t.style.cursor = this._cursor;
    if (this._animation) t.classList.add(`fa-${this._animation}`);
    t.classList.remove('fal', 'fas', 'far', 'fad');
    t.classList.add(this._faStyle);
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

  set faStyle(val: typeof this._faStyle) {
    this._faStyle = val;
    this.render();
  }
}

export interface DGIconFAT extends DGIconFA {};

// unwrap away from an element; super basic but makes it consistent across our apps
function unwrap(el: HTMLElement) {
  if (el && el.parentNode) {
    // move all children out of the element
    while (el.firstChild)
      el.parentNode.insertBefore(el.firstChild, el);

    // remove the empty element
    el.remove();
  }
}

export class DGComboPopup extends HTMLDivElement {
  private _caption: string | HTMLElement = 'Caption';
  private _items = [] as string[];
  private _isExpanded: boolean = false;

  constructor() {
    super();

    this.addEventListener('click', () => {
      this._isExpanded = !this._isExpanded;
      this.render();
    });
  }

  public set caption(val: string | HTMLElement) {
    this._caption = val;
    this.render();
  }

  public set items(val: string[]) {
    this._items = val;
    this.render();
  }

  private render() {
    ui.empty(this);
    this.classList.add('d4-combo-popup');
    this.classList.toggle('d4-combo-popup-expanded', this._isExpanded);

    const list = this._items.map((item, itemIdx) => {
      const itemLabel = ui.label(item, 'ui-label');
      itemLabel.style.marginRight = '6px';

      const listItem = ui.div(itemLabel, 'd4-list-item');

      const listItemWrapper = ui.div(listItem);
      listItemWrapper.setAttribute('name', `${item}-host`);

      listItemWrapper.addEventListener('click', ()=> {
        this.dispatchEvent(new CustomEvent(
          'selected', {detail: {item, itemIdx}}),
        );
      });

      return listItemWrapper;
    });
    const listWrapper = ui.div(list, 'd4-list show-hover');
    listWrapper.setAttribute('data-widget', 'true');
    const hideableDropdown = ui.div(listWrapper, 'd4-combo-drop-down');
    hideableDropdown.style.visibility = this._isExpanded ? 'visible' : 'hidden';
    const dropdown = ui.div(hideableDropdown, 'd4-combo-drop-down-fixed');
    this.append(dropdown);
    this.append(this._caption);
  }
}

export interface DGComboPopupT extends DGComboPopup {};
