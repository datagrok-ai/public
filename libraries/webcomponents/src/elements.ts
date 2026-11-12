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
  private innerEditor?: DG.InputBase<boolean>;

  private get editor() {
    if (!this.innerEditor) {
      const editor = this.innerEditor = ui.input.toggle('Caption');
      this.appendChild(editor.root);
      editor.onChanged.subscribe(() => this.dispatchEvent(new CustomEvent('value-changed', {detail: editor.value})));
    }
    return this.innerEditor;
  }

  set caption(caption: string) {
    this.editor.caption = caption;
  }

  set nullable(nullable: boolean) {
    this.editor.nullable = nullable;
  }

  set value(value: boolean) {
    const editor = this.editor;
    editor.notify = false;
    editor.value = value;
    editor.notify = true;
  }

  set tooltip(tooltip: string) {
    this.editor.setTooltip(tooltip);
  }
}

export interface DGToggleInputT extends DGToggleInput {};

export class DGIconFA extends HTMLElement {
  _name = 'edit';
  _cursor = 'pointer';
  _animation = null as null | string;
  _tooltip = null as null | string;
  _faStyle = 'fal' as 'fal' | 'fas' | 'far' | 'fad';

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

export class DGIconImage extends HTMLElement {
  _name = 'edit';
  _cursor = 'pointer';
  _path = null as null | string;
  _tooltip = null as null | string;

  private render() {
    ui.empty(this);
    if (!this._path)
      return;
    const t = ui.iconImage(this._name, this._path, null, this._tooltip);
    t.style.cursor = this._cursor;
    // just let the container deal with the icon size
    t.style.width = '100%';
    t.style.height = '100%';
    this.appendChild(t);
  }

  set name(val: string) {
    this._name = val;
    this.render();
  }

  set path(val: string) {
    this._path = val;
    this.render();
  }

  set cursor(val: string) {
    this._cursor = val;
    this.render();
  }

  set tooltip(val: string | null) {
    this._tooltip = val;
    this.render();
  }
}

export interface DGIconImageT extends DGIconImage {};

export class DGComboPopup extends HTMLElement {
  private _caption: string | HTMLElement = 'Caption';
  private _items = [] as string[];
  private _tooltip: string | null = null;
  private _rightAligned: boolean = false;

  public set caption(val: string | HTMLElement) {
    this._caption = val;
    this.render();
  }

  public set items(val: string[]) {
    this._items = val;
    this.render();
  }

  public set tooltip(val: string | null) {
    this._tooltip = val;
    this.render();
  }

  public set rightAligned(val: boolean) {
    this._rightAligned = val;
    this.render();
  }

  private render() {
    ui.empty(this);
    const newPopup = ui.comboPopup(
      this._caption,
      this._items,
      (item) => this.dispatchEvent(new CustomEvent(
        'selected', {detail: {item, itemIdx: this._items.findIndex((i) => i === item)}}),
      ),
    );
    if (this._rightAligned)
      newPopup.classList.add('right-aligned');
    ui.tooltip.bind(newPopup, this._tooltip);
    newPopup.style.height = '24px';
    newPopup.style.minWidth = '0px';
    newPopup.onclick = (ev) => {
      ev.stopPropagation();
    };
    this.appendChild(newPopup);
  }
}

export interface DGComboPopupT extends DGComboPopup {};
