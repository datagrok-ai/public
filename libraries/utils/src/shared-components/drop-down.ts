import * as ui from 'datagrok-api/ui';

export class DropDown {
  private _elements: string[];
  private _root?: HTMLDivElement;
  private _dropDownElements?: HTMLDivElement;
  private _isMouseOverElement: boolean;
  private _expanded: boolean;

  constructor(elements: string[]) {
    this._elements = [];
    this._isMouseOverElement = false;
    this._expanded = false;
    this._updateElements(elements);

    this._initEventListeners();
  }

  private _initEventListeners() {
    this._root!.addEventListener('mousedown', (e) => {
      // check if the button is LMB
      if (e.button !== 0) return;
      if (this._isMouseOverElement) return;

      if (!this._expanded)
        this._expand();
      else
        this._collapse();
    });

    this._dropDownElements!.addEventListener('mouseover', () => {
      this._isMouseOverElement = true;
    }, false);

    this._dropDownElements!.addEventListener('mouseleave', () => {
      this._isMouseOverElement = false;
    }, false);
  }

  private _updateElements(elements: string[]) {
    this._elements = elements;

    this._dropDownElements = ui.div(this._elements.map((elem) => ui.div(elem)), 'd4-combo-drop-down');
    this._dropDownElements.style.visibility = 'hidden';

    const dropDown = ui.div(this._dropDownElements, 'd4-combo-drop-down-fixed');
    this._root = ui.div(['DropDown', dropDown], 'd4-combo-popup');
  }

  private _expand() {
    this._expanded = true;
    this._root!.classList.add('d4-combo-popup-expanded');
    this._dropDownElements!.style.visibility = 'visible';
  }

  private _collapse() {
    this._expanded = false;
    this._root!.classList.remove('d4-combo-popup-expanded');
    this._dropDownElements!.style.visibility = 'hidden';
  }


  get root() {
    return this._root!;
  }
}
