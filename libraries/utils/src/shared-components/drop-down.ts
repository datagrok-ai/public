import * as ui from 'datagrok-api/ui';

export class DropDown<Element extends Node> {
  private _elements: Element[];
  private _dropDownElements: HTMLDivElement;
  private _isMouseOverElement: boolean;
  private _expanded: boolean;

  root: HTMLDivElement;


  constructor(elements: Element[]) {
    this._elements = [];
    this._dropDownElements = document.createElement('div');
    this._isMouseOverElement = false;
    this._expanded = false;
    this.root = document.createElement('div');
    this._updateElements(elements);

    this._initEventListeners();
  }


  private _initEventListeners() {
    this.root.addEventListener('mousedown', (e) => {
      // check if the button is LMB
      if (e.button !== 0) return;
      if (this._isMouseOverElement) return;

      this._setExpandedState(this._expanded);
    });

    this._dropDownElements.addEventListener('mouseover', () => {
      this._isMouseOverElement = true;
    }, false);

    this._dropDownElements.addEventListener('mouseleave', () => {
      this._isMouseOverElement = false;
    }, false);
  }

  private _setExpandedState(expanded: boolean) {
    this._expanded = !expanded;
    if (expanded) {
      this.root.classList.remove('d4-combo-popup-expanded');
      this._dropDownElements.style.visibility = 'hidden';
      return;
    }

    this.root.classList.add('d4-combo-popup-expanded');
    this._dropDownElements.style.visibility = 'visible';
  }

  private _updateElements(elements: Element[]) {
    this._elements = elements;

    this._dropDownElements = ui.div(this._elements, 'd4-combo-drop-down');
    this._dropDownElements.style.visibility = 'hidden';

    const dropDown = ui.div(this._dropDownElements, 'd4-combo-drop-down-fixed');
    this.root = ui.div(['DropDown', dropDown], 'd4-combo-popup');
  }
}
