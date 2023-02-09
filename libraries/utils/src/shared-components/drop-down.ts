import * as ui from 'datagrok-api/ui';
import '../../css/drop-down.css';

export class DropDown {
  private _rootElement: HTMLElement;
  private _elements: HTMLElement[];
  private _dropDownElements: HTMLDivElement;
  private _isMouseOverElement: boolean;
  private _expanded: boolean;

  root: HTMLDivElement;

  constructor(rootElement: HTMLElement, elements: HTMLElement[]) {
    this._rootElement = rootElement;
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
      this.root.classList.remove('ui-drop-down-expanded');
      this._dropDownElements.style.visibility = 'hidden';
      return;
    }

    this.root.classList.add('ui-drop-down-expanded');
    this._dropDownElements.style.visibility = 'visible';
  }

  private _updateElements(elements: HTMLElement[]) {
    this._elements = elements;

    this._dropDownElements = ui.div(this._elements.map((item) => ui.div(item, 'ui-drop-down-menu-list-item')));
    this._dropDownElements.style.visibility = 'hidden';
    this._dropDownElements.className = 'ui-drop-down-menu';

    const dropDown = ui.div(this._dropDownElements);
    dropDown.className = 'ui-drop-down-menu-fixed';

    this._rootElement.classList.add('ui-drop-down-item');

    this.root = ui.div([this._rootElement, dropDown]);
    this.root.className = 'ui-drop-down-root';
  }
}
