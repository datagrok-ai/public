import * as ui from 'datagrok-api/ui';

import '../../css/drop-down.css';

export class DropDown {
  // private _elements: HTMLElement[];
  private _element: HTMLElement;
  private _dropDownElement: HTMLDivElement;

  private _isMouseOverElement: boolean;
  isExpanded: boolean;

  root: HTMLDivElement;
  private _rootElement: HTMLElement;
  private _icon: HTMLElement;
  private _label: HTMLElement;

  constructor(iconName: string = '', label: string = '', element: HTMLElement) {
    this._element = document.createElement('div');
    this._dropDownElement = document.createElement('div');

    this._isMouseOverElement = false;
    this.isExpanded = false;

    this._rootElement = document.createElement('div');
    this.root = document.createElement('div');

    this._icon = iconName ? ui.iconFA(iconName) : document.createElement('div');
    this._label = label ? ui.divText(label, 'ui-drop-down-item-label') : document.createElement('div');

    this._updateElements(element);
    this._initEventListeners();
  }


  private _initEventListeners() {
    this.root.addEventListener('mousedown', (e) => {
      // check if the button is LMB
      if (e.button !== 0) return;
      if (this._isMouseOverElement) return;

      this._setExpandedState(this.isExpanded);
    });

    this._dropDownElement.addEventListener('mouseover', () => {
      this._isMouseOverElement = true;
    }, false);

    this._dropDownElement.addEventListener('mouseleave', () => {
      this._isMouseOverElement = false;
    }, false);
  }

  private _setExpandedState(isExpanded: boolean) {
    this.isExpanded = !isExpanded;
    if (isExpanded) {
      this.root.classList.remove('ui-drop-down-expanded');
      this._dropDownElement.style.visibility = 'hidden';
      return;
    }

    this.root.classList.add('ui-drop-down-expanded');
    this._dropDownElement.style.visibility = 'visible';
  }

  private _updateElements(element: HTMLElement) {
    this._element = element;

    // this._dropDownElements = ui.div(this._elements.map((item)=> ui.div(item, 'ui-drop-down-menu-list-item')));
    this._dropDownElement = ui.div(ui.div(this._element, 'ui-drop-down-menu-list-item'));
    this._dropDownElement.style.visibility = 'hidden';
    this._dropDownElement.className = 'ui-drop-down-menu';

    const dropDown = ui.div(this._dropDownElement);
    dropDown.className = 'ui-drop-down-menu-fixed';

    this._rootElement = ui.div([this._icon, this._label]);
    this._rootElement.className = 'ui-drop-down-item';

    this.root = ui.div([this._rootElement, dropDown]);
    this.root.className = 'ui-drop-down-root';
  }


  onExpand(callback: Function) {
    this.root.addEventListener('click', () => {
      if (!this._isMouseOverElement)
        callback();
    }, false);
  }

  onElementClick(callback: Function) {
    this._dropDownElement.addEventListener('click', (event) => {
      event.stopImmediatePropagation();
      callback();
    }, false);
  }
}
