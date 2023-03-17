import * as ui from 'datagrok-api/ui';
import {Observable, fromEvent} from 'rxjs';
import {filter} from 'rxjs/operators';

import '../../css/drop-down.css';


export class DropDown {
  private _element: HTMLElement;
  private _dropDownElement: HTMLDivElement;

  private _isMouseOverElement: boolean;
  isExpanded: boolean;

  private _label: string | Element;
  private _rootElement: HTMLElement;
  root: HTMLDivElement;

  // string | Element instead of iconName and label
  // instead of element make a creation function
  constructor(label: string | Element, elementCreation: () => HTMLElement) {
    this._element = ui.div();
    this._dropDownElement = ui.div();

    this._isMouseOverElement = false;
    this.isExpanded = false;

    this._label = label;
    this.root = ui.div();
    this._rootElement = ui.div();

    this._updateElement(elementCreation);
    this._initEventListeners();
  }

  private _updateElement(elementCreation: () => HTMLElement) {
    this._element = elementCreation();

    this._dropDownElement = ui.div(ui.div(this._element));
    this._dropDownElement.style.visibility = 'hidden';
    this._dropDownElement.className = 'ui-drop-down-menu';

    const dropDown = ui.div(this._dropDownElement);
    dropDown.className = 'ui-drop-down-menu-fixed';

    this._rootElement = ui.div([this._label]);
    this._rootElement.className = 'ui-drop-down-item';

    this.root = ui.div([this._rootElement, dropDown]);
    this.root.className = 'ui-drop-down-root';
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

  // take styles from native dropdown !!!!
  // write all events with rxjs (Observable)
  get onExpand(): Observable<MouseEvent> {
    return fromEvent<MouseEvent>(this.root, 'click')
      .pipe(
        filter(() => !this._isMouseOverElement)
      );
  }

  get onElementClick(): Observable<MouseEvent> {
    return fromEvent<MouseEvent>(this._dropDownElement, 'click');
  }
}
