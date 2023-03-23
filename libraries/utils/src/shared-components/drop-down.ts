import * as ui from 'datagrok-api/ui';
import {Observable, fromEvent} from 'rxjs';
import {filter} from 'rxjs/operators';

import '../../css/drop-down.css';


export class DropDown {
  private _element: HTMLElement;
  private _dropDownElement: HTMLDivElement;
  private _isMouseOverElement: boolean;
  private _label: string | Element;

  isExpanded: boolean;
  root: HTMLDivElement;


  constructor(label: string | Element, elementCreation: () => HTMLElement) {
    this._element = ui.div();
    this._dropDownElement = ui.div();
    this._isMouseOverElement = false;
    this._label = label;

    this.isExpanded = false;
    this.root = ui.div();

    this._updateElement(elementCreation);
    this._initEventListeners();
  }


  private _updateElement(elementCreation: () => HTMLElement) {
    this._element = elementCreation();

    this._dropDownElement = ui.div(ui.div(this._element), 'ui-drop-down-content');
    this._dropDownElement.style.visibility = 'hidden';

    this.root = ui.div([this._label, this._dropDownElement], 'ui-drop-down-root');
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

    document.addEventListener('click', (event) => {
      if (this.root.contains(event.target as Node)) return;
      if (!this.isExpanded) return;

      this._setExpandedState(this.isExpanded);
    });
  }

  private _setExpandedState(isExpanded: boolean) {
    this.isExpanded = !isExpanded;
    if (isExpanded) {
      this.root.classList.remove('ui-drop-down-root-expanded');
      this._dropDownElement.style.visibility = 'hidden';
      return;
    }

    this.root.classList.add('ui-drop-down-root-expanded');
    this._dropDownElement.style.visibility = 'visible';
  }


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
