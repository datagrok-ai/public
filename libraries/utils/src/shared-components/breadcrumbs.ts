import * as ui from 'datagrok-api/ui';
import {Observable, fromEvent} from 'rxjs';

import '../../css/breadcrumbs.css';


export class Breadcrumbs {
  path: string[];
  root: HTMLDivElement;

  constructor(path: string[]) {
    this.root = ui.div();
    this.path = path;

    this.root = ui.divH(path.map((element) =>
      ui.div(ui.span([element], `${element}`), `breadcrumbs-element`)), 'breadcrumbs');
  }

  onPathClick(pathElementName: string): Observable<MouseEvent> {
    const pathElement = this.root.getElementsByClassName(`${pathElementName}`);
    return fromEvent<MouseEvent>(pathElement, 'click');
  }
}
