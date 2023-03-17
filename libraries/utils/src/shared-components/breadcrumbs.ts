import * as ui from 'datagrok-api/ui';
import {Observable, fromEvent} from 'rxjs';

import '../../css/breadcrumbs.css';


export class Breadcrumbs {
  path: string[];
  root: HTMLDivElement;

  constructor(path: string[]) {
    this.root = ui.div();
    this.path = path;

    this.root = ui.div(path.map((elem) => ui.span([elem], 'list-elem')), 'list');
  }

  get onPathClick(): Observable<MouseEvent> {
    return fromEvent<MouseEvent>(this.root, 'click');
  }
}
