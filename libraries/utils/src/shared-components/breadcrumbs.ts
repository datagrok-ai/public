import * as ui from 'datagrok-api/ui';
import {Observable, fromEvent} from 'rxjs';
import {map} from 'rxjs/operators';

import '../../css/breadcrumbs.css';


export class Breadcrumbs {
  path: string[];
  root: HTMLDivElement;

  constructor(path: string[]) {
    this.root = ui.div();
    this.path = path;

    this.root = ui.divH(path.map((element) => ui.div(ui.span([element], `ui-breadcrumbs-text-element ${element}`),
      'ui-breadcrumbs-element')), 'ui-breadcrumbs');
  }

  get onPathClick(): Observable<string[]> {
    const pathElements = this.root.getElementsByClassName('ui-breadcrumbs-text-element');

    return fromEvent<MouseEvent>(pathElements, 'click')
      .pipe(
        map((event) => {
          const currentElement = (event.target as Element).innerHTML;
          const currentPath = this.path.slice(0, this.path.indexOf(currentElement) + 1);

          return currentPath;
        })
      );
  }
}
