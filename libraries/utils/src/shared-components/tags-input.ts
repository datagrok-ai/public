/* eslint-disable indent */
/* eslint-disable no-tabs */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, fromEvent} from 'rxjs';
import {filter, map} from 'rxjs/operators';

import '../../css/tags-input.css';


export class TagsInput extends DG.InputBase {
  tags: string[];
  _tagsDiv: HTMLDivElement;

  constructor(tags: string[]) {
    const inputElement = ui.stringInput('', '');
    super(inputElement.dart);

    this.tags = tags;
    this._tagsDiv = ui.div(tags.map((tag) => ui.span([ui.span([tag])], `ui-tag ${tag}`)), 'ui-tag-list');
    this.root.append(this._tagsDiv);

    this._initEventListeners();
  }

	_initEventListeners() {
		this.input.addEventListener('keyup', (event: KeyboardEvent) => {
      if (event.code === 'Enter')
        this.addTag((this.input as HTMLInputElement).value);
    });
	}


  addTag(tag: string) {
    if (!this.tags.includes(tag)) {
      this.tags[this.tags.length] = tag;
      this._tagsDiv.append(ui.span([ui.span([tag])], 'ui-tag'));
      (this.input as HTMLInputElement).value = '';
    }
  }

  removeTag(tag: string) {
    this.tags.splice(this.tags.indexOf(tag), 1);
		const currentTag = this._tagsDiv.getElementsByClassName(tag)[0];
		this._tagsDiv.removeChild(currentTag);
  }

  get onTagAdded(): Observable<KeyboardEvent> {
    return fromEvent<KeyboardEvent>(this.input, 'keyup').pipe(filter((event: KeyboardEvent) => event.code === 'Enter'));
  }
}
