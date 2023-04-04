import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, fromEvent} from 'rxjs';
import {filter} from 'rxjs/operators';

import '../../css/tags-input.css';


export class TagsInput extends DG.InputBase {
  tags: string[];
  _tagsDiv: HTMLDivElement;

  constructor(tags: string[]) {
    const inputElement = ui.stringInput('', '');
    super(inputElement.dart);

    const addTagBtn = ui.iconFA('plus', () => this.addTag((this.input as HTMLInputElement).value));
    this.addOptions(addTagBtn);

    this.tags = tags;
    this._tagsDiv = ui.div(tags.map((tag) => { return this._createTag(tag); }), 'ui-tag-list');
    this.root.append(this._tagsDiv);

    this._initEventListeners();
  }

  _initEventListeners() {
    this.input.addEventListener('keyup', (event: KeyboardEvent) => {
      if (event.code === 'Enter')
        this.addTag((this.input as HTMLInputElement).value);
    });
  }

  _createTag(tag: string): HTMLElement {
    const icon = ui.iconFA('times', () => this.removeTag(tag));
    const currentTag = ui.span([ui.span([tag]), icon], `ui-tag`);
    currentTag.dataset.tag = tag;

    return currentTag;
  }


  addTag(tag: string) {
    if (tag === '' || this.tags.includes(tag))
      return;

    this.tags[this.tags.length] = tag;
    const currentTag = this._createTag(tag);
    this._tagsDiv.append(currentTag);
    (this.input as HTMLInputElement).value = '';
  }

  removeTag(tag: string) {
    this.tags.splice(this.tags.indexOf(tag), 1);
    const currentTag = this._tagsDiv.querySelector(`[data-tag="${tag}"]`);
    this._tagsDiv.removeChild(currentTag!);
  }

  // TODO: immplement onclick tag add and fix it
  get onTagAdded(): Observable<KeyboardEvent> {
    return fromEvent<KeyboardEvent>(this.input, 'keyup').pipe(filter((event: KeyboardEvent) => event.code === 'Enter'));
  }

  // TODO: fix it
  get onTagRemoved(): Observable<MouseEvent> {
    const removeIcons = this._tagsDiv.getElementsByClassName('grok-icon fal fa-times');
    return fromEvent<MouseEvent>(removeIcons, 'click');
  }
}
