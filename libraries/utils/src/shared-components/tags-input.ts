import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, Subject} from 'rxjs';

import '../../css/tags-input.css';


export class TagsInput extends DG.InputBase {
  tags: string[] = [];
  _tagsDiv: HTMLDivElement = ui.div();

  _onTagAdded: Subject<string> = new Subject<string>();
  _onTagRemoved: Subject<string> = new Subject<string>();


  constructor(tags: string[]) {
    const inputElement = ui.stringInput('', '');
    super(inputElement.dart);

    this._init(tags);
  }

  _init(tags: string[]) {
    const addTagIcon = ui.iconFA('plus', () => this.addTag((this.input as HTMLInputElement).value));
    this.addOptions(addTagIcon);

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

  _isProper(tag: string): boolean {
    return !(tag === '' || this.tags.includes(tag));
  }


  addTag(tag: string) {
    if (!this._isProper(tag))
      return;

    this.tags[this.tags.length] = tag;
    const currentTag = this._createTag(tag);
    this._tagsDiv.append(currentTag);
    (this.input as HTMLInputElement).value = '';
    this._onTagAdded.next(tag);
  }

  removeTag(tag: string) {
    this.tags.splice(this.tags.indexOf(tag), 1);
    const currentTag = this._tagsDiv.querySelector(`[data-tag="${tag}"]`);
    this._tagsDiv.removeChild(currentTag!);
    this._onTagRemoved.next(tag);
  }


  get onTagAdded(): Observable<string> {
    return this._onTagAdded;
  }

  get onTagRemoved(): Observable<string> {
    return this._onTagRemoved;
  }
}
