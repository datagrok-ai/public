import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, Subject} from 'rxjs';

import '../../css/tags-input.css';


export class TagsInput extends DG.InputBase {
  private _tags: string[] = [];
  private _tagsDiv: HTMLDivElement = ui.div();

  private _onTagAdded: Subject<string> = new Subject<string>();
  private _onTagRemoved: Subject<string> = new Subject<string>();


  constructor(name: string, tags: string[], showBtn: boolean) {
    const inputElement = ui.stringInput(name, '');
    super(inputElement.dart);

    this._init(tags, showBtn);
  }

  private _init(tags: string[], showBtn: boolean) {
    if (showBtn) {
      const addTagIcon = ui.iconFA('plus', () => this.addTag((this.input as HTMLInputElement).value));
      this.addOptions(addTagIcon);
    }

    this._tags = tags;
    this._tagsDiv = ui.div(tags.map((tag) => { return this._createTag(tag); }), 'ui-tag-list');
    this.root.append(this._tagsDiv);

    this._initEventListeners();
  }

  private _initEventListeners() {
    this.input.addEventListener('keyup', (event: KeyboardEvent) => {
      if (event.code === 'Enter')
        this.addTag((this.input as HTMLInputElement).value);
    });
  }

  private _createTag(tag: string): HTMLElement {
    const icon = ui.iconFA('times', () => this.removeTag(tag));
    const currentTag = ui.span([ui.span([tag]), icon], `ui-tag`);
    currentTag.dataset.tag = tag;

    return currentTag;
  }

  private _isProper(tag: string): boolean {
    return !(tag === '' || this._tags.includes(tag));
  }


  addTag(tag: string) {
    if (!this._isProper(tag))
      return;

    this._tags[this._tags.length] = tag;
    const currentTag = this._createTag(tag);
    this._tagsDiv.append(currentTag);
    (this.input as HTMLInputElement).value = '';
    this._onTagAdded.next(tag);
  }

  removeTag(tag: string) {
    this._tags.splice(this._tags.indexOf(tag), 1);
    const currentTag = this._tagsDiv.querySelector(`[data-tag="${tag}"]`);
    this._tagsDiv.removeChild(currentTag!);
    this._onTagRemoved.next(tag);
  }

  getTags(): string[] {
    return this._tags;
  }

  setTags(tags: string[]) {
    this._tags = tags;
    this._tagsDiv = ui.div(tags.map((tag) => { return this._createTag(tag); }), 'ui-tag-list');
    const currentTagsDiv = this.root.getElementsByClassName('ui-tag-list')[0];
    currentTagsDiv.replaceWith(this._tagsDiv);
  }


  get onTagAdded(): Observable<string> {
    return this._onTagAdded;
  }

  get onTagRemoved(): Observable<string> {
    return this._onTagRemoved;
  }
}
