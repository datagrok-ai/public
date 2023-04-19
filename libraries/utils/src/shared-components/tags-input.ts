import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, Subject} from 'rxjs';

import '../../css/tags-input.css';


export class TagsInput extends DG.InputBase {
  private _tags: string[];
  private _tagsDiv: HTMLDivElement;
  private _addTagIcon: HTMLElement;

  private _onTagAdded: Subject<string> = new Subject<string>();
  private _onTagRemoved: Subject<string> = new Subject<string>();


  constructor(name: string, tags: string[], showBtn: boolean) {
    const inputElement = ui.stringInput(name, '');
    super(inputElement.dart);

    this._addTagIcon = showBtn ?
      ui.iconFA('plus', () => this.addTag((this.input as HTMLInputElement).value)) :
      ui.iconFA('');

    this._tags = tags;
    this._tagsDiv = ui.div(tags.map((tag) => { return this._createTag(tag); }), 'ui-tag-list');

    this._createRoot();
    this._initEventListeners();
  }


  private _createTag(tag: string): HTMLElement {
    const icon = ui.iconFA('times', () => this.removeTag(currentTag.innerText));
    const currentTag = ui.span([ui.span([tag]), icon], `ui-tag`);
    currentTag.dataset.tag = tag;

    currentTag.ondblclick = () => {
      const input = this._createTagEditInput(currentTag);
      (currentTag.firstElementChild as HTMLElement).innerText = '';
      currentTag.insertBefore(input, currentTag.firstElementChild);

      input.select();
      input.focus();
    };

    return currentTag;
  }

  private _createTagEditInput(currentTag: HTMLElement): HTMLInputElement {
    const tagValue = currentTag.innerText;
    const input = document.createElement('input');
    input.value = tagValue;

    let blurFlag = true;

    input.onblur = () => {
      if (!blurFlag)
        return;
      const newTag = input.value;
      input.remove();
      this._renameTag(tagValue, newTag, currentTag);
    };

    input.onkeyup = (event) => {
      if (event.key === 'Escape') {
        blurFlag = false;
        input.remove();
        (currentTag.firstElementChild as HTMLElement).innerText = tagValue;
      } else if (event.key === 'Enter') {
        blurFlag = false;
        const newTag = input.value;
        input.remove();
        this._renameTag(tagValue, newTag, currentTag);
      }
    };

    return input;
  }

  private _createRoot(): void {
    const inputContainer = ui.div([this.captionLabel, this.input, ui.div(this._addTagIcon, 'ui-input-options')],
      'ui-input-root');
    const tagContainer = ui.div([ui.label(' ', 'ui-input-label'), this._tagsDiv], 'ui-input-root');

    this.root.append(inputContainer, tagContainer);
    this.root.classList.add('ui-input-tags');
  }

  private _initEventListeners(): void {
    this.input.addEventListener('keyup', (event: KeyboardEvent) => {
      if (event.code === 'Enter')
        this.addTag((this.input as HTMLInputElement).value);
    });

    this.input.addEventListener('keydown', (event: KeyboardEvent) => {
      if ((this.input as HTMLInputElement).value.length === 0 && event.code === 'Backspace')
        this.removeTag(this._tags[0]);
    });
  }

  private _renameTag(oldTag: string, newTag: string, currentTag: HTMLElement): void {
    const currentTagSpan = currentTag.firstElementChild as HTMLElement;
    if (!this._isProper(newTag)) {
      currentTagSpan.innerText = oldTag;
      return;
    }

    currentTagSpan.innerText = newTag;
    currentTag.dataset.tag = newTag;
    this._tags[this._tags.indexOf(oldTag)] = newTag;
  }

  private _isProper(tag: string): boolean {
    return !(tag === '' || this._tags.includes(tag));
  }


  addTag(tag: string): void {
    if (!this._isProper(tag))
      return;

    this._tags.unshift(tag);
    const currentTag = this._createTag(tag);
    this._tagsDiv.insertBefore(currentTag, this._tagsDiv.firstChild);
    (this.input as HTMLInputElement).value = '';
    this._onTagAdded.next(tag);
  }

  removeTag(tag: string): void {
    const tagIndex = this._tags.indexOf(tag);
    if (tagIndex === -1)
      return;
    this._tags.splice(tagIndex, 1);
    const currentTag = this._tagsDiv.querySelector(`[data-tag="${tag}"]`);
    this._tagsDiv.removeChild(currentTag!);
    this._onTagRemoved.next(tag);
  }

  getTags(): string[] {
    return this._tags;
  }

  setTags(tags: string[]): void {
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
