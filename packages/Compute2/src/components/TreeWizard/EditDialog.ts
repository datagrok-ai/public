import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {fromEvent, Subject} from 'rxjs';
import {ItemMetadata} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/view/ViewCommunication';

type EditOptions = {
  title?: string,
  description?: string,
  isFavorite?: boolean,
  tags?: string[],
}

export class EditDialog extends DG.Dialog {
  private _onMetadataEdit = new Subject<EditOptions>();
  public onMetadataEdit = this._onMetadataEdit.asObservable();

  constructor(metadata?: ItemMetadata) {
    const dlg = ui.dialog({title: 'Save model state'});

    super(dlg.dart);

    let title = metadata?.title ?? '';
    let description = metadata?.description ?? '';
    let isFavorite = metadata?.isFavorite ?? false;
    const titleInput = ui.input.string('Title', {value: title, onValueChanged: (value) => title = value});

    const dummyInput = ui.input.string(' ', {value: ''});
    const tagsLine = DG.TagEditor.create();
    (metadata?.tags ?? []).forEach((tag: string) => {
      tagsLine.addTag(tag);
    });
    dummyInput.input.replaceWith(tagsLine.root);

    const addNewTag = () => {
      if (tagInput.value === '' ||
          // @ts-ignore
          tagsLine.tags.includes(tagInput.value)
      )
        return;

      tagsLine.addTag(tagInput.value);
      tagInput.value = '';
    };
    const tagInput = ui.input.string('Tag', {value: ''}).addOptions(ui.iconFA('plus', addNewTag, 'Add this tag'));

    const enterSub = fromEvent<KeyboardEvent>(tagInput.input, 'onkeydown').subscribe((ev) => {
      if (ev.key == 'Enter') {
        ev.stopPropagation();
        addNewTag();
      }
    });
    this.sub(enterSub);

    const descInput = ui.input.string('Description', {value: description,
      onValueChanged: (value) => description = value});

    const favInput = ui.input.bool('Is favorite', {value: isFavorite,
      onValueChanged: (value) => isFavorite = value});

    this.add(ui.form([
      titleInput,
      descInput,
      favInput,
      tagInput,
      dummyInput,
    ]));

    this.addButton('Save', async () => {
      if (tagInput.value !== '') {
        grok.shell.info(`Dialog has unsaved tags: ${tagInput.value}`);
        return;
      }

      const editOptions = {
        title: (title !== '') ? title : null,
        description: (description !== '') ? description : null,
        isFavorite,
        tags: tagsLine.tags.map((el) => el as any),
      } as EditOptions;

      this._onMetadataEdit.next(editOptions);
      this.close();
    });
  }
}
