//api: DG.ViewBase.unsavedWork, DG.UnsavedWork
// A view that names its unsaved changes is asked about when the user closes it: CLOSE discards them, the link
// saves them first. Edit the note and close its tab.

class NoteView extends DG.ViewBase {
  constructor() {
    super();
    this.name = 'Note';
    this.saved = 'Edit me, then close the tab';
    this.editor = ui.input.textArea('Note', {value: this.saved});
    this.root.append(this.editor.root);
  }

  get unsavedWork() {
    if (this.editor.value === this.saved)
      return null;
    return {
      changes: ['Edited the note'],
      saveLink: 'save the note',
      save: async () => {
        this.saved = this.editor.value;
        grok.shell.info(`Saved: ${this.saved}`);
        return true;
      },
    };
  }
}

grok.shell.addView(new NoteView());
