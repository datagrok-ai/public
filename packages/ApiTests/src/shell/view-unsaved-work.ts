import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {until} from '../ai/helpers';

const nodeSkip = typeof process !== 'undefined' ? 'NodeJS environment' : undefined;

let count = 0;

/** A view with changes while `work` is set. */
class EditorView extends DG.ViewBase {
  work: DG.UnsavedWork | null = null;

  constructor() {
    super();
    this.name = `Unsaved work ${++count}`;
  }

  get unsavedWork(): DG.UnsavedWork | null { return this.work; }
}

const isOpen = (v: DG.ViewBase): boolean => Array.from(grok.shell.views).some((x) => x.name === v.name);
const question = (v: DG.ViewBase): DG.Dialog | undefined =>
  DG.Dialog.getOpenDialogs().find((d) => d.title === `Close ${v.name}?`);

async function open(): Promise<EditorView> {
  const view = new EditorView();
  grok.shell.addView(view);
  await until(() => grok.shell.v.name === view.name && view.root.isConnected);
  return view;
}

/** Closes the current view as a user does: with its tab's close button, or the ribbon's in the simple mode. */
function userClose(view: DG.ViewBase): void {
  const tabButton = document.querySelector(`.tab-handle[name="view-handle: ${view.name}"] .tab-handle-close-button`);
  if (tabButton)
    tabButton.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
  else
    (document.querySelector('.d4-view-ribbon-controls [name="icon-times"]') as HTMLElement).click();
}

async function ask(view: EditorView, work: Partial<DG.UnsavedWork> = {}): Promise<DG.Dialog> {
  view.work = {changes: ['Edited the note', 'Renamed the note'], save: async () => true, ...work};
  userClose(view);
  await until(() => question(view) != null);
  return question(view)!;
}

category('View: Unsaved work', () => {
  test('a user\'s close asks while the view has unsaved work', async () => {
    const view = await open();
    try {
      const dlg = await ask(view, {saveLink: 'save the note'});
      const text = dlg.root.textContent!;
      expect(text.includes('Unsaved changes') && text.includes('Edited the note, renamed the note.'), true, text);
      expect(text.includes('To keep these changes, save the note first'), true, text);
      dlg.getButton('CANCEL').click();
      await until(() => question(view) == null);
      expect(isOpen(view), true, 'CANCEL keeps the view');

      view.work = null;
      userClose(view);
      await until(() => !isOpen(view));
      expect(question(view) == null, true, 'nothing to lose, nothing asked');
    } finally {
      question(view)?.close();
      view.close();
    }
  }, {skipReason: nodeSkip});

  test('CLOSE discards the changes; the link saves them first', async () => {
    const view = await open();
    let saved = false;
    let saves = 0;
    const work = {save: async () => {
      saves++;
      return saved;
    }};
    const link = (dlg: DG.Dialog): HTMLElement => dlg.root.querySelector('.grok-close-guard-save a')!;
    try {
      const refused = link(await ask(view, work));
      expect(refused.textContent, 'save the changes', 'the link without saveLink');
      refused.click();
      await until(() => saves === 1 && question(view) == null);
      expect(isOpen(view), true, 'what was not saved stays open');

      saved = true;
      link(await ask(view, work)).click();
      await until(() => !isOpen(view));
      expect(saves, 2);

      const other = await open();
      try {
        (await ask(other, work)).getButton('CLOSE').click();
        await until(() => !isOpen(other));
        expect(saves, 2, 'CLOSE saves nothing');
      } finally {
        question(other)?.close();
        other.close();
      }
    } finally {
      question(view)?.close();
      view.close();
    }
  }, {skipReason: nodeSkip});

  test('a close from code, and one whose unsaved work cannot be read, do not ask', async () => {
    const view = await open();
    view.work = {changes: ['Edited the note'], save: async () => true};
    view.close();
    expect(isOpen(view), false);
    expect(question(view) == null, true);

    const broken = await open();
    Object.defineProperty(broken, 'unsavedWork', {get: () => {throw new Error('unsaved work of a test view');}});
    try {
      userClose(broken);
      await until(() => !isOpen(broken));
      expect(question(broken) == null, true);
    } finally {
      broken.close();
    }
  }, {skipReason: nodeSkip});
}, {owner: 'drizhinashvili@datagrok.ai'});
