import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Dialog, Form, TextArea, computed, notify, signal, subscribe} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

export interface CampaignSession {
  readonly view: DG.TableView;
  readonly changed: ReadonlySignal<string | null>;
  refresh(): void;
}

export function openCampaignSession(df: DG.DataFrame): CampaignSession {
  const view = grok.shell.addTableView(df);
  const host = ViewHost.of(view);
  return host.build(() => {
    const changed = signal<string | null>(null);
    const state = host.signals;
    host.effect(() => {
      const grid = state.grid.value;
      if (grid !== null)
        subscribe(grid.onCellValueEdited, () => changed.value = 'Unsaved edits');
    });
    host.subscribe(grok.events.onFileImportRequest, () => notify.info('Import into the campaign'));
    const refresh = () => {
      changed.value = null;
      notify.info('Refreshed');
    };
    const addNote = () => {
      const note = new TextArea({label: 'Note', required: true});
      host.show(new Dialog('Add a note', {name: 'campaign-note'}).add(new Form().add(note))
        .onOK(() => {
          notify.info(`Noted: ${note.value}`);
        }));
    };
    host.ribbon({hideSave: true, name: 'campaign-actions',
      actions: [
        {icon: 'sync', tooltip: 'Refresh the campaign', run: refresh, hint: changed, hintDetail: 'New changes'},
        {icon: 'comment', tooltip: 'Add a note', text: 'Add a note...', run: addNote},
      ],
      buttons: [{text: 'Submit...', primary: true, run: () => {
        notify.info(`Submitting ${df.selection.trueCount} molecules`);
      }, disabled: computed(() => state.selection.value?.trueCount ? null : 'Select the molecules to submit')}],
      status: computed(() => `${state.selection.value?.trueCount ?? 0} of ${state.rowCount.value} selected`),
    });
    return {view, changed, refresh};
  });
}
