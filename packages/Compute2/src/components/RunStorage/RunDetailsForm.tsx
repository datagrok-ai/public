import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {Place, notifyRunsChanged, parseTags, placeName} from './run-storage';

/** Edits a saved run's name, description and tags; a name stays unique in the space the run is shown in. */
export const RunDetailsForm = Vue.defineComponent({
  name: 'RunDetailsForm',
  props: {
    run: {type: DG.FuncCall, required: true},
    place: {type: Object as Vue.PropType<Place>, default: null},
  },
  emits: {
    // with the saved run, or nothing on Cancel
    'done': (_saved?: DG.FuncCall) => true,
  },
  setup(props, {emit}) {
    const {name, description, tags} = historyUtils.getRunMetadata(props.run);
    const editName = Vue.ref(name);
    const editDescription = Vue.ref(description);
    const editTags = Vue.ref(tags.join(', '));
    const isBusy = Vue.ref(false);
    const error = Vue.ref<string>();

    const save = async () => {
      if (isBusy.value)
        return;
      const name = editName.value.trim();
      if (!name) {
        error.value = 'Enter a name';
        return;
      }
      isBusy.value = true;
      try {
        if (props.place) {
          const conflict = await historyUtils.findNameConflict(name, props.place.id);
          if (conflict.other ?? conflict.runs.find((item) => item.id !== props.run.id)) {
            error.value = `${name} is already in ${placeName(props.place)}. Choose another name.`;
            return;
          }
        }
        // the run shown keeps its details until the edit is saved
        const edited = await historyUtils.shallowLoadRun(props.run.id);
        historyUtils.setRunMetadata(edited,
          {name, description: editDescription.value.trim(), tags: parseTags(editTags.value)});
        await historyUtils.updateRunMeta(edited);
        emit('done', edited);
        notifyRunsChanged({runIds: [props.run.id]});
      } catch (e: any) {
        error.value = `Couldn't save the details: ${e?.message ?? e}`;
      } finally {
        isBusy.value = false;
      }
    };

    return () => (
      <div class='c2-runs-form'>
        { error.value && <div class='c2-runs-error'>{error.value}</div> }
        <label class='c2-runs-field'><span>Name</span>
          <input class='ui-input-editor' value={editName.value}
            onInput={(e: Event) => editName.value = (e.target as HTMLInputElement).value}/></label>
        <label class='c2-runs-field'><span>Description</span>
          <textarea class='ui-input-editor' rows={5} value={editDescription.value}
            onInput={(e: Event) => editDescription.value = (e.target as HTMLTextAreaElement).value}/></label>
        <label class='c2-runs-field'><span>Tags</span>
          <input class='ui-input-editor' placeholder='Comma-separated' value={editTags.value}
            onInput={(e: Event) => editTags.value = (e.target as HTMLInputElement).value}/></label>
        <div class='c2-runs-buttons'>
          <div class='c2-runs-button-group'>
            <button class='ui-btn' onClick={() => emit('done')}>Cancel</button>
            <button class={['ui-btn', 'ui-btn-ok', isBusy.value ? 'd4-disabled' : '']} onClick={save}>
              Save details</button>
          </div>
        </div>
      </div>
    );
  },
});
