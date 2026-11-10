import * as ui from 'datagrok-api/ui';
import {ignoreTab} from '../../composables/use-keyboard';

export function confirmAction(confirmationMessage: string | undefined, run: () => void) {
  if (confirmationMessage) {
    ui.dialog(`Action confirmation`)
      .add(ui.markdown(confirmationMessage))
      .onOK(run)
      .show({center: true, modal: true});
  } else
    run();
}

export function confirmRemove(name: string, remove: () => void) {
  ui.dialog(`Remove confirmation`)
    .add(ui.divText(`Do you want to remove "${name}"? Its unsaved inputs and results will be lost.`))
    .onOK(remove)
    .show({center: true, modal: true});
}

export function confirmSubtreeUpdate(update: (includeInfo: boolean) => void) {
  const includeInfoInput = ui.input.bool('Also reset info-typed inputs', {value: false});
  ignoreTab(ui.dialog(`Update confirmation`)
    .add(ui.markdown(`Do you want to update input values to consistent ones and rerun substeps? You will lose inconsistent values.`))
    .add(includeInfoInput.root)
    .onOK(() => update(includeInfoInput.value ?? false)))
    .show({center: true, modal: true});
}

export function confirmReturn(returnResult: () => void) {
  ui.dialog(`Action confirmation`)
    .add(ui.markdown(`Close this workflow and return its results`))
    .onOK(returnResult)
    .show({center: true, modal: true});
}
