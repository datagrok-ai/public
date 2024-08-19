import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';

export async function saveCampaignDialog(id: string): Promise<string> {
  return new Promise<string>(async (resolve) => {
    const campaignNameInput = ui.input.string('Campaign name', {value: ''});

    function onOkProxy() {
      if (!campaignNameInput.value || campaignNameInput.value === '')
        resolve(id);
      else
        resolve(campaignNameInput.value);
    }
    ui.dialog('Save Campaign')
      .add(ui.div(campaignNameInput.root))
      .onOK(onOkProxy)
      .show();
  });
}
