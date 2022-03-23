import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";

export class ErrorMarkingPanel {
  public async init(t: DG.DataFrame, friendlyName: string) {
    let id = t.rows.match({ friendly_name: friendlyName}).toDataFrame().get('id', 0);

    let error = await grok.dapi.logTypes.include('message,isError').filter(`id = "${id}"`).first();

    let isError = ui.boolInput('Is error', error.isError);
    let comment = ui.stringInput('Comment', error.comment)

    let acc = DG.Accordion.create('ErrorInfo');
    acc.addPane('ErrorInfo', () => {

      let button = ui.buttonsInput([ui.bigButton('Save', () => {
        error.isError = isError.value;
        error.comment = comment.value;
        grok.dapi.logTypes.save(error);
        grok.shell.info('Event type saved');
      })]);
      return ui.divV([isError.root, comment.root, button]);
    });

    return acc;
  }
}