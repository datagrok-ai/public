import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

export class ErrorMarkingPanel {
  public async init(t: DG.DataFrame, friendlyName: string) {
    const id = t.rows.match({friendly_name: friendlyName}).toDataFrame().get('id', 0);

    const error = await grok.dapi.logTypes.include('message,isError').filter(`id = "${id}"`).first();

    const isError = ui.input.bool('Is error', {value: error.isError});
    const comment = ui.input.string('Comment', {value: error.comment});

    const acc = DG.Accordion.create('ErrorInfo');
    acc.addPane('ErrorInfo', () => {
      const button = ui.buttonsInput([ui.bigButton('Save', () => {
        error.isError = isError.value!;
        error.comment = comment.value;
        grok.dapi.logTypes.save(error);
        grok.shell.info('Event type saved');
      })]);
      return ui.divV([isError.root, comment.root, button]);
    });

    return acc;
  }
}
