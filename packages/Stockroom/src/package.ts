/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {domains} from '@datagrok-libraries/u2/src/dg/index.js';
export * from './package.g';

// the skins of everything the u2 app renders (tokens first: every other sheet reads them)
import '@datagrok-libraries/u2/css/tokens.css';
import '@datagrok-libraries/u2/css/elements.css';
import '@datagrok-libraries/u2/css/buttons.css';
import '@datagrok-libraries/u2/css/inputs.css';
import '@datagrok-libraries/u2/css/number.css';
import '@datagrok-libraries/u2/css/date.css';
import '@datagrok-libraries/u2/css/choice.css';
import '@datagrok-libraries/u2/css/combobox.css';
import '@datagrok-libraries/u2/css/tags.css';
import '@datagrok-libraries/u2/css/typeahead.css';
import '@datagrok-libraries/u2/css/entity.css';
import '@datagrok-libraries/u2/css/form.css';
import '@datagrok-libraries/u2/css/list.css';
import '@datagrok-libraries/u2/css/menu.css';
import '@datagrok-libraries/u2/css/splitter.css';
import '@datagrok-libraries/u2/css/async.css';
import '@datagrok-libraries/u2/css/dialog.css';
import '@datagrok-libraries/u2/css/notify.css';
import '@datagrok-libraries/u2/css/tooltip.css';
import '@datagrok-libraries/u2/css/badge.css';
import '@datagrok-libraries/u2/css/file.css';
import '@datagrok-libraries/u2/css/domain.css';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Stockroom
//description: Chemical stockroom on the GHS classification — the zero-code app over databases/stockroom/schema.json
//tags: app
//meta.icon: images/flask.svg
//output: view result
export async function stockroomApp(): Promise<DG.ViewBase> {
  return (await domains.table('stockroom.substance')).app();
}
