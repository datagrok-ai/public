/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class SystemStatusWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.div());

    this.root.appendChild(ui.divText('Houston, we have a problem!'));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'current time: ');
  }
}

