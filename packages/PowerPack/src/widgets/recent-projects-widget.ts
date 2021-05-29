import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class RecentProjectsWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.div());
    this.root.appendChild(ui.divText('Recent projects'));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Recent projects');
  }
}