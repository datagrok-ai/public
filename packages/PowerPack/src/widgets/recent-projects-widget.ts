import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class RecentProjectsWidget extends DG.Widget {
  caption: string = 'Recent projects';

  constructor() {
    super(ui.div());
    this.root.appendChild(ui.divText('Recent projects'));
  }
}