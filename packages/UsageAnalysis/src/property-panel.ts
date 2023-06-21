import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaViewer} from './viewers/abstract/ua-viewer';

export class PropertyPanel {
  private acc: DG.Accordion;
  private root: HTMLElement;

  public constructor(entity: DG.Entity | null, entityViewer: UaViewer | null, viewers: UaViewer[], header: string,
    accordionKey: string) {
    this.acc = DG.Accordion.create(accordionKey);

    if (entity != null) {
      this.acc.addPane('Details', () => {
        const res = ui.divV([
          ui.render(entity),
        ]);

        if (entityViewer)
          res.append(entityViewer.root);

        return res;
      });
    }

    for (const viewer of viewers)
      this.acc.addPane(viewer.name, () => viewer.root);


    this.root = ui.block([
      ui.divV([
        ui.h1(header),
        this.acc.root,
      ]),
    ]);
  }

  public getRoot(): HTMLElement {
    return this.root;
  }
}
