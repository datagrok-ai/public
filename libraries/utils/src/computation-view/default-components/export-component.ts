/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Component, ExportService} from '../common/service-interfaces';
import {tokens} from '../common/inject-tokens';

export class DefaultExportComponent implements Component {
  public root: HTMLElement = ui.div();

  public static inject = [
    tokens.exportService,
  ] as const;

  constructor(
    private exportService: ExportService,
  ) {
    this.render();
  }

  render(): void {
    this.root = ui.comboPopup(
      ui.icons.save(null, 'Export'),
      this.exportService.supportedExportFormats,
      async (format: string) =>
        DG.Utils.download(this.exportService.filename(format), await this.exportService.export(format)),
    );
  }
}
