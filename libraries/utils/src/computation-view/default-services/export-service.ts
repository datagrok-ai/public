/* eslint-disable valid-jsdoc */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ExportService, ComputationViewStateService} from '../common/service-interfaces';
import {DefaultComputationViewStateService} from './state-service';

export class DefaultExportService implements ExportService {
  constructor(
    private state: ComputationViewStateService = new DefaultComputationViewStateService(),
  ) {}

  filename(format: string): string {
    return `${this.state.name} - ${new Date().toLocaleString()}.${this.supportedExportExtensions[format]}`;
  }

  /** Override to provide supported export formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportFormats(): string[] {
    return [
      'Excel', 'PDF', 'CSV',
    ];
  }

  /** Override to provide custom file extensions for exported formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportExtensions(): Record<string, string> {
    return {
      'Excel': 'xlsx',
      'PDF': 'pdf',
      'CSV': 'csv',
    };
  }

  /** Override to provide custom export. */
  async export(format: string): Promise<Blob> {
    throw new Error(`Format "${format}" is not supported.`);
  }
}
