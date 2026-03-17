import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {showGetRegionDialog} from './get-region';
import {ISequenceColumnInput} from '@datagrok-libraries/bio/src/utils/sequence-column-input';

/** A column input that filters to macromolecule columns and provides a
 *  "get region" button so users can extract a sub-region and use it instead. */
export class SequenceColumnInput implements ISequenceColumnInput {
  private readonly colInput: DG.InputBase<DG.Column | null>;

  private constructor(
    name: string,
    options: ui.input.IColumnInputInitOptions<DG.Column>,
  ) {
    const filter = options.filter;
    this.colInput = ui.input.column(name, {
      ...options,
      filter: (col: DG.Column) => {
        if (col.semType !== DG.SEMTYPE.MACROMOLECULE) return false;
        return filter ? filter(col) : true;
      },
    });

    const regionIcon = ui.iconFA('cut', () => this.onRegionIconClick(), 'Extract a region from the sequence');
    this.colInput.addOptions(regionIcon);
  }

  /** Creates a new SequenceColumnInput.
   *  @param name - Caption for the input.
   *  @param options - Same options as {@link ui.input.column}, table is required.
   *    The `filter` option is extended to always require semType === Macromolecule. */
  static create(
    name: string,
    options: ui.input.IColumnInputInitOptions<DG.Column>,
  ): SequenceColumnInput {
    return new SequenceColumnInput(name, options);
  }

  get root(): HTMLElement { return this.colInput.root; }
  get value(): DG.Column | null { return this.colInput.value; }
  set value(col: DG.Column | null) { this.colInput.value = col; }
  get inputBase(): DG.InputBase<DG.Column | null> { return this.colInput; }

  private onRegionIconClick(): void {
    const col = this.colInput.value;
    if (!col) {
      grok.shell.warning('Select a macromolecule column first.');
      return;
    }

    showGetRegionDialog(col, (regCol) => {
      this.colInput.value = regCol;
    });
  }
}
