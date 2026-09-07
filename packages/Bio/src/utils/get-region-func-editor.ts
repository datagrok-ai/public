import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_package} from '../package';

export interface GetRegionParams {
  table: DG.DataFrame,
  sequence: DG.Column<string>,
  start: string | null,
  end: string | null,
  /** Name for the column with sequence of the region  */ name: string | null,
}

export interface SeqRegion {
  name: string,
  description: string,
  start: string,
  end: string,
}

const regionLabel = (r: SeqRegion): string => `${r.name}: ${r.start}-${r.end}`;

export class GetRegionFuncEditor {
  inputs = new class {
    table: DG.InputBase<DG.DataFrame | null>;
    sequence: DG.InputBase<DG.Column | null>;
    /** Items are region labels (see {@link regionLabel}); the selected region is looked up in {@link regionList}. */
    region: DG.ChoiceInput<string | null>;
    start: DG.ChoiceInput<string | null>;
    end: DG.ChoiceInput<string | null>;
    name: DG.InputBase<string>;
  }();

  private regionList: SeqRegion[] = [];

  constructor(
    private readonly call: DG.FuncCall,
    private readonly seqHelper: ISeqHelper,
  ) {
    const getDesc = (paramName: string) => this.call.inputParams[paramName].property.description;

    this.inputs.table = ui.input.table('Table', {value: this.call.inputParams['table'].value ?? grok.shell.tv.dataFrame});

    //@formatter:off
    const seqColValue = this.call.inputParams['sequence'].value ??
      this.inputs.table.value!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    this.inputs.sequence = ui.input.column('Sequence', {table: grok.shell.tv.dataFrame, value: seqColValue,
      onValueChanged: this.sequenceInputChanged.bind(this), filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE});
    this.inputs.start = ui.input.choice<string>('Start', {onValueChanged: this.startInputChanged.bind(this)});
    this.inputs.end = ui.input.choice<string>('End', {onValueChanged: this.endInputChanged.bind(this)});

    this.inputs.region = ui.input.choice<string>('Region', {nullable: true,
      onValueChanged: this.regionInputChanged.bind(this)});

    this.inputs.name = ui.input.string('Column name', {value: this.getDefaultName(),
      onValueChanged: this.nameInputChanged.bind(this), clearIcon: true});
    this.inputs.name.onInput.subscribe(() => this.nameInputInput.bind(this)); // To catch clear event
    //@formatter:on

    // tooltips
    for (const paramName in this.call.inputParams) {
      // @ts-ignore
      ui.tooltip.bind(this.inputs[paramName].captionLabel, getDesc(paramName));
    }

    // initial
    this.sequenceInputChanged();
  }

  private sequenceInputChanged(): void {
    this.updateRegionItems();
    this.updateStartEndInputItems();
    this.updateRegion();
    this.updateNameInput();
  }

  private fixRegion: boolean = false;

  private regionInputChanged(): void {
    this.fixRegion = true;
    try {
      const reg = this.getSelectedRegion();
      if (reg) {
        this.inputs.start.value = reg.start;
        this.inputs.end.value = reg.end;
      } else {
        const posList = this.inputs.start.items;
        this.inputs.start.value = posList[0] ?? null;
        this.inputs.end.value = posList[posList.length - 1] ?? null;
      }
    } finally {
      this.fixRegion = false;
    }
  }

  private startInputChanged(): void {
    this.updateRegion();
    this.updateNameInput();
  }

  private endInputChanged(): void {
    this.updateRegion();
    this.updateNameInput();
  }

  private nameInputChanged(): void {
    if (!this.defaultNameUpdating)
      this.defaultName = false;
  }

  private nameInputInput(): void {
    if (!this.inputs.name.value) {
      this.defaultName = true;
      this.inputs.name.input.focus();
    }
  }

  private updateStartEndInputItems(): void {
    const seqCol = this.inputs.sequence.value;
    const posList: string[] = seqCol ? this.seqHelper.getSeqHandler(seqCol).posList : [];

    this.inputs.start.items = posList;
    this.inputs.end.items = posList;
    this.inputs.start.value = posList[0] ?? null;
    this.inputs.end.value = posList[posList.length - 1] ?? null;
  }

  private updateRegionItems(): void {
    const seqCol = this.inputs.sequence.value;
    // Read from .annotations first (new system), fall back to .regions (legacy)
    let regionList: SeqRegion[] | null = null;
    const annotationsTag: string | null = seqCol ? seqCol.getTag(bioTAGS.annotations) : null;
    if (annotationsTag) {
      try {
        const annotations = JSON.parse(annotationsTag);
        const structAnnots = annotations.filter((a: any) => a.category === 'structure' && a.start && a.end);
        if (structAnnots.length > 0) {
          regionList = structAnnots.map((a: any) => ({
            name: a.name, description: a.description ?? '', start: a.start, end: a.end,
          }));
        }
      } catch { /* ignore parse errors */ }
    }
    if (!regionList) {
      const regionsTagTxt: string | null = seqCol ? seqCol.getTag(bioTAGS.regions) : null;
      regionList = regionsTagTxt ? JSON.parse(regionsTagTxt) : null;
    }

    this.regionList = regionList ?? [];
    this.inputs.region.items = this.regionList.map(regionLabel);
    if (regionList != null)
      this.inputs.region.root.style.removeProperty('display');
    else
      this.inputs.region.root.style.display = 'none';
  }

  private getSelectedRegion(): SeqRegion | null {
    const label = this.inputs.region.value;
    return label ? this.regionList.find((r) => regionLabel(r) === label) ?? null : null;
  }

  /** Selects the region matching the current start/end positions, or none. */
  private updateRegion(): void {
    if (this.fixRegion) return;
    const startPos = this.inputs.start.value;
    const endPos = this.inputs.end.value;
    const reg = this.regionList.find((r) => r.start === startPos && r.end === endPos);
    const label = reg ? regionLabel(reg) : null;
    if (this.inputs.region.value !== label)
      this.inputs.region.value = label;
  }

  private defaultName: boolean = true;
  private defaultNameUpdating: boolean = false;

  private updateNameInput(): void {
    this.defaultNameUpdating = true;
    try {
      if (this.defaultName) this.inputs.name.value = this.getDefaultName();
    } finally {
      this.defaultNameUpdating = false;
    }
  }

  private getDefaultName(): string {
    const reg = this.getSelectedRegion();
    const seqCol: DG.Column<string> = this.inputs.sequence.value!;
    const startPos: string = this.inputs.start.value ?? '';
    const endPos: string = this.inputs.end.value ?? '';

    return reg != null ? `${seqCol.name}(${reg.name}): ${reg.start}-${reg.end}` :
      `${seqCol?.name}: (${startPos}-${endPos})`;
  }

  public getParams(): GetRegionParams {
    return {
      table: this.inputs.table.value!,
      sequence: this.inputs.sequence.value!,
      start: this.getStart(),
      end: this.getEnd(),
      name: this.getName(),
    };
  }

  private getStart(): string | null {
    return this.inputs.start.value;
  }

  private getEnd(): string | null {
    return this.inputs.end.value;
  }

  private getName(): string | null {
    const str = this.inputs.name.stringValue;
    return str == '' ? null : str;
  }

  // -- History --

  /** Serializes the region selection (start/end position codes + column name) for dialog history. */
  public getStringInput(): string {
    return JSON.stringify({
      start: this.getStart(),
      end: this.getEnd(),
      name: this.inputs.name.stringValue,
    });
  }

  /** Restores a region selection from history, applying start/end only if the positions exist in the
   * current sequence column (positions depend on the chosen column). */
  public applyStringInput(input: string): void {
    try {
      const parsed = JSON.parse(input);
      const seqCol = this.inputs.sequence.value;
      const posList: string[] = seqCol ? this.seqHelper.getSeqHandler(seqCol).posList : [];
      if (parsed.start != null && posList.includes(parsed.start))
        this.inputs.start.value = parsed.start;
      if (parsed.end != null && posList.includes(parsed.end))
        this.inputs.end.value = parsed.end;
      // set the name last: changing start/end resets it to the default via updateNameInput
      if (parsed.name)
        this.inputs.name.value = parsed.name;
    } catch (e: any) {
      _package.logger.error(e instanceof Error ? e.message : e?.toString());
    }
  }

  // -- UI --

  /** Full inputs form (table + sequence + region/start/end/name), for hosting inside the canonical
   * FuncCallParamsEditor. */
  public getEditorForm(): HTMLElement {
    return ui.inputs(Object.values(this.inputs), {style: {minWidth: '320px'}});
  }

  public dialog(): void {
    const inputsForm = ui.inputs(Object.values(this.inputs), {style: {minWidth: '320px'}});
    ui.dialog({title: 'Get Region'})
      .add(inputsForm)
      .onOK(() => {
        (async () => {
          const callParams = this.getParams();
          await this.call.func.prepare(callParams).call(true);
        })()
          .catch((err: any) => { _package.handleErrorUI(err); });
      })
      .show();
  }

  public widget(): DG.Widget {
    const inputsForm = ui.inputs(Object.entries(this.inputs)
      .filter(([inputName, _input]) => !['table', 'sequence'].includes(inputName))
      .map(([_inputName, input]) => input));
    const doBtn = ui.button('Get Region', () => {
      (async () => {
        const callParams = this.getParams();
        await this.call.func.prepare(callParams).call(true);
      })()
        .catch((err: any) => { _package.handleErrorUI(err); });
    });
    return DG.Widget.fromRoot(ui.divV([inputsForm, ui.div(doBtn)]));
  }
}
