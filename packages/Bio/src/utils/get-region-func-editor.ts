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

export class GetRegionFuncEditor {
  inputs = new class {
    table: DG.InputBase<DG.DataFrame | null>;
    sequence: DG.InputBase<DG.Column | null>;
    region: DG.InputBase<SeqRegion>;
    start: DG.InputBase<string>;
    end: DG.InputBase<string>;
    name: DG.InputBase<string>;
  }();

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
    this.inputs.start = ui.input.choice('Start', {onValueChanged: this.startInputChanged.bind(this)}) as unknown as DG.InputBase<string>;
    this.inputs.end = ui.input.choice('End', {onValueChanged: this.endInputChanged.bind(this)}) as unknown as DG.InputBase<string>;

    this.inputs.region = ui.input.choice<SeqRegion>('Region', {value: null as unknown as SeqRegion, items: [],
      onValueChanged: this.regionInputChanged.bind(this)}) as DG.InputBase<SeqRegion>;

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
    const seqCol = this.inputs.sequence.value;
    const sh = seqCol ? this.seqHelper.getSeqHandler(seqCol) : null;
    this.updateRegionItems();
    this.updateStartEndInputItems();
    this.updateRegion(true);
    this.updateNameInput();
  }

  private fixRegion: boolean = false;

  private regionInputChanged(): void {
    this.fixRegion = true;
    try {
      const regJsonStr = this.inputs.region.stringValue;
      const reg: SeqRegion | null = regJsonStr ? JSON.parse(regJsonStr) as SeqRegion : null;

      if (reg !== null) {
        this.inputs.start.value = reg?.start;
        this.inputs.end.value = reg?.end;
      } else {
        const sh = this.seqHelper.getSeqHandler(this.inputs.sequence.value!);
        this.inputs.start.value = sh.posList[0];
        this.inputs.end.value = sh.posList[sh.posList.length - 1];
      }
    } finally {
      this.fixRegion = false;
    }
  }

  private startInputChanged(): void {
    this.updateRegion(false);
    this.updateNameInput();
  }

  private endInputChanged(): void {
    this.updateRegion(false);
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
    const sh = seqCol ? this.seqHelper.getSeqHandler(seqCol) : null;

    const startSE = (this.inputs.start.input as HTMLSelectElement);
    const endSE = (this.inputs.end.input as HTMLSelectElement);
    for (let i = startSE.options.length - 1; i >= 0; --i) startSE.options.remove(i);
    for (let i = endSE.options.length - 1; i >= 0; --i) endSE.options.remove(i);
    for (const pos of sh?.posList ?? []) {
      const startPosOE = document.createElement('option');
      const endPosOE = document.createElement('option');
      startPosOE.text = endPosOE.text = pos;
      startPosOE.value = endPosOE.value = pos;
      startSE.options.add(startPosOE);
      endSE.options.add(endPosOE);
    }
    startSE.value = sh?.posList[0] ?? '';
    endSE.value = sh?.posList[sh?.posList.length - 1] ?? '';
  }

  private updateRegionItems(): void {
    const seqCol = this.inputs.sequence.value;
    const regionsTagTxt: string | null = seqCol ? seqCol.getTag(bioTAGS.regions) : null;
    const regionList: SeqRegion[] | null = regionsTagTxt ? JSON.parse(regionsTagTxt) : null;

    const regionSE = (this.inputs.region.input as HTMLSelectElement);
    for (let i = regionSE.options.length - 1; i >= 0; --i) regionSE.options.remove(i);

    const nullOE = document.createElement('option');
    nullOE.text = '';
    nullOE.value = JSON.stringify(null);
    regionSE.options.add(nullOE);

    if (regionList != null) {
      this.inputs.region.root.style.removeProperty('display');
      for (const region of regionList) {
        const regionOE = document.createElement('option');
        regionOE.text = `${region.name}: ${region.start}-${region.end}`;
        regionOE.value = JSON.stringify(region);
        regionSE.options.add(regionOE);
      }
    } else {
      this.inputs.region.root.style.display = 'none';
    }
  }

  private updateRegion(reset: boolean): void {
    const startPos: string = this.inputs.start.stringValue ?? '';
    const endPos: string = this.inputs.end.stringValue ?? '';

    if (!this.fixRegion) {
      const regionSE = (this.inputs.region.input as HTMLSelectElement);
      regionSE.selectedIndex = -1;
      for (let i = regionSE.options.length - 1; i >= 0; --i) {
        const regionOE = regionSE.options[i];
        const reg: SeqRegion = JSON.parse(regionOE.value);
        if (reg && startPos === reg.start && endPos === reg.end) regionSE.selectedIndex = i;
      }
    }
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
    const regionJsonStr = this.inputs.region.stringValue;
    const reg: SeqRegion | null = regionJsonStr ? JSON.parse(regionJsonStr) : null;

    const seqCol: DG.Column<string> = this.inputs.sequence.value!;

    const startPos: string = this.inputs.start.stringValue ?? '';
    const endPos: string = this.inputs.end.stringValue ?? '';

    return reg != null ? `${seqCol.name}(${reg.name}): ${reg.start}-${reg.end}` :
      `${seqCol?.name}: (${startPos}-${endPos})`;
  }

  private getParams(): {} {
    return {
      table: this.inputs.table.value!,
      sequence: this.inputs.sequence.value!,
      start: this.getStart(),
      end: this.getEnd(),
      name: this.getName(),
    };
  }

  private getStart(): string | null {
    return this.inputs.start.stringValue;
  }

  private getEnd(): string | null {
    return this.inputs.end.stringValue;
  }

  private getName(): string | null {
    const str = this.inputs.name.stringValue;
    return str == '' ? null : str;
  }

  // -- UI --

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
