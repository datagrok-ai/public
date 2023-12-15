import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

export enum INs {
  receptor = 'receptor',
  gridParams = 'gridParams',
  gridParamsAuto = 'gridParamsAuto',
  ligandTable = 'ligandTable',
  ligandColumn = 'ligandColumn',
}

export class RunAutodockFuncEditor {
  inputs = new class {
    [INs.receptor]: DG.InputBase<DG.FileInfo | null>;
    [INs.gridParams]: DG.InputBase<DG.FileInfo | null>;
    [INs.gridParamsAuto]: DG.InputBase<boolean>;
    [INs.ligandTable]: DG.InputBase<DG.DataFrame | null>;
    [INs.ligandColumn]: DG.ChoiceInput<DG.Column | null>;
  }();

  constructor(
    private readonly call: DG.FuncCall,
  ) {
    const getDesc = (paramName: string) => this.call.inputParams[paramName].property.description;
    const getValue = (paramName: string) => this.call.inputParams[paramName].value;

    // const receptorFileProp = DG.Property.fromOptions(
    //   {name: 'receptor', caption: 'Receptor structure file', type: 'file', /* TODO: DG.TYPE.FILE */});
    const receptorFileProp = this.call.inputParams[INs.receptor].property;
    this.inputs[INs.receptor] = DG.InputBase.forProperty(receptorFileProp);
    this.inputs[INs.receptor].captionLabel.innerText = 'Receptor structure file';
    if (this.inputs[INs.receptor].value === undefined)
      this.inputs[INs.receptor].value = null;

    // const gridParamsProp = DG.Property.fromOptions(
    //   {name: 'gridParams', caption: 'Grid parameters file', type: 'file'});
    const gridParamsProp = this.call.inputParams[INs.gridParams].property;
    this.inputs[INs.gridParams] = DG.InputBase.forProperty(gridParamsProp);
    this.inputs[INs.gridParams].captionLabel.innerText = 'Grid parameters file';
    if (this.inputs[INs.gridParams].value === undefined)
      this.inputs[INs.gridParams].value = null;

    this.inputs[INs.gridParamsAuto] = ui.boolInput('Generate grid parameters', false,
      this.gridParamsAutoInputChanged.bind(this)) as DG.InputBase<boolean>;

    this.inputs[INs.ligandTable] = DG.InputBase.forProperty(this.call.inputParams[INs.ligandTable].property);
    this.inputs[INs.ligandTable].onChanged(this.ligandTableInputChanged.bind(this));
    // TODO: DG.InputBase.forProperty() return undefined for input of type column
    //const ligandColumnInput = DG.InputBase.forProperty(this.call.inputParams[INs.ligandColumn].property);
    // const ligandColumnInput = ui.columnInput('ligandColumn',
    //   getValue(INs.ligandTable) as DG.DataFrame, getValue(INs.ligandColumn) as DG.Column);
    const ligandTableVal: DG.DataFrame = getValue(INs.ligandTable);
    const ligandColumnVal: DG.Column = getValue(INs.ligandColumn);
    const ligandColumnChoices: DG.Column[] = this.getLigandColumnChoices();

    const ligandColumnInput = ui.choiceInput<DG.Column>(INs.ligandColumn, ligandColumnVal, ligandColumnChoices);
    this.inputs[INs.ligandColumn] = ligandColumnInput;

    // initial
    this.gridParamsFiBackup = getValue(INs.gridParams);
    this.gridParamsAutoInputChanged();
  }

  // -- Handle inputs' events --

  private getLigandColumnChoices(): DG.Column[] {
    const ligandTableVal: DG.DataFrame | null = this.inputs[INs.ligandTable].value;
    const ligandColumnListRes = !ligandTableVal ? [] : ligandTableVal.columns.toList()
      .filter((c) => [DG.SEMTYPE.MOLECULE, DG.SEMTYPE.MOLECULE3D].includes(c.semType));
    return ligandColumnListRes;
  }

  private gridParamsFiBackup: DG.FileInfo | null;

  private gridParamsAutoInputChanged(): void {
    this.inputs[INs.gridParams].enabled = !this.inputs[INs.gridParamsAuto].value;
    if (this.inputs[INs.gridParamsAuto].value) {
      this.gridParamsFiBackup = this.inputs[INs.gridParams].value;
      (async () => {
        const stubFi = (await _package.files.list('samples/', false, 'empty.tmp'))[0];
        this.inputs[INs.gridParams].value = stubFi;
      })();
    } else
      this.inputs[INs.gridParams].value = this.gridParamsFiBackup;
  }

  private ligandTableInputChanged(): void {
    //(this.inputs[INS.ligandColumn] as DG.Tabl).
    this.inputs[INs.ligandColumn].items = this.getLigandColumnChoices();
    const k = 11;
  }

  private getParams(): {} {
    return {
      [INs.receptor]: this.inputs[INs.receptor].value!,
      // TODO: [INs.gridParams]: !this.inputs[INs.gridParamsAuto].value ? this.inputs[INs.gridParams].value : null,
      [INs.gridParams]: this.inputs[INs.gridParams].value,
      [INs.ligandTable]: this.inputs[INs.ligandTable].value!,
      [INs.ligandColumn]: this.inputs[INs.ligandColumn].value!,
    };
  }

  // -- UI --

  public dialog(): void {
    const inputsForm = ui.inputs(Object.values(this.inputs), {style: {minWidth: '320px'}});
    ui.dialog({title: 'AutoDock'})
      .add(inputsForm)
      .onOK(async () => {
        try {
          const callParams = this.getParams();
          await this.call.func.prepare(callParams).call(true);
        } catch (err: any) {
          _package.handleErrorUI(err);
        }
      })
      .show();
  }

  /* For column hamburger menu */
  public widget(): DG.Widget {
    const inputsForm = ui.inputs(Object.entries(this.inputs)
      .filter(([inputName, _input]) => !['ligandTable', 'ligandColumn'].includes(inputName))
      .map(([_inputName, input]) => input));
    const doBtn = ui.button('AutoDock', async () => {
      try {
        const callParams = this.getParams();
        await this.call.func.prepare(callParams).call(true);
      } catch (err: any) {
        _package.handleErrorUI(err);
      }
    });
    return DG.Widget.fromRoot(ui.divV([inputsForm, ui.div(doBtn)]));
  }
}
