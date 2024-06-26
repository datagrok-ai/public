import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerWidthMode} from '../utils/cell-renderer-consts';

import {_package} from '../package';

export class PackageSettingsEditorWidget extends DG.Widget {
  maxMonomerLengthProp: DG.Property;
  tooltipWebLogo: DG.Property;
  defaultSeparator: DG.Property;

  constructor(propList: DG.Property[]) {
    super(ui.div([], {}));

    const props: { [name: string]: DG.Property } =
      Object.assign({}, ...propList.map((p) => ({[p.name]: p})));

    this.maxMonomerLengthProp = props['MaxMonomerLength'];
    this.tooltipWebLogo = props['TooltipWebLogo'];
    this.defaultSeparator = props['DefaultSeparator'];
  }

  async init(): Promise<void> {
    const monomerWidthModeInput = ui.choiceInput('Monomer width mode',
      _package.properties.MonomerWidthMode,
      [MonomerWidthMode.short, MonomerWidthMode.long],
      (value: MonomerWidthMode) => {
        _package.properties.MonomerWidthMode = value;
      });

    const maxMonomerLengthInput = ui.intInput('Max monomer length',
      _package.properties.MaxMonomerLength,
      (value: number) => {
        // Handle user changed value
        _package.properties.MaxMonomerLength = value;
      });

    const tooltipWebLogoInput = ui.boolInput('Tooltip WebLogo',
      _package.properties.TooltipWebLogo,
      (value: boolean) => {
        _package.properties.TooltipWebLogo = value;
      });

    const defaultSeparatorInput = ui.choiceInput('Default Separator',
      _package.properties.DefaultSeparator, ['.', '/', '-'],
      (value: string) => {
        _package.properties.DefaultSeparator = value;
      });

    this.root.appendChild(ui.form([
      monomerWidthModeInput,
      maxMonomerLengthInput,
      tooltipWebLogoInput,
      defaultSeparatorInput,
    ]));
  }
}
