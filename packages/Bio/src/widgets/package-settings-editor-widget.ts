import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

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
    const maxMonomerLengthInput = ui.input.int('Max Monomer Length', {
      value: _package.properties.maxMonomerLength!,
      nullable: true, min: 0,
      onValueChanged: (value) => {
        // Handle user changed value
        _package.properties.maxMonomerLength = value ?? null;
      },
      tooltipText: this.maxMonomerLengthProp.description,
    });

    const tooltipWebLogoInput = ui.input.bool('Tooltip WebLogo', {
      value: _package.properties.tooltipWebLogo,
      nullable: false,
      onValueChanged: (value) => {
        _package.properties.tooltipWebLogo = value;
      },
      tooltipText: this.tooltipWebLogo.description,
    });

    const defaultSeparatorInput = ui.input.choice('Default Separator', {
      value: _package.properties.defaultSeparator,
      items: ['.', '/', '-'],
      onValueChanged: (value) => {
        _package.properties.defaultSeparator = value;
      },
      tooltipText: this.defaultSeparator.description,
    });

    this.root.appendChild(ui.form([
      maxMonomerLengthInput,
      tooltipWebLogoInput,
      defaultSeparatorInput,
    ]));
  }
}
