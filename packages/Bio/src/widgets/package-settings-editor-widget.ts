import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

export class PackageSettingsEditorWidget extends DG.Widget {
  maxMonomerLengthProp: DG.Property;

  constructor(propList: DG.Property[]) {
    super(ui.div([], {}));

    const props: { [name: string]: DG.Property } =
      Object.assign({}, ...propList.map((p) => ({[p.name]: p})));

    this.maxMonomerLengthProp = props['MaxMonomerLength'];
  }

  async init(): Promise<void> {
    const maxMonomerLengthInput = ui.intInput('Max monomer length',
      _package.properties.maxMonomerLength,
      (value: number) => {
        // Handle user changed value
        _package.properties.maxMonomerLength = value;
      });

    this.root.appendChild(ui.form([maxMonomerLengthInput,]));
  }
}
