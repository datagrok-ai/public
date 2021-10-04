import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from "../../../tutorial";


export class DescriptorsTutorial extends Tutorial {
  get name() {
    return 'Molecular Descriptors';
  }

  get description(): string {
    return 'Chemical structures are analyzed via molecular descriptors. ' +
      'These are numerical values that characterize properties of molecules.'+
      '<a href="https://datagrok.ai/help/domains/chem/descriptors" target="_blank" class="ui-link d4-link-external">Read more about molecular descriptors.</a>';
  }

  get steps() { return 3; }

  demoTable: string = 'chem/smiles_only.csv';

  protected async _run(): Promise<void> {
    this.header.textContent = this.name;
    this.describe('Chemical structures are analyzed via molecular descriptors. ' +
      'These are numerical values that characterize properties of molecules. ' +
      'In the next steps, we will learn how to calculate descriptors in Datagrok.');

    let pp: DG.Accordion;
    const ppDescription = 'In the property panel, you should now see all the actions ' +
      'applicable to the column under the "Actions" section. Let\'s calculate some ' +
      'molecular descriptors for the given dataset.';

    await this.action('Click on the header of a column with molecules',
      grok.events.onAccordionConstructed.pipe(filter((acc) => {
        if (acc.context instanceof DG.Column && acc.context?.name == 'canonical_smiles') {
          pp = acc;
          return true;
        }
        return false;
      })), null, ppDescription);

    pp!.getPane('Actions').expanded = true;
    const descriptorsLabel = $(pp!.root)
      .find('label.d4-link-action')
      .filter((idx, el) => el.textContent == 'Chem | Descriptors...')[0];

    const dlgDescription = 'The platform supports various groups of descriptors. You can select ' +
      'them either by category or individually. Let\'s compute the average molecular weight and ' +
      'Lipinski parameters. These are common descriptors used for drug-likeness prediction.';
    const dlg = await this.openDialog('Click on "Chem | Descriptors..."', 'Descriptors',
      descriptorsLabel, dlgDescription);

    const descriptors = ['MolWt', 'FractionCSP3', 'HeavyAtomCount', 'NHOHCount',
      'NOCount', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAliphaticRings',
      'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors',
      'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
      'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount'];

    const groupHints = $(dlg.root)
      .find('.d4-tree-view-node')
      .filter((idx, el) => el.textContent === 'Descriptors' || el.textContent === 'Lipinski')
      .get();

    const groupDescription = 'To select the "<b>MolWt</b>" expand the "<b>Descriptors</b>" group.\nThe ' +
      'results of computation will be added to the main dataframe. Take some time to examine the findings.';

    await this.action('Select "MolWt" and "Lipinski" group and press the "OK" button',
      grok.functions.onAfterRunAction.pipe(filter((call) => {
        const inputs = call.inputs.get('descriptors');
        return call.func.name === 'ChemDescriptors' &&
          descriptors.length == inputs.length &&
          descriptors.every((d) => inputs.includes(d));
      })), groupHints, groupDescription);

  }
}
