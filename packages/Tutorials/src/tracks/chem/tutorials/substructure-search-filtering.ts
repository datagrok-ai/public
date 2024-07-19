import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable, combineLatest} from 'rxjs';
import $ from 'cash-dom';
import { _package } from '../../../package';


export class SubstructureSearchFilteringTutorial extends Tutorial {
  get name() {
    return 'Substructure Search and Filtering';
  }

  get description() {
    return 'You can interactively explore datasets using filters. ' +
    'Specifically for molecules, Datagrok uses integrated sketchers to filter by substructure.';
  }

  get steps() {return 10;}
  
  helpUrl: string = 'https://datagrok.ai/help/datagrok/solutions/domains/chem/#substructure-search--filtering';
  demoTable: string = '';
  prerequisites: TutorialPrerequisites = {packages: ['Chem']};
  // manualMode = true;

  protected async _run() {
    this.header.textContent = this.name;

    this.describe(this.description + '<hr>');

    this.t = await grok.data.loadTable(`${_package.webRoot}files/demo_smiles.csv`);
    const tv = grok.shell.addTableView(this.t);

    this.title('Initiate substructure search', true);
    this.describe(`When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:<br>
    <ul>
    <li><b>Chem</b> menu (it contains all chemical tools)</li>
    </ul><br>
    Let’s use the <b>Chem</b> menu.`);

    const d = await this.openDialog('Click Chem > Search > Substructure Search…', '', this.getMenuItem('Chem', true));

    this.title('Specify substructure using sketcher', true);
    this.describe('In the sketcher, draw naphthalene.' +
    '<div class="ui-image" id="naphthalene" style="background-image: url(&quot;https://public.datagrok.ai/api/packages/published/' +
    'files/Tutorials/amuzychyna/b4MdalZhTYWUNZKAMDaN4uevTqxgauyU/6/images/naphtalene.png&quot;); height: 90px;"></div>' +
    `Note that as you draw, the chemical spreadsheet (grid) dynamically updates to show only the molecules with the
    specified substructure, highlighting it in each molecule.<br>Click <b>OK</b>.`);

    const v = [...tv.viewers].find((v) => v.type === DG.VIEWER.FILTERS);
    await this.action('Change substructure', d.onClose);

    this.title('Remove the substructure filter', true);
    this.describe(`Use the <b>Filter Panel</b> on the left to clear the substructure filter`);

    await this.buttonClickAction(v!.root, 'In the filter panel, click CLEAR', 'Clear');

    this.title('Use current molecule to filter by substructure', true);

    await this.contextMenuAction('In the grid, right-click any molecule and select Current Value > Use as flter', 'Use as filter');

    let d_: DG.Dialog;
    await this.action('On the Filter Panel, click the molecule and modify it in the sketcher',
    grok.events.onDialogShown.pipe(filter((dialog: DG.Dialog) => {
      if (dialog.title === '') {
        d_ = dialog;
        return true;
      }
      return false;
    })), v!.root.querySelector('.chem-canvas') as HTMLElement);

    await this.action('Click OK', d_!.onClose);

    this.title('Adjust the filter setting to exclude the specified substructure', true);
    this.describe(`You can choose different filtering modes to include or exclude molecules with the
    specified substructure, show similar structures, and so on.<br>
    Let’s exclude the specified substructure from the view.`);

    await this.action('Exclude the specified substructure from the view', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          //@ts-ignore
          if (m.target.value === 'Not contains') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe(v!.root.querySelector('.d4-flex-col.d4-filter')!, {subtree: true, attributes: true});
    }), v?.root.querySelector('.d4-flex-col.d4-filter') as HTMLElement,
    `<ol><li>On the <b>Filter Panel</b>, hover over a molecule and click the <b>Gear</b> icon.</li>
    <li>From the dropdown, select <b>Not contains<b>.</li><ol>`);

    this.title('Toggle the filter', true);
    this.describe(`You can toggle a filter using a checkbox to the right of the filter's name ("smiles"). Turn it off.`);

    await this.action('In the filter panel, turn off the filter by clearing the checkbox', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          if (m.attributeName === 'class') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe(v!.root.querySelector('.d4-flex-col.ui-div.chem-filter.d4-filter-element')!, {attributes: true});
    }), v?.root.querySelector('.d4-flex-col.d4-filter') as HTMLElement);

    this.title('Add more filters', true);
    this.describe(`In the top left corner of the <b>Filter Panel</b>, click the <b>Hamburger</b> icon and choose <b>Select columns...</b>.<br>
    Then, in the dialog, select the <b>NOCount</b> column. Search for the <b>NumRotatableBonds</b> column and select it too.<br>
    Click <b>OK</b>`);

    const exp = ['NOCount', 'NumRotatableBonds', 'smiles'];
    await this.action('Select columns to be used as filers', v?.onEvent('d4-filter-added').pipe(filter(() => {
      const filt = v.getOptions().look.filters;
      return filt.length === 3 && filt.every((f: any) => exp.includes(f.column));
    }))!, $('.panel-titlebar.panel-titlebar-tabhost .grok-icon.grok-font-icon-menu').get(0));

    this.title('Explore the dataset using newly added filters', true);
    this.describe(`Hover over categories or distributions in the <b>Filter Panel</b> to instantly
    highlight relevant data points.<br>
    Click any bin to select a corresponding group of rows in the grid.<br>
    Drag the range controls to filter.<br>
    To learn more about filters, complete our <b>Filters</b> tutorial.`);

    await this.action('Interact with filters by changing their values. After that select rows of your interest', combineLatest([this.t.onSelectionChanged, this.t.onFilterChanged]));
  }
}
