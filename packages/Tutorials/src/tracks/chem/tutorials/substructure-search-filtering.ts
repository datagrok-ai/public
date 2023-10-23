import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable, merge} from 'rxjs';
import $ from 'cash-dom';


export class SubstructureSearchFilteringTutorial extends Tutorial {
  get name() {
    return 'Substructure Search and Filtering';
  }

  get description() {
    return 'Datagrok offers an intuitive filtering functionality to explore and filter datasets. ' +
    'Specifically for molecules, Datagrok uses integrated sketchers to allow structure-based filtering. ' +
    'After applying the structure-based filter, queried substructures are highlighted in the filtered subset.';
  }

  get steps() { return 10; }
  
  helpUrl: string = '';

  protected async _run() {
    this.header.textContent = this.name;

    this.describe(this.description);

    const df = await grok.data.files.openTable('System:AppData/Tutorials/demo_smiles.csv');
    const tv = grok.shell.addTableView(df);

    const d = await this.openDialog('Initiate substructure search', '',
    this.getMenuItem('Chem'), `When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:
    <ul>
    <li><b>Chem</b> menu (it houses all chemical tools).</li>
    <li>Context menu (right-click for access)</li>
    <li><b>Context Panel</b> on the right.</li>
    </ul><br>
    <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">Learn more about exploring chemical data</a><br>
    Let’s use the <b>Chem</b> menu. Click it. Then, select <b>Search</b> > <b>Substructure Search</b>.`);

    const v = [...tv.viewers].find((v) => v.type === DG.VIEWER.FILTERS);
    await this.action('Specify substructure using sketcher', 
    d.onClose, undefined, 'In the sketcher, draw naphthalene.' +
    '<div class="ui-image" style="background-image: url(&quot;https://public.datagrok.ai/api/packages/published/' +
    'files/Tutorials/amuzychyna/b4MdalZhTYWUNZKAMDaN4uevTqxgauyU/6/images/naphtalene.png&quot;); height: 90px;"></div>' +
    `Note that as you draw, the chemical spreadsheet (grid) dynamically updates to show only the molecules with the
    specified substructure, with the substructure highlighted in each molecule.<br>Click <b>OK</b>.`);

    await this.buttonClickAction(v!.root, 'In the filter panel, click CLEAR', 'Clear', 
    'Remove the substructure filter.');

    await this.contextMenuAction('Use molecule to filter by substructure', 'Use as filter', undefined,
    'In the grid, right-click any molecule and select <b>Current Value</b> > <b>Use as flter</b>.');

    let d_: DG.Dialog;
    await this.action('In the filter panel, click the molecule', grok.events.onDialogShown.pipe(filter((dialog: DG.Dialog) => {
      if (dialog.title === '') {
        d_ = dialog;
        return true;
      }
      return false;
    })), v!.root.querySelector('.chem-canvas') as HTMLElement);

    await this.action('Modify it in the sketcher', d_!.onClose,
      undefined, 'In the sketcher, remove fragments<br>Click <b>OK</b>');

    await this.action('Adjust the filter setting to exclude the specified substructure', new Observable((subscriber: any) => {
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
    }), undefined, `You can switch between different filtering
    modes to include or exclude molecules with the specified substructure, show similar structures, and so on.<br>
    For this tutorial, let’s change the filter setting to exclude the specified substructure from the view:<br>
    On the <b>Filter Panel</b>, hover over a molecule and click the <b>Gear</b> icon.
    In the dropdown, select <b>Not contains</b>.`);

    await this.action('Toggle the filter', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          if (m.attributeName === 'class') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe(v!.root.querySelector('.d4-flex-col.ui-div.chem-filter.d4-filter-element')!, {attributes: true});
    }), undefined, `You can toggle a filter using a checkbox next to the filter’s name. Turn it off.`);

    const exp = ['NOCount', 'NumRotatableBonds', 'smiles'];
    await this.action('Add additional filters', v?.onEvent('d4-filter-added').pipe(filter(() => {
      return v.getOptions().look.filters.every((f: any) => exp.includes(f.column));
    }))!, $('.grok-icon.grok-font-icon-menu').get(0),
     `In the top left corner of the <b>Filter Panel</b>, click the <b>Hamburger</b> icon and select <b>Select columns...</b>.<br>
      In the dialog, select the <b>NOCount</b> column. Search for the <b>NumRotatableBonds</b> column and select it as well.<br>
      Click <b>OK</b>`);

    await this.action('Explore the dataset using newly added filters', merge(df.onSelectionChanged, df.onFilterChanged),
      undefined, `Hover over categories or distributions in the <b>Filter Panel</b> to instantly
      highlight relevant data points across all viewers.<br>
      Click any bin to select a corresponding group of rows in the grid.<br>
      Drag range controls to filter.<br>
      To learn more about filters, complete our <b>Filters</b> tutorial.`);
  }
}
