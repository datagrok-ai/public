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
    return 'You can interactively explore datasets using filters. ' +
    'Specifically for molecules, Datagrok uses integrated sketchers to filter by substructure.';
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
    Let’s use the <b>Chem</b> menu. Click it. Then, select <b>Search</b> > <b>Substructure Search</b>.`);

    const v = [...tv.viewers].find((v) => v.type === DG.VIEWER.FILTERS);
    await this.action('Specify substructure using sketcher', 
    d.onClose, undefined, 'In the sketcher, draw naphthalene.' +
    '<div class="ui-image" style="background-image: url(&quot;https://public.datagrok.ai/api/packages/published/' +
    'files/Tutorials/amuzychyna/b4MdalZhTYWUNZKAMDaN4uevTqxgauyU/6/images/naphtalene.png&quot;); height: 90px;"></div>' +
    `Note that as you draw, the chemical spreadsheet (grid) dynamically updates to show only the molecules with the
    specified substructure, highlighting it in each molecule.<br>Click <b>OK</b>.`);

    await this.buttonClickAction(v!.root, 'Remove the substructure filter', 'Clear', 
    'Use the <b>Filter Panel</b> on the left to clear the substructure filter.');

    await this.contextMenuAction('Use current molecule to filter by substructure', 'Use as filter', undefined,
    'In the grid, right-click any molecule and select <b>Current Value</b> > <b>Use as flter</b>.');

    let d_: DG.Dialog;
    await this.action('On the Filter Panel, click the molecule', grok.events.onDialogShown.pipe(filter((dialog: DG.Dialog) => {
      if (dialog.title === '') {
        d_ = dialog;
        return true;
      }
      return false;
    })), v!.root.querySelector('.chem-canvas') as HTMLElement);

    await this.action('Modify it in the sketcher', d_!.onClose,
      undefined, 'When finished, <br>click <b>OK</b>');

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
    }), v?.root.querySelector('.chem-canvas') as HTMLElement, `You can choose different filtering
    modes to include or exclude molecules with the specified substructure, show similar structures, and so on.<br>
    Let’s exclude the specified substructure from the view.<br>
    On the <b>Filter Panel</b>, hover over a molecule and click the <b>Gear</b> icon.
    From the dropdown, select <b>Not contains</b>.`);

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
    }), v?.root.querySelector('.d4-flex-row.d4-filter-header') as HTMLElement,
    `You can toggle a filter using a checkbox to the right of the filter’s name ("smiles"). Turn it off.`);

    const exp = ['NOCount', 'NumRotatableBonds', 'smiles'];
    await this.action('Add more filters', v?.onEvent('d4-filter-added').pipe(filter(() => {
      const filt = v.getOptions().look.filters;
      return filt.length === 3 && filt.every((f: any) => exp.includes(f.column));
    }))!, $('.panel-titlebar.panel-titlebar-tabhost .grok-icon.grok-font-icon-menu').get(0),
     `In the top left corner of the <b>Filter Panel</b>, click the <b>Hamburger</b> icon and choose <b>Select columns...</b>.<br>
      Then, in the dialog, select the <b>NOCount</b> column. Search for the <b>NumRotatableBonds</b> column and select it too.<br>
      Click <b>OK</b>`);

    await this.action('Explore the dataset using newly added filters', merge(df.onSelectionChanged, df.onFilterChanged),
      undefined, `Hover over categories or distributions in the <b>Filter Panel</b> to instantly
      highlight relevant data points.<br>
      Click any bin to select a corresponding group of rows in the grid.<br>
      Drag the range controls to filter.<br>
      To learn more about filters, complete our <b>Filters</b> tutorial.`);
  }
}
