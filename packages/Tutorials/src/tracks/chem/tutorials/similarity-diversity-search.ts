import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable} from 'rxjs';
import $ from 'cash-dom';


export class SimilarityDiversitySearchTutorial extends Tutorial {
  get name() {
    return 'Similarity and Diversity Search';
  }

  get description() {
    return 'Datagrok offers two analytical tools to help you analyze a collection ' +
    'of molecules based on molecular similarity: <b>Similarity Search</b> and <b>Diversity ' +
    'Search</b>. Similarity Search finds structures similar to the reference molecule, ' +
    'while Diversity Search shows N molecules of different chemical classes in the dataset.';
  }

  get steps() {return 10;}
  
  helpUrl: string = '';
  demoTable: string = '';

  protected async _run() {
    this.header.textContent = this.name;

    this.describe(this.description);

    this.t = await grok.data.files.openTable('System:AppData/Tutorials/demo_smiles.csv');
    grok.shell.addTableView(this.t);

    let sim: DG.Viewer;
    let div: DG.Viewer;
    await this.action('Start the similarity and diversity search', grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
      if (data.args.viewer.type === 'Chem Similarity Search')
        sim = data.args.viewer;
      if (data.args.viewer.type === 'Chem Diversity Search')
        div = data.args.viewer;
      return Boolean(sim && div);
    })), this.getMenuItem('Chem'), `When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:
    <ul>
    <li><b>Chem</b> menu (it houses all chemical tools).</li>
    <li>Context menu (right-click for access)</li>
    <li><b>Context Panel</b> on the right.</li>
    </ul><br>
    Let’s search for similar and diverse structures.<br>
    Click <b>Chem</b> > <b>Search</b> > <b>Similarity Search...</b><br>
    Click <b>Chem</b> > <b>Search</b> > <b>Diversity Search...</b>`);

    let i = 0;
    await this.action('Explore the dataset using similarity and diversity viewers', this.t.onCurrentRowChanged.pipe(filter(() => {
      i++;
      return i >= 2;
    })), undefined, `The similarity and diversity viewers are interactive and are synchronized with each other,
    chemical spreadsheet (grid), and other viewers.<br>On the <b>Most similar structures</b> viewer, click the
    molecule next to the reference molecule<br>
    Note how it’s now the current molecule in the grid.<br>
    Now, click any molecule in the diversity viewer.<br>
    Note the change in the similarity viewer.`);

    const d = await this.openDialog('Specify a reference molecule', '',
      $('.grok-icon.fal.fa-pen.similarity-search-edit.chem-mol-view-icon').get(0));

    await this.action('Specify a reference molecule', 
      d.onClose, undefined, `To specify a reference molecule, click the <b>Edit</b> icon on the reference molecule card.<br>
    For this tutorial, let's paste this identifier in the sketcher and press <b>Enter</b>: <br><b>RCKUZCZZUOYUBE-DHDCSXOGSA-N</b><br> 
    The sketcher updates to show the new structure.<br>
    Click <b>OK</b> to apply.<br>
    If you want, you can also draw the structure manually.`);

    await this.contextMenuAction('Lock in a reference molecule', 'Properties...',
      $('.d4-chem-similarity-search .ui-div.d4-flex-col.d4-current').get(0),
      `By default, a reference molecule follows the current row. However, you can lock it in.<br>
      In the similarity viewer, right-click a reference molecule and select <b>Properties...</b>`);

    await this.action('Under <b>Misc</b>, clear the <b>Follow Current Row</b> checkbox', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          //@ts-ignore
          if (m.attributeName === 'class' && m.target.innerText === 'Follow Current Row') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {subtree: true, attributes: true});
    }));

    let j = 0;
    await this.action('Click anywhere in the viewer and in the grid', this.t.onCurrentRowChanged.pipe(filter(() => {
      j++;
      return j >= 2;
    })), undefined, `Click anywhere in the viewer. Now click any molecule in the grid.
    The current molecule changes, but the reference molecule stays the same.`);

    await this.contextMenuAction('Get insights using <b>Context Panel</b>', 'Explore',
      $('.d4-chem-similarity-search .grok-icon.fal.fa-ellipsis-v.chem-mol-view-icon.pep-more-icon').get(0),
      `In the similarity viewer, hover over the reference molecule and then click the <b>More</b> icon. Then, click <b>Explore</b>.
       The <b>Context Panel</b> on the right now shows molecule-specific information.`);
    
    // await this.action('Get insights using Context Panel', new Observable((subscriber: any) => {
    //   console.log($('.grok-prop-panel .d4-flex-wrap.ui-div'));
    //   $('.grok-prop-panel .d4-flex-wrap.ui-div').on('scroll', () => {console.log('SCROLL'); subscriber.next(true)});
    // }), undefined, `Under <b>Biology</b>, click <b>Toxicity</b>.<br>
    // Now click <b>Drug Likeness</b>.<br>
    // On the <b>Drug Likeness</b> info pane, scroll through the substructure list.<br>
    // These info panes are interactive and fully customizable.<br>
    // <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">
    //Learn about chemical info panes here</a>`);

    await this.action('Access the viewer settings', new Observable((subscriber: any) => {
      $('.d4-chem-diversity-search .ui-link').on('click', () => subscriber.next(true));
    }), $('.d4-chem-diversity-search .ui-link').get(0), `You can change search parameters and other settings anytime.
    To do so, use the viewer’s settings. In the top right corner of the diversity viewer, click <b>Tanimoto, Morgan</b><br>`);

    await this.action('Add more data', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          //@ts-ignore
          if (m.target.innerText === '1 / 31') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {subtree: true, attributes: true, childList: true, characterData: true});
    }), undefined, `Let’s add some data on the molecule card:<br>
      Under <b>Misc</b>, next to <b>Molecule Properties</b>, click the column selector icon.<br>
      In the dialog, select the <b>NumValenceElectrons</b> column.<br>
      Click <b>OK</b>.`);

    await this.action('Color-code for quick profiling', this.t.onMetadataChanged.pipe(filter((data: DG.EventData) => {
      return data.args.change === 'set' && data.args.key === '.color-coding-type' && data.args.value === 'Linear' &&
        data.args.source.name === 'NumValenceElectrons';
    })));
  }
}
