import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable} from 'rxjs';
import $ from 'cash-dom';


export class SimilarityDiversitySearchTutorial extends Tutorial {
  get name() {
    return 'Similarity and Diversity Search';
  }

  get description() {
    return 'Datagrok offers two analytical tools to help you analyze a collection ' +
    'of molecules based on molecular similarity: Similarity Search and Diversity ' +
    'Search. Similarity Search finds structures similar to the reference molecule, ' +
    'while Diversity Search shows N molecules of different chemical classes in the dataset.';
  }

  get steps() {return 12;}
  
  helpUrl: string = '';
  demoTable: string = '';
  prerequisites: TutorialPrerequisites = {packages: ['Chem']};
  // manualMode = true;

  protected async _run() {
    this.header.textContent = this.name;

    this.describe('Datagrok offers two analytical tools to help you analyze a collection ' +
    'of molecules based on molecular similarity: <b>Similarity Search</b> and <b>Diversity ' +
    'Search</b>. Similarity Search finds structures similar to the reference molecule, ' +
    'while Diversity Search shows N molecules of different chemical classes in the dataset.<hr>');

    this.t = await grok.data.files.openTable('System:AppData/Tutorials/demo_smiles.csv');
    grok.shell.addTableView(this.t);

    this.describe(`<h3>Start the similarity and diversity search</h3>
    When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:<br>
    <ul>
    <li><b>Chem</b> menu (it houses all chemical tools).</li>
    <li>Context menu (right-click for access)</li>
    <li><b>Context Panel</b> on the right.</li>
    </ul><br>
    Let’s search for similar and diverse structures. On the <b>Top Menu</b>, click <b>Chem</b> > <b>Search</b>
    > <b>Acitvity Cliffs...</b> > <b>Similarity Search...</b> Repeat to initiate the diversity search.`, true);

    let sim: DG.Viewer;
    let div: DG.Viewer;
    await this.action('Click Chem > Search > Similarity/Diversity Search...', grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
      if (data.args.viewer.type === 'Chem Similarity Search')
        sim = data.args.viewer;
      if (data.args.viewer.type === 'Chem Diversity Search')
        div = data.args.viewer;
      return Boolean(sim && div);
    })), this.getMenuItem('Chem'));

    this.describe(`<h3>Explore the dataset using similarity and diversity viewers</h3>
    The similarity and diversity viewers are interactive and are synchronized with each other,
    chemical spreadsheet (grid), and other viewers.<br>On the <b>Most similar structures</b> viewer, click the
    molecule next to the reference molecule<br>
    Note how it’s now the current molecule in the grid.<br>
    Now, click any molecule in the diversity viewer.<br>
    Note the change in the similarity viewer.`, true);

    let i = 0;
    await this.action('Select different molecules', this.t.onCurrentRowChanged.pipe(filter(() => {
      i++;
      return i >= 2;
    })));

    this.describe(`<h3>Lock in a reference molecule</h3>
    By default, a reference molecule follows the current row. However, you can lock it in.`, true);

    await this.contextMenuAction('In the similarity viewer, right-click a reference molecule and select Properties...', 'Properties...',
      $('.d4-chem-similarity-search .ui-div.d4-flex-col.d4-current').get(0));

    await this.action('Under Misc, clear the Follow Current Row checkbox', new Observable((subscriber: any) => {
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

    this.describe(`<h3>Specify a custom reference molecule</h3>
    Another option to set a reference molecule is to use a sketcher. Click the <b>Edit</b> icon on the reference molecule card.`, true);

    const d = await this.openDialog('In the Most similar structures viewer > reference molecule tile, click the edit (pencil) icon', '',
      $('.grok-icon.fal.fa-pen.similarity-search-edit.chem-mol-view-icon').get(0));

    await this.action('Set new reference molecule', 
      d.onClose, undefined, `For this tutorial, let's paste this identifier in the sketcher and press
      <b>Enter</b>: <br><b>RCKUZCZZUOYUBE-DHDCSXOGSA-N</b><br> 
      The sketcher updates to show the new structure.<br>
      Click <b>OK</b> to apply.<br>
      If you want, you can also draw the structure manually.`);

    this.describe(`<h3>Get insights using Context Panel</h3>
      In the similarity viewer, hover over the reference molecule and then click the <b>More</b> icon. Then, click Explore.
      The <b>Context Panel</b> on the right now shows molecule-specific information.<br>
      Under <b>Biology</b>, click <b>Toxicity</b>.<br>
      Now click <b>Drug Likeness</b>.<br>
      On the <b>Drug Likeness</b> info pane, scroll through the substructure list.<br>
      These info panes are interactive and fully customizable.<br>
      <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">
      Learn about chemical info panes here</a>`, true);
    
    await this.contextMenuAction('Click the More icon > Explore', 'Explore',
      $('.d4-chem-similarity-search .grok-icon.fal.fa-ellipsis-v.chem-mol-view-icon.pep-more-icon').get(0));
    
    let k = 0;
    await this.action('Explore Context Panel', new Observable((subscriber: any) => {
      const func = () => {
        console.log(k);
        k++;
        if (k >= 3) {
          document.querySelector('.grok-prop-panel')?.removeEventListener('click', func);
          subscriber.next(true);
        }
      };
      $('.grok-prop-panel').on('click', func);
    }));

    this.describe(`<h3>Add more data</h3>
    You can change search parameters and other settings anytime. To do so, use the viewer’s settings.
    For now, let’s add some data on the molecule cards:`, true);

    await this.action('In the Most diverse structures viewer, click the Tanimoto, Morgan link', new Observable((subscriber: any) => {
      $('.d4-chem-diversity-search .ui-link').on('click', () => subscriber.next(true));
    }), $('.d4-chem-diversity-search .ui-link').get(0));

    await this.action('Select NumValenceElectrons column', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          //@ts-ignore
          if (m.target.innerText === '1 / 31') {
            console.log(m);
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {subtree: true, attributes: true, childList: true, characterData: true});
    }), undefined, `<ol>
    <li>Next to <b>Molecule Properties</b>, click the column selector icon ('...').</li>
    <li>Select these column: NumValenceElectrons.</li>
    <li>Click <b>OK</b>.</li></ol>`);

    this.describe(`<h3>Color-code for quick profiling</h3>
    Datagrok viewers can pick up color-coding from the grid. Let’s color-code the newly added values:
    <ol>
    <li>In the grid, locate the <b>NumValenceElectrons</b> column and right-click its header.</li>
    <li>Select <b>Color Coding</b> > <b>Linear</b>.</li>
    </ol>
    You can further adjust how the color setting in the viewer’s properties.`, true);

    await this.action('Add Color Coding for NumValenceElectrons column', this.t.onMetadataChanged.pipe(filter((data: DG.EventData) => {
      return data.args.change === 'set' && data.args.key === '.color-coding-type' && data.args.value === 'Linear' &&
        data.args.source.name === 'NumValenceElectrons';
    })));
  }
}
