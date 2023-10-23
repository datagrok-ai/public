import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import {fromEvent, interval, Observable, Subject, merge} from 'rxjs';
import $ from 'cash-dom';


export class SimilarityDiversitySearchTutorial extends Tutorial {
  get name() {
    return 'Similarity and Diversity Search';
  }

  get description() {
    return 'Datagrok offers two analytical tools to help you analyze a collection ' +
    'of molecules based on molecular similarity: Similarity Search and Diversity ' +
    'Search viewers. Similarity search finds structures similar to the reference molecule, ' +
    'while diversity search shows N molecules of different chemical classes in the dataset.';
  }

  get steps() { return 10; }
  
  helpUrl: string = '';

  protected async _run() {
    this.header.textContent = this.name;

    this.describe(this.description);

    const df = await grok.data.files.openTable('System:AppData/Tutorials/demo_smiles.csv');
    const tv = grok.shell.addTableView(df);

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
    <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">Learn more about exploring chemical data</a><br>
    Let’s search for similar and diverse structures.<br>
    Click <b>Chem</b> > <b>Search</b> > <b>Similarity Search...</b><br>
    Click <b>Chem</b> > <b>Search</b> > <b>Diversity Search...</b>`);

    let i = 0;
    await this.action('Explore the dataset using similarity and diversity viewers', df.onCurrentRowChanged.pipe(filter(() => {
      i++;
      return i >= 2;
    })), undefined, `The similarity and diversity viewers are interactive and are synchronized with each other,
    chemical spreadsheet (grid), and other viewers.<br>In the <b>Most similar structures</b> viewer, click a tile next to the reference molecule<br>
    In the diversity viewer, click any tile`);

    const d = await this.openDialog('Most similar structures viewer > Reference molecule tile, click the edit (pencil) icon.', '',
      $('.grok-icon.fal.fa-pen.similarity-search-edit.chem-mol-view-icon').get(0));

    await this.action('Specify a reference molecule', 
    d.onClose, undefined, `In the sketcher, paste the molecule’s identifier. For this tutorial, paste<br><b>RCKUZCZZUOYUBE-DHDCSXOGSA-N</b><br>and press <b>Enter</b>.
    The sketcher updates to show the new structure.<br>
    Click <b>OK</b> to apply.<br>
    If you want, you can also draw the structure manually.`);

    await this.contextMenuAction('Lock in a reference molecule', 'Properties...',
      $('.ui-div.d4-flex-col.d4-current').get(0),
      `By default, a reference molecule follows the current row. If you click a different molecule, the similarity viewer updates accordingly.
      You can lock in a specific reference molecule.<br>
      In the similarity viewer, right-click a reference molecule and select <b>Properties...</b>`);

    await this.action('Under Misc info pane, clear the Follow Current Row checkbox', new Observable((subscriber: any) => {
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
    await this.action('Click anywhere in the viewer and in the grid', df.onCurrentRowChanged.pipe(filter(() => {
      j++;
      return j >= 2;
    })), undefined, `Click anywhere in the viewer. Now click any molecule in the grid.
    The current molecule changes, but the reference molecule stays the same.`);

    await this.contextMenuAction('Click More > Explore', 'Explore',
    $('.grok-icon.fal.fa-ellipsis-v.chem-mol-view-icon.pep-more-icon').get(0),
      `In the similarity viewer, hover over the reference molecule and then click the <b>More</b> icon. Click <b>Explore</b>.
       The <b>Context Panel</b> on the right now shows molecule-specific information.`);
    
    // await this.action('Get insights using Context Panel', new Observable((subscriber: any) => {
    //   console.log($('.grok-prop-panel .d4-flex-wrap.ui-div'));
    //   $('.grok-prop-panel .d4-flex-wrap.ui-div').on('scroll', () => {console.log('SCROLL'); subscriber.next(true)});
    // }), undefined, `Under <b>Biology</b>, click <b>Toxicity</b>.<br>
    // Now click <b>Drug Likeness</b>.<br>
    // On the <b>Drug Likeness</b> info pane, scroll through the substructure list.<br>
    // <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">Learn about chemical info panes here</a>`);

    await this.action('In the top right corner of the diversity viewer, click Tanimoto, Morgan', new Observable((subscriber: any) => {
      $('.d4-chem-diversity-search .ui-link').on('click', () => subscriber.next(true));
    }), $('.d4-chem-diversity-search .ui-link').get(0), 'You can change search parameters and other settings anytime. To do so, use the viewer’s settings.');

    await this.action('Change search parameters and show additional data', new Observable((subscriber: any) => {
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
    }), undefined, `For this tutorial, let’s show additional data on molecule cards:<br>
      Under <b>Misc</b>, next to <b>Molecule Properties</b>, click the column selector icon.<br>
      In the dialog, select the <b>NumValenceElectrons</b> column.<br>
      Click <b>OK</b>.`);
  }
}
