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

    await this.action('Explore the dataset using similarity and diversity viewers', merge(new Observable((subscriber: any) => {
      sim.root.querySelector('.ui-div.d4-flex-col:nth-child(2)')
        ?.addEventListener('click', () => subscriber.next(true));
    }), new Observable((subscriber: any) => {
      div.root.querySelector('.d4-flex-wrap.ui-div')
        ?.addEventListener('click', () => subscriber.next(true));
    })), undefined, `The similarity and diversity viewers are interactive and are synchronized with each other,
    chemical spreadsheet (grid), and other viewers.<br>In the <b>Most similar structures</b> viewer, click a tile next to the reference molecule<br>
    In the diversity viewer, click any tile`);
  }
}
