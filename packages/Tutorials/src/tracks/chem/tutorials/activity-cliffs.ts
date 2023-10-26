import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable} from 'rxjs';
import $ from 'cash-dom';


export class ActivityCliffsTutorial extends Tutorial {
  get name() {
    return 'Activity Cliffs';
  }

  get description() {
    return 'The Activity Cliffs tool detects and visualizes pairs of molecules ' +
    'with highly similar structures but significantly different activity levels. ' +
    'It uses distance-based dimensionality reduction algorithms (such as tSNE and UMAP) ' +
    'to convert cross-similarities into a 2D scatterplot.';
  }

  get steps() {return 9;}
  
  helpUrl: string = '';
  demoTable: string = '';

  protected async _run() {
    this.header.textContent = this.name;

    this.describe(`The <b>Activity Cliffs</b> tool detects and visualizes pairs of molecules
      with highly similar structures but significantly different activity levels.
      It uses distance-based dimensionality reduction algorithms (such as tSNE and UMAP)
      to convert cross-similarities into a 2D scatterplot.<hr>
      <h3>Start the Activity Cliffs tool</h3>
      When you open a chemical dataset, Datagrok automatically detects molecules
      and shows molecule-specific tools, actions, and information. Access them through:
      <ul>
      <li><b>Chem</b> menu (it houses all chemical tools).</li>
      <li>Context menu (right-click for access)</li>
      <li><b>Context Panel</b> on the right.</li>
      </ul><br>
      Let's launch the Activity Cliffs tool.`);

    this.t = await grok.data.files.openTable('System:AppData/Tutorials/activity_cliffs.csv');
    const tv = grok.shell.addTableView(this.t);

    const d = await this.openDialog('Click Chem > Analyze > Activity Cliffs...', 'Activity Cliffs',
      this.getMenuItem('Chem'), ``);

    // <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">
    //Learn more about exploring chemical data</a><br>

    this.describe(`<h3>Set Parameters</h3>
    In the <b>Activity Cliffs</b> dialog, you can specify parameters like the similarity
    cutoff or the dimensionality reduction algorithm. For this tutorial, let's continue with the default settings.`);

    await this.action(' Click OK', d.onClose, $(d.root).find('button.ui-btn.ui-btn-ok')[0]);

    let v: DG.ScatterPlotViewer;
    await this.action('Wait for analysis to complete',
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
        const found = data.args.viewer.type === DG.VIEWER.SCATTER_PLOT;
        if (found)
          v = data.args.viewer;
        return found;
      })));

    this.describe(`<h3>Start analyzing the results</h3>
    Activity cliffs are visualized on an interactive scatterplot,
    where the proximity of the points indicates structural similarity.`);

    await this.action('Hover over data points for molecule information', new Observable((subscriber: any) => {
      const onMousemove = () => {
        console.log('UHFIUSDHFIUSDFUD');
        if ($('.d4-tooltip').css('display') === 'block') {
          subscriber.next(true);
          v.root.removeEventListener('mousemove', onMousemove);
        }
      };
      v.root.addEventListener('mousemove', onMousemove);
    }));

    await this.action('To view only the cliffs, toggle the Show only cliffs control.', new Observable((subscriber: any) => {
      $('.ui-input-switch').one('click', () => subscriber.next(true));
    }), $('.ui-input-switch').get(0));

    await this.action('Zoom in on the area of interest', v!.onZoomed, undefined,
      `On the scatterplot, the marker color corresponds to the activity level, and the size represents
    the maximum detected activity cliff for that molecule. The pairs with larger red markers may be
    particularly interesting as they indicate molecules with high activity levels and significant detected activity cliffs.<br>
    Let's zoom in. Use <b>Alt</b>+<b>Mouse Drag</b>.`);
    
    // await this.action('Step 5. [Explore the pairs of molecules]');
    
    //`The opacity of the green line connecting molecule pairs corresponds to the size of the activity cliff. 
    // Hover over the line to see the pair. Note the difference in their structures. 
    // Now, click the line. The <b>Context Panel</b> on the right updates to show the molecule pair.  
    //Click on the bottom molecule. This molecule becomes the current molecule in the grid.`

    await this.buttonClickAction(v!.root, 'Add a summary table with cliffs', '15 cliffs', 
      `To view all pairs of molecules identified as cliffs simultaneously, click <b>15 CLIFFS</b> at
    the top right corner of the scatterplot. This opens a table.`);

    let initH: number;
    const grid: DG.Grid = [...tv.viewers].find((v) => v.dataFrame.columns.length === 6) as DG.Grid;
    await this.action('Drag the top border of that table upwards to create more space for it.',
      new Observable((subscriber: any) => {
        const ro = new ResizeObserver((entries) => {
          if (initH === undefined) {
            initH = entries[0].contentRect.height;
            return;
          }
          if (entries[0].contentRect.height > initH) {
            subscriber.next(true);
            ro.unobserve(grid.root);
          }
        });
        ro.observe(grid.root);
      }), grid.root);

    await this.action('Now, click the first row.',
      new Observable((subscriber: any) => {
        const sub = grid.onCellClick.subscribe((cell) => {
          if (cell.gridRow === 0) {
            subscriber.next(true);
            sub.unsubscribe();
          }
        });
      }), undefined, `Note the changes in the scatterplot's view and the source table.`);
  }
}
