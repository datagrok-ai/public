import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable, fromEvent} from 'rxjs';
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

  get steps() {return 12;}
  
  helpUrl: string = '';
  demoTable: string = '';
  prerequisites: TutorialPrerequisites = {packages: ['Chem']};
  // manualMode = true;

  protected async _run() {
    this.header.textContent = this.name;
    grok.shell.windows.showContextPanel = true;

    this.describe(`The <b>Activity Cliffs</b> tool detects and visualizes pairs of molecules
      with highly similar structures but significantly different activity levels.
      It uses distance-based dimensionality reduction algorithms (such as tSNE and UMAP)
      to convert cross-similarities into a 2D scatterplot.<hr>`);

    this.t = await grok.data.files.openTable('System:AppData/Tutorials/activity_cliffs.csv');
    const tv = grok.shell.addTableView(this.t);

    this.describe(`<h3>Start the Activity Cliffs tool</h3>
    When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:<br>
    <ul>
    <li><b>Chem</b> menu (it houses all chemical tools).</li>
    <li>Context menu (right-click for access)</li>
    <li><b>Context Panel</b> on the right.</li>
    </ul><br>
    Let's launch the Activity Cliffs tool.`, true);

    const d = await this.openDialog('Click Chem > Analyze > Activity Cliffs...', 'Activity Cliffs',
      this.getMenuItem('Chem'), ``);

    // <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">
    //Learn more about exploring chemical data</a><br>

    this.describe(`<h3>Set Parameters</h3>
    In the <b>Activity Cliffs</b> dialog, you can specify parameters like the similarity
    cutoff or the dimensionality reduction algorithm. For this tutorial, let's continue with the default settings.`);

    await this.action(' Click OK', d.onClose, $(d.root).find('button.ui-btn.ui-btn-ok')[0], '', false);

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
    where the proximity of the points indicates structural similarity.`, true);

    await this.action('Hover over data points for molecule information', new Observable((subscriber: any) => {
      const onMousemove = () => {
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

    this.describe(`<h3>Zoom in on the area of interest</h3>
    On the scatterplot, the marker color corresponds to the activity level, and the size represents
    the maximum detected activity cliff for that molecule. The pairs with larger red markers may be
    particularly interesting as they indicate molecules with high activity levels and significant detected activity cliffs.`, true);

    await this.action('Let’s zoom in. Use Alt + Mouse Drag', v!.onZoomed);
    
    this.describe(`<h3>Explore the pairs of molecules</h3>
    The opacity of the green line connecting molecule pairs corresponds to the size of the activity cliff.`, true);

    await this.action('Hover over the line to see the pair',
      grok.events.onTooltipShown.pipe(filter((_) => {
        return $('.d4-tooltip').text().startsWith('smilesActivity');
    })), undefined, 'Note the difference in their structures');

    await this.action('Now, click the line', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          if (m.addedNodes.length && m.addedNodes[0].textContent?.startsWith('Activity cliffs')) {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {childList: true, subtree: true});
    }));

    await this.action('Click on the bottom molecule', this.t.onCurrentRowChanged,
      $('.grok-prop-panel').get(0), 'This molecule becomes the current molecule in the grid');

    this.describe(`<h3>Add a summary table with cliffs</h3>
    To view all pairs of molecules identified as cliffs simultaneously, click <b>15 CLIFFS</b>
    at the top right corner of the scatterplot. A table is added to the view.<br>
    Drag the top border of that table upwards to create more space for it.<br>
    Now, click the first row. Note the changes in the scatterplot’s view and the source table.`, true);

    await this.buttonClickAction(v!.root, 'Click the 15 CLIFFS link', '15 cliffs');

    let initH: number;
    const grid: DG.Grid = [...tv.viewers].find((v) => v.dataFrame.columns.length === 6) as DG.Grid;
    await this.action(`Drag the bottom table's border up`,
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

    await this.action('In the cliffs table, in the first row, click any cell in the first row',
      new Observable((subscriber: any) => {
        const sub = grid.onCellClick.subscribe((cell) => {
          if (cell.gridRow === 0) {
            subscriber.next(true);
            sub.unsubscribe();
          }
        });
    }));
  }
}
