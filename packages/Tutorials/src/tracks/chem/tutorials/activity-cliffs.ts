import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable} from 'rxjs';
import $ from 'cash-dom';
import { _package } from '../../../package';


export class ActivityCliffsTutorial extends Tutorial {
  get name() {
    return 'Activity Cliffs';
  }

  get description() {
    return 'The Activity Cliffs tool detects and visualizes pairs of molecules ' +
    'with highly similar structures but significantly different activity levels. ' +
    'It uses distance-based dimensionality reduction algorithms ' +
    'to convert cross-similarities into a 2D scatterplot.';
  }

  get steps() {return 12;}
  
  helpUrl: string = 'https://datagrok.ai/help/datagrok/solutions/domains/chem/#activity-cliffs';
  demoTable: string = '';
  prerequisites: TutorialPrerequisites = {packages: ['Chem']};
  // manualMode = true;

  protected async _run() {
    this.header.textContent = this.name;
    grok.shell.windows.showContextPanel = true;

    this.describe(`The <b>Activity Cliffs</b> tool detects and visualizes pairs of molecules
      with highly similar structures but significantly different activity levels.
      It uses distance-based dimensionality reduction algorithms 
      to convert cross-similarities into a 2D scatterplot.<hr>`);

    this.t = await grok.data.loadTable(`${_package.webRoot}files/activity_cliffs.csv`);
    const tv = grok.shell.addTableView(this.t);

    this.title('Start the Activity Cliffs tool', true);
    this.describe(`When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:<br>
    <ul>
    <li><b>Chem</b> menu (it contains all chemical tools)</li>
    </ul><br>
    Let's launch the Activity Cliffs tool.`);

    const d = await this.openDialog('On the Top Menu, click Chem > Analyze > Activity Cliffs...',
    'Activity Cliffs', this.getMenuItem('Chem', true));

    // <a href="https://datagrok.ai/help/datagrok/solutions/domains/chem/#exploring-chemical-data">
    //Learn more about exploring chemical data</a><br>

    this.title('Set Parameters', true);
    this.describe(`In the <b>Activity Cliffs</b> dialog, you can specify parameters like the similarity
    cutoff or the dimensionality reduction algorithm. For this tutorial, let's continue with the default settings.`);

    await this.action('Click OK', d.onClose, $(d.root).find('button.ui-btn.ui-btn-ok')[0]);

    let v: DG.ScatterPlotViewer;
    await this.action('Wait for analysis to complete',
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
        const found = data.args.viewer.type === DG.VIEWER.SCATTER_PLOT;
        if (found)
          v = data.args.viewer;
        return found;
      })));

    this.title('Start analyzing the results', true);
    this.describe(`Activity cliffs are visualized on an interactive scatterplot,
    where the proximity of the points indicates structural similarity.`);

    await this.action('Hover over data points for molecule information', new Observable((subscriber: any) => {
      const onMousemove = () => {
        if ($('.d4-tooltip').css('display') === 'block') {
          subscriber.next(true);
          v.root.removeEventListener('mousemove', onMousemove);
        }
      };
      v.root.addEventListener('mousemove', onMousemove);
    }));

    await this.action('To view only the cliffs, toggle Show only cliffs.', new Observable((subscriber: any) => {
      $('.ui-input-switch').one('click', () => subscriber.next(true));
    }), $('.ui-input-switch').get(0));

    this.title('Zoom in on the area of interest', true);
    this.describe(`On the scatterplot, the marker color corresponds to the activity level, and the size represents
    the maximum detected activity cliff for that molecule. The pairs with larger red markers may be
    particularly interesting as they indicate molecules with high activity levels and significant detected activity cliffs.<br>
    Let’s zoom in. Use <b>Alt + Mouse Drag</b>.`);

    await this.action('Press Use Alt + Mouse Drag to zoom in', v!.onZoomed);
    
    this.title('Explore the pairs of molecules', true);
    this.describe(`The opacity of the green line connecting molecules corresponds to the magnitude of the activity cliff.
    Hover over it to view structural differences between a molecule pair, or click it to see the pair in the <b>Context Panel</b>`);

    await this.action('Hover over the green line to see the pair of molecules',
      grok.events.onTooltipShown.pipe(filter((_) => {
        return $('.d4-tooltip').text().startsWith('smilesActivity');
    })), undefined, 'Note the difference in their structures');

    await this.action('Click on the green line connecting that molecule pair', new Observable((subscriber: any) => {
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

    await this.action('On the Context Panel, click any molecule', this.t.onCurrentRowChanged,
      undefined, 'Note changes in the grid');

    this.title('Add a summary table with cliffs', true);
    this.describe(`Let’s help our analysis and include a table listing all pairs identified as cliffs.
    In Datagrok, all viewers are synchronized. This means they respond to the same filtering, selection,
    and other interactions.`);

    await this.buttonClickAction(v!.root, 'At the top right corner of the scatterplot, click 15 CLIFFS', '15 cliffs',
    'A table is added to the view');

    let initH: number;
    const grid: DG.Grid = [...tv.viewers].find((v) => v.dataFrame.columns.length === 6) as DG.Grid;
    await this.action(`Drag the top border of that table upwards to create more space for it`,
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

    await this.action('In the cliffs table, click any cell in the first row',
      new Observable((subscriber: any) => {
        const sub = grid.onCellClick.subscribe((cell) => {
          if (cell.gridRow === 0) {
            subscriber.next(true);
            sub.unsubscribe();
          }
        });
    }), undefined, 'Note the changes on the scatterplot’s view and the source table');
  }
}
