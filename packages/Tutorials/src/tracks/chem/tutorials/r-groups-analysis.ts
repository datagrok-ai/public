import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable, combineLatest, interval} from 'rxjs';
import $, {Cash} from 'cash-dom';
import { _package } from '../../../package';


export class RGroupsAnalysisTutorial extends Tutorial {
  get name() {
    return 'R-Groups Analysis';
  }

  get description() {
    return 'R-Groups Analysis lets you identify all R-group ' +
    'variations around a scaffold, analyze substitution patterns, evaluate their ' +
    'impact on crucial compound properties, and find gaps.';
  }

  get steps() {return 14;}
  
  helpUrl: string = 'https://datagrok.ai/help/datagrok/solutions/domains/chem/#r-groups-analysis';
  prerequisites: TutorialPrerequisites = {packages: ['Chem'], jupyter: true};
  demoTable: string = '';
  // manualMode = true;

  protected async _run() {
    this.header.textContent = this.name;

    this.describe('<b>R-Groups Analysis</b> lets you identify all R-group ' +
    'variations around a scaffold, analyze substitution patterns, evaluate their ' +
    'impact on crucial compound properties, and find gaps.<hr>');

    this.t = await grok.data.loadTable(`${_package.webRoot}files/sar_small-R-groups.csv`);
    const tv = grok.shell.addTableView(this.t);

    this.title('Start the R-Groups Analysis tool', true);
    this.describe(`When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:<br>
    <ul>
    <li><b>Chem</b> menu (it contains all chemical tools)</li>
    </ul><br>
    Let’s launch the RGA tool.`);

    const d = await this.openDialog('On the Top Menu, click Chem > Analyze > R-Groups Analysis...', 'R-Groups Analysis',
      this.getMenuItem('Chem', true));

    this.title('Specify the scaffold', true);
    this.describe(`In the sketcher, you have two options to specify the scaffold:<br>
      <ul>
      <li>Manually draw or paste a scaffold.</li>
      <li>Click <b>MCS</b> to find the most common substructure.</li>
      </ul><br>
      For this tutorial, let’s choose <b>MCS</b>. Click it, then click <b>OK</b>`);

    await this.action('Click MCS', new Observable((subscriber: any) => {
      $('.chem-mcs-button').one('click', () => subscriber.next(true));
    }));

    await this.action('Click OK', d.onClose, $(d.root).find('button.ui-btn.ui-btn-ok')[1]);

    let v: DG.Viewer;
    await this.action('Wait for the analysis to complete',
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
        const found = data.args.viewer.type === DG.VIEWER.TRELLIS_PLOT;
        if (found)
          v = data.args.viewer;
        return found;
      })));
    
    this.title('Set up the visualization', true);
    this.describe(`Once the analysis is complete, the R-group columns are added to the table,
    along with a trellis plot for visual exploration.<br>Let’s set up the visualization.`);

    await this.action('In the trellis plot, click the gear icon for the embedded viewer',
      new Observable((subscriber: any) => {
        $('.grok-icon.grok-font-icon-settings').one('click', () => subscriber.next(true));
      }), v!.root.parentElement?.parentElement?.getElementsByClassName('grok-font-icon-settings')[0] as HTMLElement,
      `The <b>Context Panel</b> on the right now shows the settings for the trellis plot and the pie chart.`);
   
    grok.shell.windows.showContextPanel = true;
    await this.action('Under Pie chart tab > Data, set Category to LC/MS', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          if (m.previousSibling?.textContent === 'LC/MS') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {childList: true, subtree: true});
    }));

    this.title('Analyze and explore', true);
    this.describe(`All datagrok viewers are synchronized, and share the same filter and selection.
    The <b>Context Panel</b> provides information and actions relevant to your selection.<br>
    Let’s explore.`);

    await this.action('Click any segment on a pie chart', this.t.onSelectionChanged, undefined, 'Scroll to see the selected rows in the grid.');

    await this.action('Press Escape', new Observable((subscriber: any) => {
      document.addEventListener("keydown", ({key}) => {
        if (key === "Escape")
          subscriber.next(true);
      }, {once: true});
    }));

    await this.action('In the grid, press Shift+Drag Mouse Down', combineLatest([this.t.onSelectionChanged,
      new Observable((subscriber: any) => {
        const observer = new MutationObserver((mutationsList, observer) => {
          mutationsList.forEach((m) => {
            //@ts-ignore
            if (m.target.innerText.includes('Distributions')) {
              subscriber.next(true);
              observer.disconnect();
            }
          });
        });
        observer.observe($('.grok-prop-panel').get(0)!, {childList: true, subtree: true});
    })]), undefined, 'Note changes on the pie charts');

    // await this.action('In the grid, press Shift + Drag Mouse Down', this.t.onSelectionChanged.pipe(filter(() => {
    //   return grok.shell.o.constructor.name === 'RowGroup' && 
    //     $('.d4-accordion-pane-header').filter((_, el) => el.textContent === 'Distributions').length > 0;
    // })), undefined, 'Note changes on the pie charts');

    let pane: Cash;
    await this.action('On the Context Panel, expand the Distributions pane',
      new Observable((subscriber: any) => {
        pane = $('.d4-accordion-pane-header').filter((_, el) => el.textContent === 'Distributions');
        if (pane.hasClass('expanded'))
          subscriber.next(true);
        pane.on('click', () => subscriber.next(true));
      }), $('.d4-accordion-pane-header').filter((_, el) => el.textContent === 'Distributions').get(0));
    
    await this.action('In the pane, hover over line charts', new Observable((subscriber: any) => {
      const onMousemove = () => {
        if ($('.d4-tooltip').css('display') === 'block') {
          subscriber.next(true);
          document.querySelector('.grok-prop-panel')!.removeEventListener('mousemove', onMousemove);
        }
      };
      document.querySelector('.grok-prop-panel')!.addEventListener('mousemove', onMousemove);
    }));

    this.title('Get a different view', true);
    this.describe(`The trellis plot initially shows pie charts, but you can change visualizations to
    get a different view of your data. Let’s change a pie chart to a histogram and set it up.`);

    await this.action('In the top-left corner of the trellis plot, select Histogram.',
      v!.onEvent('d4-trellis-plot-viewer-type-changed').pipe(filter((s: string) => s === 'Histogram')),
      v!.root.querySelector('.d4-combo-popup') as HTMLElement);

    await this.action('Set Value to In-Vivo Activity', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          if (m.previousSibling?.textContent === 'In-vivo Activity') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {childList: true, subtree: true});
    }), undefined, `Use the <b>Gear</b> icon next to the <b>Viewer</b> control to access
      the histogram’s settings, and under <b>Histogram</b> tab set <b>Value</b> to <b>In-vivo Activity</b>.`);

    this.title('Switch axes', true);
    this.describe(`Finally, let’s change the R-groups used on the plot.`);

    const trellis = Array.from(grok.shell.tv.viewers).find((v) => v.type === DG.VIEWER.TRELLIS_PLOT)! as DG.Viewer<DG.ITrellisPlotSettings>;
    await this.action('Set the value for the X axis to R4',
      interval(1000).pipe(filter(() => {
        return trellis.props.xColumnNames[0] === 'R4';
    })));
  }
}
