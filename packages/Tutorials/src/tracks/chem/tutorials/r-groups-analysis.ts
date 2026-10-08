import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable, combineLatest, fromEvent, interval} from 'rxjs';
import $ from 'cash-dom';
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

  get steps() {return 15;}

  get icon() {
    return '🧬🧩';
  }

  helpUrl: string = 'https://datagrok.ai/help/datagrok/solutions/domains/chem/#r-groups-analysis';
  prerequisites: TutorialPrerequisites = {packages: ['Chem']};
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
      () => this.getMenuItem('Chem', true), '',
      () => grok.shell.v.ribbonMenu.find('Chem | Analyze | R-Groups Analysis...').click());

    this.title('Specify the scaffold', true);
    this.describe(`In the sketcher, you have two options to specify the scaffold:<br>
      <ul>
      <li>Manually draw or paste a scaffold.</li>
      <li>Click <b>MCS</b> to find the most common substructure.</li>
      </ul><br>
      For this tutorial, let’s choose <b>MCS</b>. Click it, then click <b>OK</b>`);

    await this.action('Click MCS', new Observable((subscriber: any) => {
      $('.chem-mcs-button').one('click', () => subscriber.next(true));
    }), undefined, '', () => $(d.root).find('.chem-mcs-button')[0]!.click());

    await this.action('Click OK', d.onClose, $(d.root).find('button.ui-btn.ui-btn-ok')[1], '', async () => {
      await Tutorial.waitFor(() => $(d.root).find('.d4-update-shadow')[0] ? null : true);
      d.getButton('OK').click();
    });

    let v: DG.Viewer;
    await this.action('Wait for the analysis to complete',
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
        const found = data.args.viewer.type === DG.VIEWER.TRELLIS_PLOT;
        if (found)
          v = data.args.viewer;
        return found;
      })), null, '', null);

    this.title('Set up the visualization', true);
    this.describe(`Once the analysis is complete, the R-group columns are added to the table,
    along with a trellis plot for visual exploration.<br>Let’s set up the visualization.`);

    // the gear sits in the trellis title bar, outside its root; only this viewer's gear counts
    const trellisGear = () => v!.root.parentElement?.parentElement?.getElementsByClassName('grok-font-icon-settings')[0] as HTMLElement ?? null;
    await this.action('In the trellis plot, click the gear icon for the embedded viewer',
      fromEvent(document, 'click').pipe(filter((e) => trellisGear() != null && trellisGear().contains(e.target as Node))), trellisGear,
      `The <b>Context Panel</b> on the right now shows the settings for the trellis plot and the pie chart.`,
      () => trellisGear()!.click());

    grok.shell.windows.showContextPanel = true;
    const pieChartTab = () => document.querySelector('.grok-prop-panel [name="Pie chart"]') as HTMLElement ?? null;
    await this.action('Go to Pie chart tab',
      interval(200).pipe(filter(() => pieChartTab()?.classList.contains('selected') ?? false)), pieChartTab, '',
      () => pieChartTab()!.click());

    const innerLook = () => (v!.getOptions().look as any).innerViewerLook ?? {};
    await this.action('Under Pie chart tab > Data, set Category to LC/MS',
      interval(200).pipe(filter(() => innerLook().categoryColumnName === 'LC/MS')), undefined, '',
      () => v!.setOptions({innerViewerLook: {...innerLook(), categoryColumnName: 'LC/MS'}}));

    this.title('Analyze and explore', true);
    this.describe(`All datagrok viewers are synchronized, and share the same filter and selection.
    The <b>Context Panel</b> provides information and actions relevant to your selection.<br>
    Let’s explore.`);

    await this.action('Click any segment on a pie chart', this.t.onSelectionChanged, undefined, 'Scroll to see the selected rows in the grid.',
      () => this.t!.selection.init((i) => this.t!.get('LC/MS', i) === 'pass'));

    await this.action('Press Escape', new Observable((subscriber: any) => {
      document.addEventListener("keydown", ({key}) => {
        if (key === "Escape")
          subscriber.next(true);
      }, {once: true});
    }), undefined, '', Tutorial.apiSkip(() => this.t!.selection.setAll(false)));

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

    // the pane exists only while a RowGroup is the current object, which the previous step
    // establishes — so it is resolved per tick, not once when this step is built
    const distrHeader = (): HTMLElement | null => $('.d4-accordion-pane-header')
      .filter((_, el) => el.textContent === 'Distributions').get(0) ?? null;

    await this.action('On the Context Panel, expand the Distributions pane',
      new Observable((subscriber: any) => {
        const header = distrHeader();
        if (header != null && $(header).hasClass('expanded')) {
          subscriber.next(true);
          return;
        }
        // delegated, so it does not matter whether the pane is in the DOM yet
        const onClick = (e: Event) => {
          const el = (e.target as HTMLElement).closest('.d4-accordion-pane-header');
          if (el != null && el.textContent === 'Distributions')
            subscriber.next(true);
        };
        document.addEventListener('click', onClick, true);
        return () => document.removeEventListener('click', onClick, true);
      }), distrHeader, '', Tutorial.apiSkip(() => distrHeader()?.click()));

    const distrLineChart = () => document.querySelector('.d4-pane-distributions.expanded canvas') as HTMLElement ?? null;
    // the pane is rebuilt as the selection changes, so the pointer is followed on the document; the
    // tooltip shows once the pointer rests, after the last move
    let overPane = false;
    const paneMoves = fromEvent(document, 'mousemove')
      .subscribe((e) => overPane = (e.target as HTMLElement).closest?.('.d4-pane-distributions') != null);
    await this.action('In the pane, hover over line charts to see distributions',
      interval(200).pipe(filter(() => overPane && $('.d4-tooltip').css('display') === 'block')), distrLineChart);
    paneMoves.unsubscribe();

    this.title('Get a different view', true);
    this.describe(`The trellis plot initially shows pie charts, but you can change visualizations to
    get a different view of your data. Let’s change a pie chart to a histogram and set it up.`);

    await this.action('In the top-left corner of the trellis plot, select Histogram.',
      interval(200).pipe(filter(() => v!.props.viewerType === 'Histogram')),
      v!.root.querySelector('.d4-combo-popup') as HTMLElement, '', () => v!.setOptions({viewerType: 'Histogram'}));

    await this.action('Set Value to In-Vivo Activity',
      interval(200).pipe(filter(() => innerLook().valueColumnName === 'In-vivo Activity')), undefined, `Use the <b>Gear</b> icon next to the <b>Viewer</b> control to access
      the histogram’s settings, and under <b>Histogram</b> tab set <b>Value</b> to <b>In-vivo Activity</b>.`,
      () => v!.setOptions({innerViewerLook: {...innerLook(), valueColumnName: 'In-vivo Activity'}}));

    this.title('Switch axes', true);
    this.describe(`Finally, let’s change the R-groups used on the plot.`);

    const trellis = Array.from(grok.shell.tv.viewers).find((v) => v.type === DG.VIEWER.TRELLIS_PLOT)! as DG.Viewer<DG.ITrellisPlotSettings>;
    await this.action('Set the value for the X axis to R4',
      interval(1000).pipe(filter(() => {
        return trellis.props.xColumnNames[0] === 'R4';
      })), undefined, '', () => trellis.setOptions({xColumnNames: ['R4']}));
  }
}
