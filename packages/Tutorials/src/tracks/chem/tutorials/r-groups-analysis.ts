import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import {Observable} from 'rxjs';
import $, { Cash } from 'cash-dom';


export class RGroupsAnalysisTutorial extends Tutorial {
  get name() {
    return 'R-Groups Analysis';
  }

  get description() {
    return 'R-Groups Analysis lets you identify all R-group ' +
    'variations around a scaffold, analyze substitution patterns, evaluate their ' +
    'impact on crucial compound properties, and find gaps.';
  }

  get steps() { return 12; }
  
  helpUrl: string = '';

  protected async _run() {
    this.header.textContent = this.name;

    this.describe(this.description);

    const df = await grok.data.files.openTable('System:AppData/Tutorials/sar_small-R-groups.csv');
    const tv = grok.shell.addTableView(df);

    const d = await this.openDialog('Start the R-Groups Analysis tool', 'R-Groups Analysis',
    this.getMenuItem('Chem'), `When you open a chemical dataset, Datagrok automatically detects molecules
    and shows molecule-specific tools, actions, and information. Access them through:
    <ul>
    <li><b>Chem</b> menu (it houses all chemical tools).</li>
    <li>Context menu (right-click for access).</li>
    <li><b>Context Panel</b> on the right.</li>
    </ul><br>
    Let’s launch the R-Groups Analysis tool. On the <b>Top Menu</b>, click <b>Chem</b> > <b>Analyze</b> > <b>R-Groups Analysis...</b>`);

    await this.action('Specify the scaffold', new Observable((subscriber: any) => {
        $('.chem-mcs-button').one('click', () => subscriber.next(true));
      }), undefined, `In the sketcher, you have two options to specify the scaffold:
      <ul>
      <li>Manually draw or paste a scaffold.</li>
      <li>Click <b>MCS</b> to automatically detect the most common substructure.</li>
      </ul>
      For this tutorial, let’s choose <b>MCS</b>.`);

    await this.action('Click OK', d.onClose, $(d.root).find('button.ui-btn.ui-btn-ok')[1]);

    let v: DG.Viewer;
    await this.action('Wait for the analysis to complete',
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
        const found = data.args.viewer.type === DG.VIEWER.TRELLIS_PLOT;
        if (found)
          v = data.args.viewer;
        return found;
    })));
    
    await this.action('Set up the visualization',
    new Observable((subscriber: any) => {
      $('.grok-icon.grok-font-icon-settings.d4-viewer-icon').one('click', () => subscriber.next(true));
    }), $('.grok-icon.grok-font-icon-settings.d4-viewer-icon').get(0), `
    Once the analysis is complete, the R-group columns are added to the table (grid), along with a trellis plot for visual exploration.<br>
    Let’s set up the visualization. In the top-left corner, click the <b>Gear</b> icon.<br>
    The <b>Context Panel</b> on the right now shows the settings for the pie chart.<br>`);

   
    await this.action('Set Category to LC/MS', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          if (m.previousSibling?.textContent === 'LC/MS') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {childList: true, subtree: true});
    }), undefined,
    'On the Context Panel, under <b>Data</b> set <b>Category</b> to <b>LC/MS</b>.');

    await this.action('Analyze and explore', df.onSelectionChanged.pipe(filter(() => {
      return grok.shell.o.constructor.name === 'RowGroup' && 
        $('.d4-accordion-pane-header').filter((ind, el) => el.textContent === 'Distributions').length > 0;
    })), v!.root, `All datagrok viewers are synchronized. This means they respond to the same filtering, selection, and so on.
    Click any segment on the pie chart. Note the selected items are instantly highlighted on the grid.
    <a href="https://datagrok.ai/help/visualize/viewers/pie-chart">Learn more about pie charts</a>`);

    let pane: Cash;
    await this.action('On the Context Panel expand the Distributions pane',
    new Observable((subscriber: any) => {
      pane = $('.d4-accordion-pane-header').filter((ind, el) => el.textContent === 'Distributions');
      if (pane.hasClass('expanded'))
        subscriber.next(true);
      pane.on('click', () => subscriber.next(true));
    }), $('.d4-accordion-pane-header').filter((ind, el) => el.textContent === 'Distributions').get(0));
    
    await this.action('Hover over any line chart in that info pane',
    new Observable((subscriber: any) => {
      const content = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Distributions')?.nextElementSibling;
      content?.querySelectorAll('.d4-flex-col')[2].addEventListener('mouseover', () => subscriber.next(true), {once: true});
    }));

    await this.action('Get a different view',
      v!.onEvent('d4-trellis-plot-viewer-type-changed').pipe(filter((s: string) => s === 'Histogram')),
      v!.root.querySelector('.d4-combo-popup') as HTMLElement, `The trellis plot initially shows pie charts, but you can change
      visualizations to get a different view of your data. Let’s change a pie chart to a histogram. In the top-left corner,
      click the <b>Viewer</b> control and select <b>Histogram</b>.`);

    await this.action('Set <b>Value</b> to <b>In-Vivo Activity</b>', new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          if (m.previousSibling?.textContent === 'In-vivo Activity') {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe($('.grok-prop-panel').get(0)!, {childList: true, subtree: true});
    }), undefined, `Now, use the <b>Gear</b> icon next to the <b>Viewer</b> control to access
      the histogram’s settings, and set <b>Value</b> to <b>In-vivo Activity</b>.`);

    await this.action('Switch axes', new Observable((subscriber: any) => {
      const s = v.root.querySelectorAll('.d4-columns-host')[0].children[1] as HTMLSelectElement;
      s.addEventListener('change', (e) => {
        //@ts-ignore
        if (e.target.value === 'R4')
          subscriber.next(true)
      });
    }), undefined, `Finally, let’s change the R-groups used on the plot. In the top left corner,
    set the <b>X</b> axis to <b>R4</b>.<br>`);
  }
}
