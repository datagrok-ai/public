/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView} from './function-view';
import '../css/computation-view.css';

/**
 * Base class for handling Compute models (see https://github.com/datagrok-ai/public/blob/master/help/compute/compute.md).
 * In most cases, the computation is a simple {@link Func}
 * Extend it in cases where a behavior or UI not supported by the {@link FunctionView} is needed.
 *
 * It provides the following functionality out-of-the-box, where each section could be customized:
 * - a structured way to represent input and output parameters: {@link parameters}
 * - generic way to generate UI for inputs, outputs, and interactivity (running the model, etc)
 *   - persisting historical results to the db (via {@link parameters})
 * - export (to Excel and PDF): {@link export}
 * - easy loading of historical runs
 * - routing
 * - entering the real, measured (as opposed to predicted) values manually
 * - notifications for changing inputs, completion of computations, etc: {@link onInputChanged}
 * */
export class ComputationView extends FunctionView {
  /** Find the function by fully specified name, link it to the view and constructs the view.
    * If function name is not specified, calls {@link init} and {@link build} without FuncCall linkage.
    * @stability Stable
  */
  constructor(funcName?: string) {
    const url = new URL(grok.shell.startUri);
    const runId = url.searchParams.get('id');
    console.log(grok.shell.startUri);
    console.log('url ', url);
    console.log('runId ', runId);

    super();

    this.parentCall = grok.functions.getCurrentCall();
    this.parentView = grok.functions.getCurrentCall().parentCall.aux['view'];
    this.basePath = `/${grok.functions.getCurrentCall()?.func.name}`;

    ui.setUpdateIndicator(this.root, true);
    if (runId) {
      setTimeout(async () => {
        this.linkFunccall(await this.loadRun(runId));
        await this.init();
        this.build();
        ui.setUpdateIndicator(this.root, false);
      }, 0);
    } else {
      if (funcName) {
        grok.functions.eval(funcName).then(async (func: DG.Func) => {
          const funccall = func.prepare({});
          this.linkFunccall(funccall);
          await this.init();
          this.build();
        }).finally(() => {
          ui.setUpdateIndicator(this.root, false);
        });
      } else {
        setTimeout(async () => {
          await this.init();
          this.build();
          ui.setUpdateIndicator(this.root, false);
        }, 0);
      }
    }

    grok.shell.o = this.historyRoot;
  }

  /** Override to customize getting mocks
    * @stability Experimental
  */
  getMocks: ({mockName: string, action: () => Promise<void>}[]) | null = null;

  /** Override to customize getting templates
    * @stability Experimental
  */
  getTemplates: ({name: string, action: () => Promise<void>}[]) | null = null;

  /** Override to customize getting help feature
    * @stability Stable
  */
  getHelp: (() => Promise<void>) | null = null;

  /** Override to customize bug reporting feature
    * @stability Stable
  */
  reportBug: (() => Promise<void>) | null = null;

  override buildRibbonPanels() {
    super.buildRibbonPanels();

    const historyButton = ui.iconFA('history', () => {
      grok.shell.windows.showProperties = !grok.shell.windows.showProperties;
      historyButton.classList.toggle('d4-current');
    });

    historyButton.classList.add('d4-toggle-button');

    const newRibbonPanels = [
      ...this.getRibbonPanels(),
      [ui.div(historyButton)] as HTMLElement[]
    ];
    this.setRibbonPanels(newRibbonPanels);
    return newRibbonPanels;
  }

  /**
   * Looks for {@link reportBug}, {@link getHelp} and {@link exportConfig} members and creates model menus
   * @stability Stable
  */
  override buildRibbonMenu() {
    super.buildRibbonMenu();

    if (!this.exportConfig && !this.reportBug && !this.getHelp && !this.getMocks && !this.getTemplates) return;

    const ribbonMenu = this.ribbonMenu.group('Model');

    if (this.getMocks && this.getMocks.length > 0) {
      if (this.getMocks.length === 1) {
        ribbonMenu.item('Input data mock', this.getMocks[0].action);
      } else {
        const dataGroup = ribbonMenu.group('Input data mocks');
        this.getMocks.forEach((val) => {
          dataGroup.item(val.mockName, val.action);
        });
        ribbonMenu.endGroup();
      }
    }

    if (this.getTemplates && this.getTemplates.length > 0) {
      if (this.getTemplates.length === 1) {
        ribbonMenu.item('Input data template', this.getTemplates[0].action);
      } else {
        const dataGroup = ribbonMenu.group('Input data templates');
        this.getTemplates.forEach((val) => {
          dataGroup.item(val.name, val.action);
        });
        ribbonMenu.endGroup();
      }
    }

    if (this.exportConfig && this.exportConfig.supportedFormats.length > 0) {
      ribbonMenu
        .group('Export')
        .items(this.exportConfig.supportedFormats, async (format: string) => DG.Utils.download(this.exportConfig!.filename(format), await this.exportConfig!.export(format)))
        .endGroup();
    }

    if (this.reportBug)
      ribbonMenu.item('Report a bug', () => this.reportBug!());

    if (this.getHelp)
      ribbonMenu.item('Help', () => this.getHelp!());
  }

  override buildHistoryBlock(): HTMLElement {
    const mainAcc = ui.accordion();
    mainAcc.root.style.width = '100%';
    mainAcc.addTitle(ui.h1('History'));
    const dateInput = ui.stringInput('Date', 'Any time');
    dateInput.addPatternMenu('datetime');

    mainAcc.addPane('Filter', () => {
      const form =ui.divV([
        ui.choiceInput('User', 'Current user', ['Current user']),
        dateInput,
      ], 'ui-form-condensed ui-form');
      form.style.marginLeft = '0px';

      return form;
    });

    const showPinDialog = (funcCall: DG.FuncCall) => {
      let title = funcCall.options['title'] ?? '';
      let annotation = funcCall.options['annotation'] ?? '';

      ui.dialog({title: 'Pin run'})
        .add(ui.form([
          ui.stringInput('Title', title, (s: string) => { title = s; }, {clearIcon: true, placeholder: 'Enter run title here...'}),
          ui.stringInput('Annotation', annotation, (s: string) => { annotation = s; }, {clearIcon: true, placeholder: 'Enter annotation title here...'}),
        ]))
        .onOK(async () => {
          if (title.length > 0) {
            funcCall.options['title'] = title;
            funcCall.options['annotation'] = annotation;
            await this.pinRun(funcCall);
          } else { showPinDialog(funcCall); }
        })
        .show();
    };

    const renderPinnedCard = async (funcCall: DG.FuncCall) => {
      // FIXME: black magic to get author
      const author = DG.toJs(funcCall.dart.x2.a);

      const card = ui.divH([
        ui.divV([
          ui.divText(funcCall.options['title'] ?? 'Default title', 'title'),
          ...(funcCall.options['annotation']) ? [ui.divText(funcCall.options['annotation'], 'description')]: [],
          ui.divH([ui.render(author), ui.span([new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})], 'date')]),
        ]),
        ui.divH([
          ui.iconFA('pen', async (ev) => {
            showPinDialog(funcCall);
            ev.stopPropagation();
          }, 'Edit run metadata')
        ], 'funccall-card-icons')
      ], 'funccall-card');

      return card;
    };

    const renderHistoryCard = async (funcCall: DG.FuncCall) => {
      // FIXME: black magic to get author
      const author: DG.User = DG.toJs(funcCall.dart.x2.a);
      const icon = author.picture as HTMLElement;
      icon.style.width = '25px';
      icon.style.height = '25px';
      icon.style.fontSize = '20px';
      icon.style.marginRight = '3px';
      icon.style.alignSelf = 'center';
      const userLabel = ui.label(author.friendlyName, 'd4-link-label');
      ui.bind(author, userLabel);

      const card = ui.divH([
        ui.divH([
          icon,
          ui.divV([
            userLabel,
            ui.span([new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})])
          ]),
        ]),
        ui.divH([
          ui.iconFA('star', async (ev) => {
            showPinDialog(funcCall);
            ev.stopPropagation();
          }, 'Pin the run'),
          ui.iconFA('share-alt', async (ev) => {
            const url = grok.shell.startUri.substring(0, grok.shell.startUri.indexOf(`?id`));
            await navigator.clipboard.writeText(`${url}?id=${funcCall.id}`);
            ev.stopPropagation();
          }, 'Copy link to the run'),
          ui.iconFA('trash-alt', async (ev) => {
            await this.deleteRun(funcCall);
            ev.stopPropagation();
          }, 'Delete the run'),
        ], 'funccall-card-icons')
      ], 'funccall-card');

      card.addEventListener('click', async () => {
        await this.loadRun(funcCall.id);
      });

      return card;
    };

    mainAcc.addPane('Pinned', () => ui.wait(async () => {
      const historicalRuns = await this.pullRuns(this.func!.id);
      const pinnedRuns = historicalRuns.filter((run) => run.options['isPinned']);
      if (pinnedRuns.length > 0)
        return ui.divV(pinnedRuns.map((run) => ui.wait(() => renderPinnedCard(run))));
      else
        return ui.divText('Pin run in history panel to have quick access to it', 'description');
    }));

    mainAcc.addPane('History', () => ui.wait(async () => {
      const historicalRuns = await this.pullRuns(this.func!.id);
      return ui.divV(historicalRuns.map((run) => ui.wait(() => renderHistoryCard(run))));
    }));

    const newHistoryBlock = mainAcc.root;
    ui.empty(this.historyRoot);
    this.historyRoot.style.removeProperty('justify-content');
    this.historyRoot.style.width = '100%';
    this.historyRoot.append(newHistoryBlock);
    return newHistoryBlock;
  }
}
