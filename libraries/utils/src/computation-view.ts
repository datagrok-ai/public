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

    super();

    this.parentCall = grok.functions.getCurrentCall();
    this.parentView = grok.functions.getCurrentCall().parentCall.aux['view'];
    this.basePath = `/${grok.functions.getCurrentCall()?.func.name}`;

    ui.setUpdateIndicator(this.root, true);
    if (runId) {
      setTimeout(async () => {
        await this.init();
        this.linkFunccall(await this.loadRun(runId));
        this.build();
        this.linkFunccall(await this.loadRun(runId));
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
      grok.shell.o = this.historyRoot;
    });

    historyButton.classList.add('d4-toggle-button');
    if (grok.shell.windows.showProperties) historyButton.classList.add('d4-current');


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
    mainAcc.addTitle(ui.span(['History']));

    const buildFilterPane = () => {
      const dateInput = ui.stringInput('Date', 'Any time');
      dateInput.addPatternMenu('datetime');
      const form =ui.divV([
        ui.stringInput('Search', '', () => {}),
        ui.choiceInput('User', 'Current user', ['Current user']),
        dateInput,
      ], 'ui-form-condensed ui-form');
      form.style.marginLeft = '0px';

      return form;
    };
    let filterPane = mainAcc.addPane('Filter', buildFilterPane);
    const updateFilterPane = () => {
      const isExpanded = filterPane.expanded;
      mainAcc.removePane(filterPane);
      filterPane = mainAcc.addPane('Filter', buildFilterPane, isExpanded, favoritesListPane);
    };

    const showAddToFavoritesDialog = (funcCall: DG.FuncCall) => {
      let title = funcCall.options['title'] ?? '';
      let annotation = funcCall.options['annotation'] ?? '';
      const titleInput = ui.stringInput('Title', title, (s: string) => {
        title = s;
        if (s.length === 0) {
          titleInput.setTooltip('Title cannot be empty');
          setTimeout(() => titleInput.input.classList.add('d4-invalid'), 100);
        } else {
          titleInput.setTooltip('');
          setTimeout(() => titleInput.input.classList.remove('d4-invalid'), 100);
        }
      });

      ui.dialog({title: 'Add to favorites'})
        .add(ui.form([
          titleInput,
          ui.stringInput('Annotation', annotation, (s: string) => { annotation = s; }),
        ]))
        .onOK(async () => {
          if (title.length > 0) {
            funcCall.options['title'] = title;
            funcCall.options['annotation'] = annotation;
            await this.addRunToFavorites(funcCall);
            updateHistoryPane();
            updateFavoritesPane();
          } else {
            grok.shell.warning('Title cannot be empty');
          }
        })
        .show({center: true});
    };

    const showDeleteFavoritesDialog = (funcCall: DG.FuncCall) => {
      ui.dialog({title: 'Delete run'})
        .add(ui.divText('The deleted run is impossible to restore. Are you sure?'))
        .onOK(async () => {
          await this.deleteRun(funcCall);
          updateHistoryPane();
        })
        .show({center: true});
    };

    let historyCards = [] as HTMLElement[];
    let favoriteCards = [] as HTMLElement[];
    const renderFavoriteCards = async (funcCalls: DG.FuncCall[]) => {
      favoriteCards = funcCalls.map((funcCall) => {
        const solidStar = ui.iconFA('star', async (ev) => {
          ev.stopPropagation();
          await this.removeRunFromFavorites(funcCall);
          updateHistoryPane();
          updateFavoritesPane();
        }, 'Unfavorite the run');
        solidStar.classList.add('fas');

        const card = ui.divH([
          ui.divV([
            ui.divText(funcCall.options['title'] ?? 'Default title', 'title'),
            ...(funcCall.options['annotation']) ? [ui.divText(funcCall.options['annotation'], 'description')]: [],
            ui.divH([ui.render(funcCall.author), ui.span([new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})], 'date')]),
          ]),
          ui.divH([
            ui.iconFA('pen', async (ev) => {
              ev.stopPropagation();
              showAddToFavoritesDialog(funcCall);
            }, 'Edit run metadata'),
            solidStar
          ], 'funccall-card-icons')
        ], 'funccall-card');

        card.addEventListener('click', async () => {
          this.linkFunccall(await this.loadRun(funcCall.id));
          card.classList.add('clicked');
        });
        return card;
      });

      const allCards =[...historyCards, ...favoriteCards];
      allCards.forEach((card) => card.addEventListener('click', () => allCards.forEach((c) => c.classList.remove('clicked'))));

      return ui.divV(favoriteCards);
    };

    const renderHistoryCards = async (funcCalls: DG.FuncCall[]) => {
      historyCards = funcCalls.map((funcCall) => {
        const icon = funcCall.author.picture as HTMLElement;
        icon.style.width = '25px';
        icon.style.height = '25px';
        icon.style.fontSize = '20px';
        icon.style.marginRight = '3px';
        icon.style.alignSelf = 'center';
        const userLabel = ui.label(funcCall.author.friendlyName, 'd4-link-label');
        ui.bind(funcCall.author, userLabel);

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
              ev.stopPropagation();
              showAddToFavoritesDialog(funcCall);
            }, 'Pin the run'),
            ui.iconFA('link', async (ev) => {
              ev.stopPropagation();
              await navigator.clipboard.writeText(`${window.location.href}`);
            }, 'Copy link to the run'),
            ui.iconFA('trash-alt', async (ev) => {
              ev.stopPropagation();
              showDeleteFavoritesDialog(funcCall);
            }, 'Delete the run'),
          ], 'funccall-card-icons')
        ], 'funccall-card');

        card.addEventListener('click', async () => {
          this.linkFunccall(await this.loadRun(funcCall.id));
          card.classList.add('clicked');
        });

        return card;
      });

      const allCards =[...historyCards, ...favoriteCards];
      allCards.forEach((card) => card.addEventListener('click', () => allCards.forEach((c) => c.classList.remove('clicked'))));

      return ui.divV(historyCards);
    };

    const buildFavoritesList = () => ui.wait(async () => {
      const historicalRuns = await this.pullRuns(this.func!.id);
      const pinnedRuns = historicalRuns.filter((run) => run.options['isFavorite']);
      if (pinnedRuns.length > 0)
        return ui.wait(() => renderFavoriteCards(pinnedRuns));
      else
        return ui.divText('No runs are marked as favorites', 'description');
    });
    let favoritesListPane = mainAcc.addPane('Favorites', buildFavoritesList, true);
    const updateFavoritesPane = () => {
      const isExpanded = favoritesListPane.expanded;
      mainAcc.removePane(favoritesListPane);
      favoritesListPane = mainAcc.addPane('Favorites', buildFavoritesList, isExpanded, historyPane);
    };

    const buildHistoryPane = () => ui.wait(async () => {
      const historicalRuns = (await this.pullRuns(this.func!.id)).filter((run) => !run.options['isFavorite']);
      if (historicalRuns.length > 0)
        return ui.wait(() => renderHistoryCards(historicalRuns));
      else
        return ui.divText('No runs are found in history', 'description');
    });
    let historyPane = mainAcc.addPane('History', buildHistoryPane, true);
    const updateHistoryPane = () => {
      const isExpanded = historyPane.expanded;
      mainAcc.removePane(historyPane);
      historyPane = mainAcc.addPane('History', buildHistoryPane, isExpanded);
    };

    const newHistoryBlock = mainAcc.root;
    ui.empty(this.historyRoot);
    this.historyRoot.style.removeProperty('justify-content');
    this.historyRoot.style.width = '100%';
    this.historyRoot.append(newHistoryBlock);
    return newHistoryBlock;
  }
}
