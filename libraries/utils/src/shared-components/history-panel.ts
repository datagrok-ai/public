/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {BehaviorSubject, Subject} from 'rxjs';
import {defaultUsersIds} from '../function-view';
import {historyUtils} from '../history-utils';

export class HistoryPanel {
  // Emitted when FuncCall should is chosen. Contains FuncCall ID
  public onRunChosen = new Subject<string>();

  // Emitted when FuncCall is added to favorites
  public onRunAddToFavorites = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is added to shared
  public onRunAddToShared = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is removed from favorites
  public onRunRemoveFromFavorites = new Subject<string>();

  // Emitted when FuncCall is removed from shared
  public onRunRemoveFromShared = new Subject<string>();

  // Emitted when FuncCall is deleted
  public onRunDeleted = new Subject<string>();

  private store = {
    filteringOptions: {text: ''} as {
      text: string,
      author?: DG.User,
      startedAfter?: dayjs.Dayjs,
    },
    allRuns: new Subject<DG.FuncCall[]>,

    myRuns: new BehaviorSubject<DG.FuncCall[]>([]),
    favoriteRuns: new BehaviorSubject<DG.FuncCall[]>([]),
    sharedRuns: new BehaviorSubject<DG.FuncCall[]>([]),

    filteredMyRuns: new BehaviorSubject<DG.FuncCall[]>([]),
    filteredFavoriteRuns: new BehaviorSubject<DG.FuncCall[]>([]),
    filteredSharedRuns: new BehaviorSubject<DG.FuncCall[]>([]),
  };

  private myRunsFilter = new Subject<true>();
  private favRunsFilter = new Subject<true>();
  private sharedRunsFilter = new Subject<true>();

  public allRunsFetch = new Subject<true>();
  public myRunsFetch = new Subject<true>();
  public favRunsFetch = new Subject<true>();
  public sharedRunsFetch = new Subject<true>();

  historyTab = ui.div();
  favTab = ui.div();
  sharedTab = ui.div();
  tabs = ui.tabControl({
    'HISTORY': ui.box(this.historyTab),
    'FAVORITES': ui.box(this.favTab),
    'SHARED': ui.box(this.sharedTab),
  });
  filterPane = ui.div();

  private _root = ui.divV([
    ui.panel([
      ui.h2('History'),
      this.filterPane,
    ]),
    ui.element('div', 'splitbar-horizontal'),
    this.tabs.root,
  ], {style: {width: '100%'}});

  myCards = [] as HTMLElement[];
  favoriteCards = [] as HTMLElement[];
  sharedCards = [] as HTMLElement[];

  updateFilterPane() {
    const buildFilterPane = () => {
      return ui.wait(async () => {
        const filteringText = new Subject();

        const textInput = ui.stringInput('Search', '', (v: string) => filteringText.next(v));
        DG.debounce(filteringText.asObservable(), 600).subscribe(() => {
          this.store.filteringOptions.text = textInput.stringValue;
          this.favRunsFilter.next();
          this.sharedRunsFilter.next();
        });

        const dateInput = ui.dateInput('Started after', dayjs().subtract(1, 'week'), (v: dayjs.Dayjs) => {
          this.store.filteringOptions.startedAfter = v;
          this.myRunsFilter.next();
          this.favRunsFilter.next();
          this.sharedRunsFilter.next();
        });

        const defaultUsers = Object.values(defaultUsersIds);
        const allUsers = await grok.dapi.users.list() as DG.User[];
        const filteredUsers = allUsers.filter((user) => !defaultUsers.includes(user.id));

        const authorInput = ui.choiceInput<DG.User | string>('Author', 'Anyone', ['Anyone', ...filteredUsers], (v: DG.User | string) => {
          this.store.filteringOptions.author = (v === 'Anyone') ? undefined : v as DG.User;
          this.sharedRunsFilter.next();
        });
        dateInput.addPatternMenu('datetime');
        const form = ui.divV([
          textInput,
          dateInput,
          authorInput,
        ], 'ui-form-condensed ui-form');
        form.style.padding = '0px';

        this.tabs.onTabChanged.subscribe(() => {
          const currentTabName = this.tabs.currentPane.name;
          if (currentTabName === 'HISTORY') {
            textInput.root.style.display = 'none';
            dateInput.root.style.removeProperty('display');
            authorInput.root.style.display = 'none';
          }
          if (currentTabName === 'FAVORITES') {
            textInput.root.style.removeProperty('display');
            dateInput.root.style.removeProperty('display');
            authorInput.root.style.display = 'none';
          }
          if (currentTabName === 'SHARED') {
            textInput.root.style.removeProperty('display');
            dateInput.root.style.removeProperty('display');
            authorInput.root.style.removeProperty('display');
          }
        });

        return form;
      });
    };

    const newFilterPane = buildFilterPane();
    this.filterPane.replaceWith(newFilterPane);
    this.filterPane = newFilterPane;
  };

  showAddToFavoritesDialog(funcCall: DG.FuncCall) {
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
          funcCall = await historyUtils.loadRun(funcCall.id);
          funcCall.options['title'] = title;
          funcCall.options['annotation'] = annotation;

          this.onRunAddToFavorites.next(funcCall);
          this.myRunsFetch.next();
          this.favRunsFetch.next();
        } else {
          grok.shell.warning('Title cannot be empty');
        }
      })
      .show({center: true});
  };

  showDeleteRunDialog(funcCall: DG.FuncCall) {
    ui.dialog({title: 'Delete run'})
      .add(ui.divText('The deleted run is impossible to restore. Are you sure?'))
      .onOK(async () => {
        this.onRunDeleted.next(funcCall.id);
        this.myRunsFetch.next();
      })
      .show({center: true});
  };

  showAddToSharedDialog(funcCall: DG.FuncCall) {
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

    ui.dialog({title: 'Add to shared'})
      .add(ui.form([
        titleInput,
        ui.stringInput('Annotation', annotation, (s: string) => { annotation = s; }),
      ]))
      .onOK(async () => {
        if (title.length > 0) {
          funcCall = await historyUtils.loadRun(funcCall.id);
          funcCall.options['title'] = title;
          funcCall.options['annotation'] = annotation;
          this.onRunAddToShared.next(funcCall);
        } else {
          grok.shell.warning('Title cannot be empty');
        }
      })
      .show({center: true});
  };

  updateSharedPane(sharedRuns: DG.FuncCall[]) {
    const newTab = (sharedRuns.length > 0) ? this.renderSharedCards(sharedRuns) : ui.divText('No runs are marked as shared', 'no-elements-label');
    this.sharedTab.replaceWith(newTab);
    this.sharedTab = newTab;
  };

  updateFavoritesPane(favoriteRuns: DG.FuncCall[]) {
    const newTab = (favoriteRuns.length > 0) ? this.renderFavoriteCards(favoriteRuns) : ui.divText('No runs are marked as favorites', 'no-elements-label');
    this.favTab.replaceWith(newTab);
    this.favTab = newTab;
  };

  updateMyPane(myRuns: DG.FuncCall[]) {
    const newTab = (myRuns.length > 0) ? this.renderHistoryCards(myRuns) : ui.divText('No runs are found in history', 'no-elements-label');
    this.historyTab.replaceWith(newTab);
    this.historyTab = newTab;
  };

  renderFavoriteCards(funcCalls: DG.FuncCall[]) {
    this.favoriteCards = funcCalls.map((funcCall) => {
      const unstarIcon = ui.iconFA('star', async (ev) => {
        ev.stopPropagation();
        this.onRunRemoveFromFavorites.next(funcCall.id);
      }, 'Unfavorite the run');
      unstarIcon.classList.add('fas');

      const shareIcon = ui.iconFA('eye', async (ev) => {
        ev.stopPropagation();
        this.showAddToSharedDialog(funcCall);
      }, 'Add to shared');
      shareIcon.classList.add('fal');

      const card = ui.divH([
        ui.divV([
          ui.divText(funcCall.options['title'] ?? 'Default title', 'title'),
          ...(funcCall.options['annotation']) ? [ui.divText(funcCall.options['annotation'], 'description')]: [],
          ui.divH([ui.render(funcCall.author), ui.span([new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})], 'date')]),
        ], 'cv-card-content'),
        ui.divH([
          shareIcon,
          ui.iconFA('pen', async (ev) => {
            ev.stopPropagation();
            this.showAddToFavoritesDialog(funcCall);
          }, 'Edit run metadata'),
          unstarIcon
        ], 'cv-funccall-card-icons')
      ], 'cv-funccall-card');

      card.addEventListener('click', async () => {
        this.onRunChosen.next(funcCall.id);
        card.classList.add('clicked');
      });
      return card;
    });

    const allCards = [...this.myCards, ...this.favoriteCards, ...this.sharedCards];
    allCards.forEach((card) => card.addEventListener('click', () => allCards.forEach((c) => c.classList.remove('clicked'))));

    return ui.divV(this.favoriteCards);
  };

  renderSharedCards(funcCalls: DG.FuncCall[]) {
    this.sharedCards = funcCalls.map((funcCall) => {
      const unshareIcon = ui.iconFA('eye-slash', async (ev) => {
        ev.stopPropagation();
        this.onRunRemoveFromShared.next(funcCall.id);
      }, 'Hide from shared');

      const card = ui.divH([
        ui.divV([
          ui.divText(funcCall.options['title'] ?? 'Default title', 'title'),
          ...(funcCall.options['annotation']) ? [ui.divText(funcCall.options['annotation'], 'description')]: [],
          ui.divH([ui.render(funcCall.author), ui.span([new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})], 'date')]),
        ], 'cv-card-content'),
        ui.divH([
          ui.iconFA('link', async (ev) => {
            ev.stopPropagation();
            await navigator.clipboard.writeText(`${window.location.origin}${window.location.pathname}?id=${funcCall.id}`);
          }, 'Copy link to the run'),
          ...(funcCall.author.id === grok.shell.user.id) ? [
            ui.iconFA('pen', async (ev) => {
              ev.stopPropagation();
              this.showAddToSharedDialog(funcCall);
            }, 'Edit run metadata'),
            unshareIcon]: [],
        ], 'cv-funccall-card-icons')
      ], 'cv-funccall-card');

      card.addEventListener('click', async () => {
        this.onRunChosen.next(funcCall.id);
        card.classList.add('clicked');
      });
      return card;
    });


    const allCards = [...this.myCards, ...this.favoriteCards, ...this.sharedCards];
    allCards.forEach((card) => card.addEventListener('click', () => allCards.forEach((c) => c.classList.remove('clicked'))));

    return ui.divV(this.sharedCards);
  };

  renderHistoryCards(funcCalls: DG.FuncCall[]) {
    this.myCards = funcCalls.map((funcCall) => {
      const icon = funcCall.author.picture as HTMLElement;
      icon.style.width = '25px';
      icon.style.height = '25px';
      icon.style.fontSize = '20px';
      icon.style.marginRight = '3px';
      icon.style.alignSelf = 'center';
      const userLabel = ui.label(funcCall.author.friendlyName, 'd4-link-label');
      ui.bind(funcCall.author, icon);

      const shareIcon = ui.iconFA('eye', async (ev) => {
        ev.stopPropagation();
        this.showAddToSharedDialog(funcCall);
      }, 'Add to shared');
      shareIcon.classList.add('fal');

      const unshareIcon = ui.iconFA('eye-slash', async (ev) => {
        ev.stopPropagation();
        this.onRunRemoveFromShared.next(funcCall.id);
      }, 'Hide from shared');

      const unstar = ui.iconFA('star', async (ev) => {
        ev.stopPropagation();
        this.onRunRemoveFromFavorites.next(funcCall.id);
      }, 'Unfavorite the run');
      unstar.classList.add('fas');

      const card = ui.divH([
        ui.divH([
          icon,
          ui.divV([
            userLabel,
            ui.span([new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})])
          ], 'cv-card-content'),
        ]),
        ui.divH([
          ...(funcCall.options['isShared']) ? [unshareIcon]: [shareIcon],
          ...(funcCall.options['isFavorite']) ? [unstar] : [ui.iconFA('star', async (ev) => {
            ev.stopPropagation();
            this.showAddToFavoritesDialog(funcCall);
          }, 'Add to favorites')],
          ui.iconFA('trash-alt', async (ev) => {
            ev.stopPropagation();
            this.showDeleteRunDialog(funcCall);
          }, 'Delete the run'),
        ], 'cv-funccall-card-icons')
      ], 'cv-funccall-card');

      card.addEventListener('click', async () => {
        this.onRunChosen.next(funcCall.id);
        card.classList.add('clicked');
      });

      return card;
    });

    const allCards = [...this.myCards, ...this.favoriteCards, ...this.sharedCards];
    allCards.forEach((card) => card.addEventListener('click', () => allCards.forEach((c) => c.classList.remove('clicked'))));

    return ui.divV(this.myCards);
  };

  constructor(
    private func: DG.Func
  ) {
    this.tabs.root.style.width = '100%';
    this.tabs.header.style.justifyContent = 'space-between';

    this.store.myRuns.subscribe((myRuns) => this.store.filteredMyRuns.next(myRuns));
    this.store.favoriteRuns.subscribe((favoriteRuns) => this.store.filteredFavoriteRuns.next(favoriteRuns));
    this.store.sharedRuns.subscribe((sharedRuns) => this.store.filteredSharedRuns.next(sharedRuns));

    this.store.filteredMyRuns.subscribe((myRuns) => this.updateMyPane(myRuns));
    this.store.filteredFavoriteRuns.subscribe((newFavoriteRuns) => this.updateFavoritesPane(newFavoriteRuns));
    this.store.filteredSharedRuns.subscribe((sharedRuns) => this.updateSharedPane(sharedRuns));

    this.myRunsFilter.subscribe(() => {
      const filteredMyRuns = this.store.myRuns.value.filter((val) => {
        const startedAfter = (this.store.filteringOptions.startedAfter ? val.started > this.store.filteringOptions.startedAfter : true);

        return startedAfter;
      });
      this.store.filteredMyRuns.next(filteredMyRuns);
    });

    this.favRunsFilter.subscribe(() => {
      const filteredFavRuns = this.store.favoriteRuns.value.filter((val) => {
        const startedAfter = (this.store.filteringOptions.startedAfter ? val.started > this.store.filteringOptions.startedAfter : true);
        const titleContainsText = (this.store.filteringOptions.text.length && !!val.options['title']) ? val.options['title'].includes(this.store.filteringOptions.text): true;
        const descContainsText = (this.store.filteringOptions.text.length && !!val.options['description']) ? val.options['description'].includes(this.store.filteringOptions.text): true;

        return (titleContainsText || descContainsText) && startedAfter;
      });
      this.store.filteredFavoriteRuns.next(filteredFavRuns);
    });

    this.sharedRunsFilter.subscribe(() => {
      const filteredSharedRuns = this.store.sharedRuns.value.filter((val) => {
        const isAuthored = this.store.filteringOptions.author ? val.author.id === this.store.filteringOptions.author.id: true;
        const startedAfter = (this.store.filteringOptions.startedAfter ? val.started > this.store.filteringOptions.startedAfter : true);
        const titleContainsText = (this.store.filteringOptions.text.length && !!val.options['title']) ? val.options['title'].includes(this.store.filteringOptions.text): true;
        const descContainsText = (this.store.filteringOptions.text.length && !!val.options['description']) ? val.options['description'].includes(this.store.filteringOptions.text): true;

        return (titleContainsText || descContainsText) && startedAfter && isAuthored;
      });
      this.store.filteredSharedRuns.next(filteredSharedRuns);
    });

    this.store.allRuns.subscribe((allRuns) => {
      this.store.myRuns.next(allRuns.filter((run) => run.author.id === grok.shell.user.id));
      this.store.favoriteRuns.next(allRuns.filter((run) => run.options['isFavorite'] && !run.options['isImported']));
      this.store.sharedRuns.next(allRuns.filter((run) => run.options['isShared']));
    });

    this.allRunsFetch.subscribe(async () => {
      const allRuns = (await historyUtils.pullRunsByName(this.func.name, [], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.store.allRuns.next(allRuns);
    });

    this.myRunsFetch.subscribe(async () => {
      const myRuns = (await historyUtils.pullRunsByName(this.func.name, [{author: grok.shell.user}], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.store.myRuns.next(myRuns);
    });

    this.favRunsFetch.subscribe(async () => {
      const myRuns = (await historyUtils.pullRunsByName(this.func.name, [{author: grok.shell.user}], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.store.favoriteRuns.next(myRuns.filter((run) => run.options['isFavorite'] && !run.options['isImported']));
    });

    this.sharedRunsFetch.subscribe(async () => {
      const sharedRuns = (await historyUtils.pullRunsByName(this.func.name, [], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.store.sharedRuns.next(sharedRuns.filter((run) => run.options['isShared']));
    });

    this.allRunsFetch.next();
    this.updateFilterPane();
  }

  get root() {
    return this._root;
  }
}
