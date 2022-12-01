/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {BehaviorSubject, Observable, Subject} from 'rxjs';
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

  mainAcc = ui.accordion();
  filterPane = this.mainAcc.addPane('Filter', () => ui.div());
  sharedListPane = this.mainAcc.addPane('Shared', () => ui.div(), true);
  favoritesListPane = this.mainAcc.addPane('My favorites', () => ui.div(), true);
  historyPane = this.mainAcc.addPane('My history', () => ui.div(), true);

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
        const allUsers = await grok.dapi.users.list();
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
        form.style.marginLeft = '0px';

        return form;
      });
    };

    const isExpanded = this.filterPane.expanded;
    this.mainAcc.removePane(this.filterPane);
    this.filterPane = this.mainAcc.addPane('Filter', () => buildFilterPane(), isExpanded, this.sharedListPane);
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
    const isExpanded = this.sharedListPane.expanded;
    this.mainAcc.removePane(this.sharedListPane);
    this.sharedListPane = this.mainAcc.addPane('Shared', () => {
      if (sharedRuns.length > 0)
        return this.renderSharedCards(sharedRuns);
      else
        return ui.divText('No runs are marked as shared', 'description');
    }, isExpanded, this.favoritesListPane);
  };

  updateFavoritesPane(favoriteRuns: DG.FuncCall[]) {
    const isExpanded = this.favoritesListPane.expanded;
    this.mainAcc.removePane(this.favoritesListPane);
    this.favoritesListPane = this.mainAcc.addPane('My favorites', () => {
      if (favoriteRuns.length > 0)
        return this.renderFavoriteCards(favoriteRuns);
      else
        return ui.divText('No runs are marked as favorites', 'description');
    }, isExpanded, this.historyPane);
  };

  updateMyPane(myRuns: DG.FuncCall[]) {
    const isExpanded = this.historyPane.expanded;
    this.mainAcc.removePane(this.historyPane);
    this.historyPane = this.mainAcc.addPane('My history', () => {
      if (myRuns.length > 0)
        return this.renderHistoryCards(myRuns);
      else
        return ui.divText('No runs are found in history', 'description');
    }, isExpanded);
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
        ]),
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
        ]),
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
          ]),
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
    this.mainAcc.root.style.width = '100%';
    this.mainAcc.addTitle(ui.span(['History']));

    this.store.myRuns.subscribe((myRuns) => this.store.filteredMyRuns.next(myRuns));
    this.store.favoriteRuns.subscribe((favoriteRuns) => this.store.filteredFavoriteRuns.next(favoriteRuns));
    this.store.sharedRuns.subscribe((sharedRuns) => this.store.filteredSharedRuns.next(sharedRuns));

    this.store.filteredMyRuns.subscribe((myRuns) => this.updateMyPane(myRuns));
    this.store.filteredFavoriteRuns.subscribe((newFavoriteRuns) => this.updateFavoritesPane(newFavoriteRuns));
    this.store.filteredSharedRuns.subscribe((sharedRuns) => this.updateSharedPane(sharedRuns));

    this.myRunsFilter.subscribe(() => {
      const filteredMyRuns = this.store.myRuns.value.filter((val) => {
        const startedAfter = (this.store.filteringOptions.startedAfter ? val.started > this.store.filteringOptions.startedAfter : true);
        const titleContainsText = !!val.options['title'] && val.options['title'].includes(this.store.filteringOptions.text);
        const descContainsText = !!val.options['description'] && val.options['description'].includes(this.store.filteringOptions.text);

        return (titleContainsText || descContainsText) && startedAfter;
      });
      this.store.filteredMyRuns.next(filteredMyRuns);
    });

    this.favRunsFilter.subscribe(() => {
      const filteredFavRuns = this.store.favoriteRuns.value.filter((val) => {
        const startedAfter = (this.store.filteringOptions.startedAfter ? val.started > this.store.filteringOptions.startedAfter : true);
        const titleContainsText = !!val.options['title'] && val.options['title'].includes(this.store.filteringOptions.text);
        const descContainsText = !!val.options['description'] && val.options['description'].includes(this.store.filteringOptions.text);

        return (titleContainsText || descContainsText) && startedAfter;
      });
      this.store.filteredFavoriteRuns.next(filteredFavRuns);
    });

    this.sharedRunsFilter.subscribe(() => {
      const filteredSharedRuns = this.store.sharedRuns.value.filter((val) => {
        const startedAfter = (this.store.filteringOptions.startedAfter ? val.started > this.store.filteringOptions.startedAfter : true);
        const isAuthored = this.store.filteringOptions.author ? val.author.id === this.store.filteringOptions.author.id: true;
        const titleContainsText = !!val.options['title'] && val.options['title'].includes(this.store.filteringOptions.text);
        const descContainsText = !!val.options['description'] && val.options['description'].includes(this.store.filteringOptions.text);

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
      this.store.filteredMyRuns.next(myRuns);
    });

    this.favRunsFetch.subscribe(async () => {
      const myRuns = (await historyUtils.pullRunsByName(this.func.name, [{author: grok.shell.user}], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.store.filteredFavoriteRuns.next(myRuns.filter((run) => run.options['isFavorite'] && !run.options['isImported']));
    });

    this.sharedRunsFetch.subscribe(async () => {
      const sharedRuns = (await historyUtils.pullRunsByName(this.func.name, [], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.store.filteredSharedRuns.next(sharedRuns.filter((run) => run.options['isShared']));
    });

    this.allRunsFetch.next();
    this.updateFilterPane();
  }

  get root() {
    return this.mainAcc.root;
  }
}
