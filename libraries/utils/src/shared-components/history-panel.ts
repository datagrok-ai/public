/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import wu from 'wu';
import {BehaviorSubject, Subject} from 'rxjs';
import {historyUtils} from '../history-utils';
import '../../css/history-panel.css';

export const defaultUsersIds = {
  'Test': 'ca1e672e-e3be-40e0-b79b-d2c68e68d380',
  'Admin': '878c42b0-9a50-11e6-c537-6bf8e9ab02ee',
  'System': '3e32c5fa-ac9c-4d39-8b4b-4db3e576b3c3',
};

export const defaultGroupsIds = {
  'All users': 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6',
  'Developers': 'ba9cd191-9a50-11e6-9cc9-910bf827f0ab',
  'Need to create': '00000000-0000-0000-0000-000000000000',
  'Test': 'ca1e672e-e3be-40e0-b79b-8546b8bf62e6',
  'Admin': 'a4b45840-9a50-11e6-c537-6bf8e9ab02ee',
  'System': 'a4b45840-ac9c-4d39-8b4b-4db3e576b3c3',
  'Administrators': '1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5',
};

class HistoryPanelStore {
  filteringOptions = {text: ''} as {
    text: string,
    author?: DG.User,
    startedAfter?: dayjs.Dayjs,
  };

  allRuns = new BehaviorSubject<DG.FuncCall[]>([]);

  get myRuns() {
    return this.allRuns.value.filter((run) => run.author.id === grok.shell.user.id);
  }

  get favoriteRuns() {
    return this.allRuns.value.filter((run) => run.options['isFavorite'] && !run.options['isImported']);
  }

  get sharedRuns() {
    return this.allRuns.value.filter((run) => run.options['isShared']);
  }

  get filteredMyRuns() {
    return this.myRuns.filter((val) => (this.filteringOptions.startedAfter ? val.started > this.filteringOptions.startedAfter : true));
  }

  get filteredFavoriteRuns() {
    return this.favoriteRuns.filter((val) => {
      const startedAfter = (this.filteringOptions.startedAfter ? val.started > this.filteringOptions.startedAfter : true);
      const titleContainsText = (!this.filteringOptions.text.length) ? true : (!!val.options['title']) ? val.options['title'].includes(this.filteringOptions.text): false;
      const descContainsText = (!this.filteringOptions.text.length) ? true : (!!val.options['description']) ? val.options['description'].includes(this.filteringOptions.text): false;

      return (titleContainsText || descContainsText) && startedAfter;
    });
  }

  get filteredSharedRuns() {
    return this.sharedRuns.filter((val) => {
      const isAuthored = this.filteringOptions.author ? val.author.id === this.filteringOptions.author.id: true;
      const startedAfter = (this.filteringOptions.startedAfter ? val.started > this.filteringOptions.startedAfter : true);
      const titleContainsText = (!this.filteringOptions.text.length) ? true : (!!val.options['title']) ? val.options['title'].includes(this.filteringOptions.text): false;
      const descContainsText = (!this.filteringOptions.text.length) ? true : (!!val.options['description']) ? val.options['description'].includes(this.filteringOptions.text): false;

      return (titleContainsText || descContainsText) && startedAfter && isAuthored;
    });
  }
}

export class HistoryPanel {
  // Emitted when FuncCall should is chosen. Contains FuncCall ID
  public onRunChosen = new Subject<string>();

  // Emitted when FuncCall is added to favorites
  public beforeRunAddToFavorites = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is added to shared
  public beforeRunAddToShared = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is removed from favorites
  public beforeRunRemoveFromFavorites = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is removed from shared
  public beforeRunRemoveFromShared = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is deleted
  public beforeRunDeleted = new Subject<string>();

  // Emitted when FuncCall is added to favorites
  public afterRunAddToFavorites = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is added to shared
  public afterRunAddToShared = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is removed from favorites
  public afterRunRemoveFromFavorites = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is removed from shared
  public afterRunRemoveFromShared = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is deleted
  public afterRunDeleted = new Subject<string>();

  private store = new HistoryPanelStore();

  private myRunsFilter = new Subject<true>();
  private favRunsFilter = new Subject<true>();
  private sharedRunsFilter = new Subject<true>();

  public allRunsFetch = new Subject<true>();

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
  ], {style: {width: '100%', height: '100%'}});

  myCards = [] as HTMLElement[];
  favoriteCards = [] as HTMLElement[];
  sharedCards = [] as HTMLElement[];

  updateFilterPane() {
    const buildFilterPane = () => {
      return ui.wait(async () => {
        const filteringText = new Subject();

        const textInput = ui.stringInput('Search', '', (v: string) => filteringText.next(v));
        textInput.root.style.display = 'none';
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
        authorInput.root.style.display = 'none';

        dateInput.addPatternMenu('datetime');
        const form = ui.divV([
          textInput,
          dateInput,
          authorInput,
        ], 'ui-form ui-form-wide ui-form-left');
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

  showDeleteRunDialog(funcCall: DG.FuncCall) {
    ui.dialog({title: 'Delete run'})
      .add(ui.divText('The deleted run is impossible to restore. Are you sure?'))
      .onOK(async () => {
        this.beforeRunDeleted.next(funcCall.id);
        await historyUtils.deleteRun(funcCall);
        this.afterRunDeleted.next(funcCall.id);
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

  renderCard(funcCall: DG.FuncCall) {
    const icon = funcCall.author.picture as HTMLElement;
    icon.style.width = '25px';
    icon.style.height = '25px';
    icon.style.fontSize = '20px';
    const cardLabel = ui.label(funcCall.options['title'] ?? funcCall.author.friendlyName, {style: {'color': 'var(--blue-1)'}});
    ui.bind(funcCall.author, icon);

    const editIcon = ui.iconFA('edit', (ev) => {
      ev.stopPropagation();
      this.showEditDialog(funcCall);
    });
    editIcon.classList.add('cv-funccall-card-hover-icon');

    const deleteIcon = ui.iconFA('trash-alt', async (ev) => {
      ev.stopPropagation();
      this.showDeleteRunDialog(funcCall);
    }, 'Delete the run');
    deleteIcon.classList.add('cv-funccall-card-hover-icon');

    const favoritedIcon = ui.iconFA('star', null, 'Unfavorite the run');
    favoritedIcon.classList.add('fas', 'cv-funccall-card-def-icon');

    const sharedIcon = ui.iconFA('share-alt', null, 'Add to shared');
    sharedIcon.classList.add('fas', 'cv-funccall-card-def-icon');
    sharedIcon.style.color = 'var(--blue-1)';

    const dateStarted = new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'});

    const card = ui.divH([
      ui.divH([
        icon,
        ui.divV([
          cardLabel,
          ui.span([dateStarted]),
          ...(funcCall.options['description']) ? [ui.divText(funcCall.options['description'], 'description')]: [],
        ], 'cv-card-content'),
      ]),
      ui.divH([
        editIcon, deleteIcon,
        ...(funcCall.options['isFavorite']) ? [favoritedIcon] : [],
        ...(funcCall.options['isShared']) ? [sharedIcon]: [],
      ])
    ], 'cv-funccall-card');

    card.addEventListener('click', async () => {
      this.onRunChosen.next(funcCall.id);

      const allCards = [...this.myCards, ...this.favoriteCards, ...this.sharedCards];
      allCards.forEach((c) => c.classList.remove('clicked'));

      card.classList.add('clicked');
    });
    ui.tooltip.bind(card, () => ui.tableFromMap({
      Author: grok.shell.user.toMarkup(),
      Date: dateStarted,
      ...(funcCall.options['title']) ? {'Title': funcCall.options['title']}:{},
      ...(funcCall.options['description']) ? {'Description': funcCall.options['description']}:{},
    }));
    return card;
  }

  showEditDialog(funcCall: DG.FuncCall) {
    const dlg = ui.dialog({title: 'Edit run'});

    let title = funcCall.options['title'] ?? '';
    let description = funcCall.options['description'] ?? '';
    const titleInput = ui.stringInput('Title', title, (s: string) => {
      title = s;
      if (s.length === 0) {
        titleInput.setTooltip('Title cannot be empty');
        setTimeout(() => titleInput.input.classList.add('d4-invalid'), 100);
        dlg.getButton('OK').disabled = true;
      } else {
        titleInput.setTooltip('');
        setTimeout(() => titleInput.input.classList.remove('d4-invalid'), 100);
        dlg.getButton('OK').disabled = false;
      }
    });
    let isShared = funcCall.options['isShared'] ?? false;
    let isFavorite = funcCall.options['isFavorite'] ?? false;

    dlg.add(ui.form([
      titleInput,
      ui.stringInput('Description', description, (s: string) => { description = s; }),
      ui.boolInput('Favorites', isFavorite, (b: boolean) => isFavorite = b),
      ui.boolInput('Share', isShared, (b: boolean) => isShared = b),
    ]))
      .onOK(async () => {
        funcCall = await historyUtils.loadRun(funcCall.id);
        funcCall.options['title'] = title;
        funcCall.options['description'] = description;

        if (isShared === funcCall.options['isShared'] && isFavorite === funcCall.options['isFavorite']) {
          await historyUtils.saveRun(funcCall);

          const editedRun = this.store.allRuns.value.find((call) => call.id === funcCall.id);
          editedRun!.options['title'] = funcCall.options['title'];
          editedRun!.options['description'] = funcCall.options['description'];

          this.store.allRuns.next(this.store.allRuns.value);
          return;
        }

        if (isShared && !funcCall.options['isShared']) {
          funcCall.options['isShared'] = isShared;
          ui.setUpdateIndicator(this.tabs.root, true);
          await this.addRunToShared(funcCall);
          ui.setUpdateIndicator(this.tabs.root, false);
        }

        if (!isShared && funcCall.options['isShared']) {
          funcCall.options['isShared'] = isShared;
          ui.setUpdateIndicator(this.tabs.root, true);
          await this.removeRunFromShared(funcCall);
          ui.setUpdateIndicator(this.tabs.root, false);
        }

        if (isFavorite && !funcCall.options['isFavorite']) {
          funcCall.options['isFavorite'] = isFavorite;
          ui.setUpdateIndicator(this.tabs.root, true);
          await this.addRunToFavorites(funcCall);
          ui.setUpdateIndicator(this.tabs.root, false);
        }

        if (!isFavorite && funcCall.options['isFavorite']) {
          funcCall.options['isFavorite'] = isFavorite;
          ui.setUpdateIndicator(this.tabs.root, true);
          await this.removeRunFromFavorites(funcCall);
          ui.setUpdateIndicator(this.tabs.root, false);
        }
      })
      .show({center: true});
    titleInput.fireChanged();
  }

  renderFavoriteCards(funcCalls: DG.FuncCall[]) {
    this.favoriteCards = funcCalls.map((funcCall) => this.renderCard(funcCall));
    return ui.divV(this.favoriteCards);
  };

  renderSharedCards(funcCalls: DG.FuncCall[]) {
    this.sharedCards = funcCalls.map((funcCall) => this.renderCard(funcCall));
    return ui.divV(this.sharedCards);
  };

  renderHistoryCards(funcCalls: DG.FuncCall[]) {
    this.myCards = funcCalls.map((funcCall) => this.renderCard(funcCall));
    return ui.divV(this.myCards, {style: {height: '100%'}});
  };

  constructor(
    private func: DG.Func
  ) {
    this.tabs.root.style.width = '100%';
    this.tabs.root.style.height = '100%';

    this.myRunsFilter.subscribe(() => this.updateMyPane(this.store.filteredMyRuns));
    this.favRunsFilter.subscribe(() => this.updateFavoritesPane(this.store.filteredFavoriteRuns));
    this.sharedRunsFilter.subscribe(() => this.updateSharedPane(this.store.filteredSharedRuns));

    this.store.allRuns.subscribe(() => {
      this.updateMyPane(this.store.filteredMyRuns);
      this.updateFavoritesPane(this.store.filteredFavoriteRuns);
      this.updateSharedPane(this.store.filteredSharedRuns);
    });

    this.allRunsFetch.subscribe(async () => {
      ui.setUpdateIndicator(this.root, true);
      const allRuns = (await historyUtils.pullRunsByName(this.func.name, [], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.store.allRuns.next(allRuns);
      ui.setUpdateIndicator(this.root, false);
    });

    const clearMetadata = (funcCall: DG.FuncCall) => {
      if (!funcCall.options['isFavorite'] && !funcCall.options['isShared']) {
        funcCall.options['title'] = null;
        funcCall.options['description'] = null;
      }
    };

    this.afterRunAddToFavorites.subscribe((added) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === added.id);
      editedRun!.options['title'] = added.options['title'];
      editedRun!.options['description'] = added.options['description'];
      editedRun!.options['isFavorite'] = true;
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunAddToShared.subscribe((added) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === added.id);
      editedRun!.options['title'] = added.options['title'];
      editedRun!.options['description'] = added.options['description'];
      editedRun!.options['isShared'] = true;
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunRemoveFromFavorites.subscribe((removed) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === removed.id)!;
      editedRun.options['isFavorite'] = false;
      clearMetadata(editedRun);
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunRemoveFromShared.subscribe((removed) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === removed.id)!;
      editedRun.options['isShared'] = false;
      clearMetadata(editedRun);
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunDeleted.subscribe((id) => {
      this.store.allRuns.next(this.store.allRuns.value.filter((call) => call.id !== id));
    });

    this.allRunsFetch.next();
    this.updateFilterPane();
  }

  private async addRunToFavorites(callToFavorite: DG.FuncCall): Promise<DG.FuncCall> {
    callToFavorite.options['isFavorite'] = true;

    this.beforeRunAddToFavorites.next(callToFavorite);
    const savedFavorite = await grok.dapi.functions.calls.allPackageVersions().save(callToFavorite);
    this.afterRunAddToFavorites.next(savedFavorite);
    return savedFavorite;
  }

  private async removeRunFromFavorites(callToUnfavorite: DG.FuncCall): Promise<DG.FuncCall> {
    callToUnfavorite.options['title'] = null;
    callToUnfavorite.options['description'] = null;
    callToUnfavorite.options['isFavorite'] = false;

    this.beforeRunRemoveFromFavorites.next(callToUnfavorite);
    const favoriteSave = await grok.dapi.functions.calls.allPackageVersions().save(callToUnfavorite);
    this.afterRunRemoveFromFavorites.next(favoriteSave);
    return favoriteSave;
  }

  private async addRunToShared(callToShare: DG.FuncCall): Promise<DG.FuncCall> {
    callToShare.options['isShared'] = true;

    this.beforeRunAddToShared.next(callToShare);

    const allGroup = await grok.dapi.groups.find(defaultGroupsIds['All users']);

    const dfOutputs = wu(callToShare.outputParams.values() as DG.FuncCallParam[])
      .filter((output) => output.property.propertyType === DG.TYPE.DATA_FRAME);

    for (const output of dfOutputs) {
      const df = callToShare.outputs[output.name] as DG.DataFrame;
      await grok.dapi.permissions.grant(df.getTableInfo(), allGroup, false);
    }

    const dfInputs = wu(callToShare.inputParams.values() as DG.FuncCallParam[])
      .filter((input) => input.property.propertyType === DG.TYPE.DATA_FRAME);
    for (const input of dfInputs) {
      const df = callToShare.inputs[input.name] as DG.DataFrame;
      await grok.dapi.permissions.grant(df.getTableInfo(), allGroup, false);
    }

    const savedShared = await grok.dapi.functions.calls.allPackageVersions().save(callToShare);

    this.afterRunAddToShared.next(savedShared);
    return savedShared;
  }

  private async removeRunFromShared(callToUnshare: DG.FuncCall): Promise<DG.FuncCall> {
    callToUnshare.options['title'] = null;
    callToUnshare.options['description'] = null;
    callToUnshare.options['isShared'] = false;

    this.beforeRunRemoveFromShared.next(callToUnshare);
    const savedShared = await grok.dapi.functions.calls.allPackageVersions().save(callToUnshare);
    this.afterRunRemoveFromShared.next(savedShared);
    return savedShared;
  }

  public get root() {
    return this._root;
  }
}
