/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import wu from 'wu';
import {BehaviorSubject, Subject} from 'rxjs';
import {historyUtils} from '../../history-utils';
import '../css/history-panel.css';

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
    tags?: string[]
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
    return this.myRuns.filter((run) => {
      const startedAfter = (this.filteringOptions.startedAfter ? run.started > this.filteringOptions.startedAfter : true);
      const titleContainsText = (!this.filteringOptions.text.length) ? true : (!!run.options['title']) ? run.options['title'].includes(this.filteringOptions.text): false;
      const descContainsText = (!this.filteringOptions.text.length) ? true : (!!run.options['description']) ? run.options['description'].includes(this.filteringOptions.text): false;


      const hasTags = (!this.filteringOptions.tags || this.filteringOptions.tags.length === 0) ? true :
        (!!run.options['tags']) ? this.filteringOptions.tags.every((searchedTag) => run.options['tags'].includes(searchedTag)): false;

      return (titleContainsText || descContainsText) && startedAfter && hasTags;
    });
  }

  get filteredFavoriteRuns() {
    return this.favoriteRuns.filter((run) => {
      const startedAfter = (this.filteringOptions.startedAfter ? run.started > this.filteringOptions.startedAfter : true);
      const titleContainsText = (!this.filteringOptions.text.length) ? true : (!!run.options['title']) ? run.options['title'].includes(this.filteringOptions.text): false;
      const descContainsText = (!this.filteringOptions.text.length) ? true : (!!run.options['description']) ? run.options['description'].includes(this.filteringOptions.text): false;
      const hasTags = (!this.filteringOptions.tags || this.filteringOptions.tags.length === 0) ? true :
        (!!run.options['tags']) ? this.filteringOptions.tags.every((searchedTag) => run.options['tags'].includes(searchedTag)): false;

      return (titleContainsText || descContainsText) && startedAfter && hasTags;
    });
  }

  get filteredSharedRuns() {
    return this.sharedRuns.filter((run) => {
      const isAuthored = this.filteringOptions.author ? run.author.id === this.filteringOptions.author.id: true;
      const startedAfter = (this.filteringOptions.startedAfter ? run.started > this.filteringOptions.startedAfter : true);
      const titleContainsText = (!this.filteringOptions.text.length) ? true : (!!run.options['title']) ? run.options['title'].includes(this.filteringOptions.text): false;
      const descContainsText = (!this.filteringOptions.text.length) ? true : (!!run.options['description']) ? run.options['description'].includes(this.filteringOptions.text): false;
      const hasTags = (!this.filteringOptions.tags || this.filteringOptions.tags.length === 0) ? true :
        (!!run.options['tags']) ? this.filteringOptions.tags.every((searchedTag) => run.options['tags'].includes(searchedTag)): false;


      return (titleContainsText || descContainsText) && startedAfter && isAuthored && hasTags;
    });
  }
}

const MY_PANE_LABEL = 'HISTORY' as const;
const FAVORITES_LABEL = 'FAVORITES' as const;
const SHARED_LABEL = 'SHARED' as const;

export class HistoryPanel {
  // Emitted when FuncCall should is chosen. Contains FuncCall ID
  public onRunChosen = new Subject<string>();

  // Emitted when FuncCalls are called for comparison. Contains FuncCalls' IDs
  public onComparison = new Subject<string[]>();

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

  private selectedCallsSet = new Set<DG.FuncCall>;

  public allRunsFetch = new Subject<true>();

  historyTab = ui.div();
  favTab = ui.div();
  sharedTab = ui.div();
  filterPane = ui.div();
  tabs = ui.tabControl({
    [MY_PANE_LABEL]: ui.box(this.historyTab),
    [FAVORITES_LABEL]: ui.box(this.favTab),
    [SHARED_LABEL]: ui.box(this.sharedTab),
  });
  actionsSection = this.buildActionsSection();

  private _root = ui.divV([
    this.tabs.root,
  ], {style: {width: '100%', height: '100%'}});

  myCards = [] as HTMLElement[];
  favoriteCards = [] as HTMLElement[];
  sharedCards = [] as HTMLElement[];

  buildActionsSection() {
    return ui.divH([
      ui.span([`Selected: ${this.selectedCallsSet.size}`], {style: {'align-self': 'center'}}),
      ui.divH([
        (() => {
          const t = ui.iconFA('exchange', async () => {
            this.onComparison.next([...wu(this.selectedCallsSet.keys()).map((selected) => selected.id)]);
          }, 'Compare selected runs');
          t.style.margin = '5px';
          if (this.selectedCallsSet.size < 2)
            t.classList.add('hp-disabled');
          return t;
        })(),
        (() => {
          const t = ui.iconFA('trash-alt', () => this.showDeleteRunDialog(this.selectedCallsSet), 'Delete selected runs');
          t.style.margin = '5px';
          if (this.selectedCallsSet.size === 0)
            t.classList.add('hp-disabled');
          return t;
        })(),
        ...this.selectedCallsSet.size === 0 ? [(() => {
          const t = ui.iconFA('square', () => {
            switch (this.tabs.currentPane.name) {
            case MY_PANE_LABEL:
              this.store.myRuns.forEach((run) => this.selectedCallsSet.add(run));
              this.myRunsFilter.next();
              break;
            case FAVORITES_LABEL:
              this.store.favoriteRuns.forEach((run) => this.selectedCallsSet.add(run));
              this.favRunsFilter.next();
              break;
            case SHARED_LABEL:
              this.store.sharedRuns.forEach((run) => this.selectedCallsSet.add(run));
              this.sharedRunsFilter.next();
              break;
            }
            this.updateActionsSection();
          }, 'Select all'); t.style.margin = '5px'; return t;
        })()]: [
          (() => {
            let fullListCount = 0;
            switch (this.tabs.currentPane.name) {
            case MY_PANE_LABEL:
              fullListCount = this.store.myRuns.length;
              break;
            case FAVORITES_LABEL:
              fullListCount = this.store.favoriteRuns.length;
              break;
            case SHARED_LABEL:
              fullListCount = this.store.sharedRuns.length;
              break;
            }

            const iconType = this.selectedCallsSet.size === fullListCount? 'check-square': 'minus-square';
            const t = ui.iconFA(iconType, () => {
              switch (this.tabs.currentPane.name) {
              case MY_PANE_LABEL:
                this.store.myRuns.forEach((run) => this.selectedCallsSet.delete(run));
                this.myRunsFilter.next();
                break;
              case FAVORITES_LABEL:
                this.store.favoriteRuns.forEach((run) => this.selectedCallsSet.delete(run));
                this.favRunsFilter.next();
                break;
              case SHARED_LABEL:
                this.store.sharedRuns.forEach((run) => this.selectedCallsSet.delete(run));
                this.sharedRunsFilter.next();
                break;
              }
              this.updateActionsSection();
            }, 'Unselect all');

            t.style.margin = '5px';
            return t;
          })(),
        ],
      ]),
    ], {style: {
      'justify-content': 'space-between',
      'padding': '0 12px',
    }});
  }

  buildFilterPane(currentTabName: string) {
    return ui.wait(async () => {
      const filteringText = new Subject();

      const textInput = ui.stringInput('Search', '', (v: string) => filteringText.next(v));
      DG.debounce(filteringText.asObservable(), 600).subscribe(() => {
        this.store.filteringOptions.text = textInput.stringValue;
        this.myRunsFilter.next();
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

      let tagInput = ui.choiceInput<string>('Tag', 'Choose tag to filter', ['Choose tag to filter']);

      const filterTagEditor = DG.TagEditor.create();
      const dummyInput = ui.stringInput(' ', '');
      dummyInput.input.replaceWith(filterTagEditor.root);
      dummyInput.root.style.display = 'none';

      filterTagEditor.onChanged(() => {
        //@ts-ignore
        this.store.filteringOptions.tags = filterTagEditor.tags.filter((tag) => !!tag);

        if (this.store.filteringOptions.tags!.length === 0)
          dummyInput.root.style.display = 'none';
        else
          dummyInput.root.style.removeProperty('display');

        this.myRunsFilter.next();
        this.favRunsFilter.next();
        this.sharedRunsFilter.next();
      });

      this.store.allRuns.subscribe(() => {
        const tags = this.store.allRuns.value
          .map((run) => run.options['tags'] as string[])
          .reduce((acc, runTags) => {
            if (!!runTags) runTags.forEach((runTag) => {if (!acc.includes(runTag)) acc.push(runTag);});
            return acc;
          }, [] as string[]);

        const newTagInput = ui.choiceInput<string>('Tag', 'Choose tag to filter', ['Choose tag to filter', ...tags]);
        tagInput.root.replaceWith(newTagInput.root);
        tagInput = newTagInput;

        tagInput.onChanged(() => {
          if (
            tagInput.value &&
            tagInput.value !== 'Choose tag to filter' &&
            //@ts-ignore
            !filterTagEditor.tags.includes(tagInput.value)
          )
            filterTagEditor.addTag(tagInput.value);
        });
      });

      dateInput.addPatternMenu('datetime');
      const form = ui.divV([
        textInput,
        dateInput,
        authorInput,
        tagInput,
        dummyInput,
      ], 'ui-form ui-form-wide ui-form-left');

      if (currentTabName === MY_PANE_LABEL) {
        dateInput.root.style.removeProperty('display');
        authorInput.root.style.display = 'none';
      }
      if (currentTabName === FAVORITES_LABEL) {
        dateInput.root.style.removeProperty('display');
        authorInput.root.style.display = 'none';
      }
      if (currentTabName === SHARED_LABEL) {
        dateInput.root.style.removeProperty('display');
        authorInput.root.style.removeProperty('display');
      }

      return form;
    });
  };

  myPaneFilter = this.buildFilterPane(MY_PANE_LABEL) as HTMLElement;
  favoritesPaneFilter =this.buildFilterPane(FAVORITES_LABEL) as HTMLElement;
  sharedPaneFilter = this.buildFilterPane(SHARED_LABEL)as HTMLElement;

  showDeleteRunDialog(funcCalls: Set<DG.FuncCall>) {
    ui.dialog({title: `Delete  ${funcCalls.size > 1 ? 'runs': 'run'}`})
      .add(ui.divText(`The deleted ${funcCalls.size > 1 ? 'runs': 'run'} is impossible to restore. Are you sure?`))
      .onOK(async () => {
        funcCalls.forEach(async (funcCall) => {
          this.beforeRunDeleted.next(funcCall.id);
          this.selectedCallsSet.delete(funcCall);
          await historyUtils.deleteRun(funcCall);
          this.afterRunDeleted.next(funcCall.id);
        });
      })
      .show({center: true});
  };

  updateActionsSection() {
    const newActionsSection = this.buildActionsSection();
    this.actionsSection.replaceWith(newActionsSection);
    this.actionsSection = newActionsSection;
  }

  updateSharedPane(sharedRuns: DG.FuncCall[]) {
    const sharedCards = (sharedRuns.length > 0) ? this.renderSharedCards(sharedRuns) : ui.divText('No runs are marked as shared', 'hi-no-elements-label');
    const sharedTab = ui.divV([
      this.sharedPaneFilter,
      ui.element('div', 'splitbar-horizontal'),
      sharedCards,
    ]);
    this.sharedTab.replaceWith(sharedTab);
    this.sharedTab = sharedTab;
  };

  updateFavoritesPane(favoriteRuns: DG.FuncCall[]) {
    const favCards = (favoriteRuns.length > 0) ? this.renderFavoriteCards(favoriteRuns) : ui.divText('No runs are marked as favorites', 'hi-no-elements-label');
    const favTab = ui.divV([
      this.favoritesPaneFilter,
      ui.element('div', 'splitbar-horizontal'),
      favCards,
    ]);
    this.favTab.replaceWith(favTab);
    this.favTab = favTab;
  };

  updateMyPane(myRuns: DG.FuncCall[]) {
    const myCards = (myRuns.length > 0) ? this.renderHistoryCards(myRuns) : ui.divText('No runs are found in history', 'hi-no-elements-label');
    const newTab = ui.divV([
      this.myPaneFilter,
      this.actionsSection,
      ui.element('div', 'splitbar-horizontal'),
      myCards,
    ]);
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

    const addToCompare = ui.iconFA('square', (ev) => {
      ev.stopPropagation();
      this.selectedCallsSet.add(funcCall);

      addToCompare.style.display = 'none';
      removeFromCompare.style.removeProperty('display');
      this.updateActionsSection();
    }, 'Select this run');
    addToCompare.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    const removeFromCompare = ui.iconFA('check-square', (ev) => {
      ev.stopPropagation();
      this.selectedCallsSet.delete(funcCall);
      removeFromCompare.style.display = 'none';

      this.updateActionsSection();
      addToCompare.style.removeProperty('display');
    }, 'Unselect this run');
    removeFromCompare.classList.add('hp-funccall-card-icon');

    if (!this.selectedCallsSet.has(funcCall)) {
      addToCompare.style.removeProperty('display');
      removeFromCompare.style.display = 'none';
    } else {
      addToCompare.style.display = 'none';
      removeFromCompare.style.removeProperty('display');
    }

    const editIcon = ui.iconFA('edit', (ev) => {
      ev.stopPropagation();
      this.showEditDialog(funcCall);
    }, 'Edit selected run');
    editIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    const deleteIcon = ui.iconFA('trash-alt', async (ev) => {
      ev.stopPropagation();
      const tempSet = new Set<DG.FuncCall>();
      tempSet.add(funcCall);
      this.showDeleteRunDialog(tempSet);
    }, 'Delete selected runs');
    deleteIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    const favoritedIcon = ui.iconFA('star', null, 'Unfavorite the run');
    favoritedIcon.classList.add('fas', 'hp-funccall-card-def-icon');

    const sharedIcon = ui.iconFA('share-alt', null, 'Add to shared');
    sharedIcon.classList.add('fas', 'hp-funccall-card-def-icon');
    sharedIcon.style.color = 'var(--blue-1)';

    const dateStarted = new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'});

    const card = ui.divH([
      ui.divH([
        icon,
        ui.divV([
          cardLabel,
          ui.span([dateStarted]),
          ...(funcCall.options['description']) ? [ui.divText(funcCall.options['description'], 'description')]: [],
          ...(funcCall.options['tags'] && funcCall.options['tags'].length > 0) ?
            [ui.div(funcCall.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')))]:[],
        ], 'hp-card-content'),
      ]),
      ui.divH([
        editIcon, deleteIcon,
        ...(funcCall.options['isFavorite']) ? [favoritedIcon] : [],
        ...(funcCall.options['isShared']) ? [sharedIcon]: [],
        addToCompare, removeFromCompare,
      ]),
    ], 'hp-funccall-card');

    card.addEventListener('click', async () => {
      this.onRunChosen.next(funcCall.id);

      const allCards = [...this.myCards, ...this.favoriteCards, ...this.sharedCards];
      allCards.forEach((c) => c.classList.remove('clicked')); ``;

      card.classList.add('clicked');
    });
    ui.tooltip.bind(card, () => ui.tableFromMap({
      Author: grok.shell.user.toMarkup(),
      Date: dateStarted,
      ...(funcCall.options['title']) ? {'Title': funcCall.options['title']}:{},
      ...(funcCall.options['description']) ? {'Description': funcCall.options['description']}:{},
      ...(funcCall.options['tags'] && funcCall.options['tags'].length > 0) ?
        {'Tags': ui.div(funcCall.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')), 'd4-tag-editor')}:{},
    }));
    return card;
  }

  showEditDialog(funcCall: DG.FuncCall) {
    const dlg = ui.dialog({title: 'Edit run'});

    let title = funcCall.options['title'] ?? '';
    let description = funcCall.options['description'] ?? '';
    const titleInput = ui.stringInput('Title', title, (s: string) => title = s);
    let isShared = funcCall.options['isShared'] ?? false;
    let isFavorite = funcCall.options['isFavorite'] ?? false;

    const tagsLine = DG.TagEditor.create();
    (funcCall.options['tags'] || []).forEach((tag: string) => {
      tagsLine.addTag(tag);
    });

    const dummyInput = ui.stringInput(' ', '');
    dummyInput.input.replaceWith(tagsLine.root);

    const tagInput = ui.stringInput('Tag', '');
    const addNewTag = () => {
      if (tagInput.value === '' ||
        // @ts-ignore
        tagsLine.tags.includes(tagInput.value)
      )
        return;

      tagsLine.addTag(tagInput.value);
      tagInput.value = '';
    };

    tagInput.addOptions(ui.iconFA('plus', addNewTag, 'Add this tag'));

    tagInput.input.onkeydown = async (ev) => {
      if (ev.key == 'Enter') {
        ev.stopPropagation();
        addNewTag();
      }
    };

    dlg.add(ui.form([
      titleInput,
      ui.stringInput('Description', description, (s: string) => {description = s;}),
      ui.boolInput('Favorites', isFavorite, (b: boolean) => isFavorite = b),
      ui.boolInput('Share', isShared, (b: boolean) => isShared = b),
      tagInput,
      dummyInput,
    ]))
      .onOK(async () => {
        funcCall.options['title'] = title !== '' ? title : null;
        funcCall.options['description'] = description !== '' ? description : null;
        funcCall.options['tags'] = [...tagsLine.tags];

        if (!!isShared === !!funcCall.options['isShared'] && !!isFavorite === !!funcCall.options['isFavorite']) {
          ui.setUpdateIndicator(this.tabs.root, true);
          await historyUtils.saveRun(funcCall);

          const editedRun = this.store.allRuns.value.find((call) => call.id === funcCall.id)!;
          editedRun.options['title'] = title !== '' ? title : null;
          editedRun.options['description'] = description !== '' ? description : null;
          editedRun.options['tags'] = [...funcCall.options['tags']];

          this.store.allRuns.next(this.store.allRuns.value);
          ui.setUpdateIndicator(this.tabs.root, false);
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
      .show({center: true, width: 500});
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
    private func: DG.Func,
  ) {
    this.tabs.root.style.width = '100%';
    this.tabs.root.style.height = '100%';

    this.tabs.onTabChanged.subscribe(() => {
      switch (this.tabs.currentPane.name) {
      case MY_PANE_LABEL:
        this.myPaneFilter.append(this.actionsSection);
        break;
      case FAVORITES_LABEL:
        this.favoritesPaneFilter.append(this.actionsSection);
        break;
      case SHARED_LABEL:
        this.sharedPaneFilter.append(this.actionsSection);
        break;
      }
    });

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

    this.afterRunAddToFavorites.subscribe((added) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === added.id)!;
      editedRun.options['title'] = added.options['title'];
      editedRun.options['description'] = added.options['description'];
      editedRun.options['tags'] = added.options['tags'];
      editedRun.options['isFavorite'] = true;
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunAddToShared.subscribe((added) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === added.id)!;
      editedRun.options['title'] = added.options['title'];
      editedRun.options['description'] = added.options['description'];
      editedRun.options['tags'] = added.options['tags'];
      editedRun.options['isShared'] = true;
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunRemoveFromFavorites.subscribe((removed) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === removed.id)!;
      editedRun.options['title'] = removed.options['title'];
      editedRun.options['description'] = removed.options['description'];
      editedRun.options['tags'] = removed.options['tags'];
      editedRun.options['isFavorite'] = false;
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunRemoveFromShared.subscribe((removed) => {
      const editedRun = this.store.allRuns.value.find((call) => call.id === removed.id)!;
      editedRun.options['title'] = removed.options['title'];
      editedRun.options['description'] = removed.options['description'];
      editedRun.options['tags'] = removed.options['tags'];
      editedRun.options['isShared'] = false;
      this.store.allRuns.next(this.store.allRuns.value);
    });

    this.afterRunDeleted.subscribe((id) => {
      this.store.allRuns.next(this.store.allRuns.value.filter((call) => call.id !== id));
      this.updateActionsSection();
    });

    this.allRunsFetch.next();
  }

  private async addRunToFavorites(callToFavorite: DG.FuncCall): Promise<DG.FuncCall> {
    callToFavorite.options['isFavorite'] = true;

    this.beforeRunAddToFavorites.next(callToFavorite);
    const savedFavorite = await grok.dapi.functions.calls.allPackageVersions().save(callToFavorite);
    this.afterRunAddToFavorites.next(savedFavorite);
    return savedFavorite;
  }

  private async removeRunFromFavorites(callToUnfavorite: DG.FuncCall): Promise<DG.FuncCall> {
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
