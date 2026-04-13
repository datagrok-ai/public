# Анализ UI-архитектуры приложений HitTriage и MPO Profiles

Результаты анализа двух подходов к созданию приложений с UI в Datagrok.

## 1. HitTriage — многоприложенческая архитектура

### Общая информация

- **Пакет:** `packages/HitTriage`
- **Назначение:** оценка качества молекулярных и пептидных «хитов», управление кампаниями
- **4 приложения:** Hit Triage, Hit Design, PeptiHit, PepTriage
- **Объём:** ~34 TypeScript-файла, глубокая иерархия классов

### Регистрация приложений

В `src/package.ts` приложения регистрируются декоратором `@grok.decorators.app`:

```typescript
@grok.decorators.app({name: 'Hit Triage'})
hitTriageApp() { return new HitTriageApp(c).multiView; }

@grok.decorators.app({name: 'Hit Design'})
hitDesignApp() { return new HitDesignApp(c).multiView; }
```

Возвращается `DG.MultiView` — контейнер для переключения вкладок.

### Иерархия классов

```
HitAppBase<Template>              (абстрактный базовый класс)
├── HitTriageApp
│   └── PepTriageApp
├── HitDesignApp<Template>
│   └── PeptiHitApp

DG.ViewBase
└── HitBaseView<Template, App>
    ├── InfoView (Hit Triage)     — список кампаний и шаблонов
    ├── SubmitView (Hit Triage)   — отправка результатов
    ├── HitDesignInfoView         — управление кампаниями дизайна
    ├── HitDesignSubmitView       — отправка молекул
    ├── TilesView                 — прогресс-трекер по этапам
    └── Pep*InfoView              — пептидные варианты
```

### Инициализация

```
package.ts: @grok.decorators.app
  → new HitTriageApp(c)
    → constructor: создаёт InfoView, оборачивает в DG.MultiView
      → InfoView.constructor: checkCampaign() → init()  ← рендерит UI
```

### Структура главной страницы (InfoView.init)

```typescript
this.root.appendChild(ui.div([
  ui.divV([appHeader, sortingHeader]),   // ① Шапка с иконкой + заголовок "Continue Campaign"
  this.campaignsTableRoot,               // ② Таблица кампаний
  createNewCampaignHeader,               // ③ Заголовок "New Campaign"
  templatesDiv,                          // ④ Выбор шаблона
  campaignAccordionDiv,                  // ⑤ Форма создания кампании + кнопка START
]));
```

### Шапка с иконкой и описанием

Используется утилита `u2.appHeader()` из `@datagrok-libraries/utils`:

```typescript
u2.appHeader({
  iconPath: _package.webRoot + '/images/icons/hit-triage-icon.png',
  learnMoreUrl: '...',
  description: '- Configure your own workflow...\n- Calculate different molecular properties...',
});
```

### Таблица "Continue Campaign"

Колонки определены в `src/app/consts.ts`:

```typescript
const HTDefaultCampaignTableInfoGetters = {
  'Name': (info) => info.friendlyName ?? info.name,
  'Created': (info) => info.createDate,
  'Total': (info) => info.rowCount,
  'Selected': (info) => info.filteredRowCount,
  'Status': (info) => info.status,
  // Hit Design добавляет: 'Code', 'Author', 'Last Modified by', 'Molecules'
};
```

Рендеринг:

```typescript
ui.table(campaignsInfo, (info) => ([
  ui.link(info.name, () => this.setCampaign(info.name)),  // кликабельное имя
  ...shownColumns.map(col => getter(info)),               // остальные колонки
  shareIcon, deleteIcon,                                   // иконки управления
]));
```

Иконки управления над таблицей:
- `ui.iconFA('layer-group')` — группировка
- `ui.iconFA('eye')` — выбор видимых колонок
- `ui.iconFA('sync')` — обновление

### Форма "New Campaign" (аккордеон)

**Hit Triage** (`src/app/accordeons/new-campaign-accordeon.ts`):

```typescript
const dfInput = ui.input.table('Dataframe');              // загрузка файла
const dataSourceFunctionInput = ui.input.choice('Source'); // или выбор функции
const campaignPropsForm = ui.input.form(obj, props);      // кастомные поля из шаблона
const startButton = ui.bigButton('START', () => resolve({df, campaignProps}));
```

**Hit Design** (`src/app/accordeons/new-hit-design-campaign-accordeon.ts`):
- Проще: нет выбора источника данных
- Создаёт пустой DataFrame с колонками `molecule`, `stage`, `V-iD`

### Навигация

- **MultiView** с вкладками для переключения между видами
- **Breadcrumbs** (`src/app/utils.ts`): `ui.breadcrumbs(['Hit Triage', campaignName])`
- Клик по кампании → `setCampaign()` → загрузка JSON → `setTemplate()` → TableView

### Хранение данных

- JSON (метаданные кампании) + CSV (обогащённая таблица) в файловой системе
- Путь: `System.AppData/HitTriage/{appName}/campaigns/{id}/`
- Шаблоны: JSON-файлы в `files/`

### Сильные стороны

1. Элегантная мультиприложенческая архитектура через наследование
2. Template-driven формы — поля генерируются из JSON-шаблона
3. Система прав доступа (edit/view по группам)
4. Кеширование SMILES-конвертаций
5. Механизм блокировок для конкурентного редактирования

### Проблемы

1. Большие классы (HitTriageApp, HitDesignApp) — смешение ответственностей
2. `console.error()` без обратной связи пользователю
3. Утечки подписок RxJS (`.subscribe()` без cleanup)
4. Несколько `// @ts-ignore`
5. Хардкод путей
6. Отсутствие тестов

---

## 2. MPO Profiles — компактная CRUD-архитектура

### Общая информация

- **Пакет:** `packages/Chem` (подпапка `src/mpo/`)
- **Назначение:** создание и управление MPO-профилями для оценки молекул
- **Объём:** ~5 TypeScript-файлов, плоская структура

### Регистрация приложения

В `src/package.ts`:

```typescript
@grok.decorators.app({
  'name': 'MPO profiles',
  'meta': {browsePath: 'Chem', icon: 'images/mpo.png'},
})
static async mpoProfilesApp(path?: string): Promise<DG.View> {
  // Роутинг по URL:
  if (path && url.pathname.endsWith('/create-profile'))
    return new MpoProfileCreateView().tableView!;       // создание

  if (path && profileId)
    return new MpoProfileCreateView(profile).view;      // редактирование

  const infoView = new MpoProfilesView();               // список
  await infoView.render();
  return infoView.view;
}
```

### Ключевые файлы

| Файл | Назначение |
|------|-----------|
| `mpo-profiles-view.ts` | Главный вид со списком профилей |
| `mpo-profile-manager.ts` | Singleton для CRUD-операций |
| `mpo-profile-handler.ts` | `DG.ObjectHandler` — превью, контекстное меню |
| `mpo-create-profile.ts` | Вид создания/редактирования профиля |
| `utils.ts` | Типы, загрузка данных, breadcrumbs |

### Главный класс — MpoProfilesView

```typescript
export class MpoProfilesView {
  root = ui.divV([]);
  view: DG.View;

  constructor() {
    this.view = DG.View.fromRoot(this.root);   // View из обычного div
    this.view.name = 'MPO Profiles';
  }
}
```

Используется `DG.View.fromRoot()` — создание View из произвольного HTML-элемента (без наследования от DG.ViewBase).

### Метод render() — сборка страницы

```typescript
async render(): Promise<void> {
  ui.empty(this.root);
  this.root.append(
    this.buildHeader(),                              // ① Шапка с иконкой
    ui.h1('Manage Profiles'),                        // ② Заголовок
    this.tableContainer,                             // ③ Таблица профилей
    ui.divH([createButton, importButton]),            // ④ Кнопки внизу
  );
  await this.reloadProfiles();
  this.listenForChanges();
}
```

### Шапка (тот же паттерн, что и HitTriage)

```typescript
private buildHeader(): HTMLElement {
  return u2.appHeader({
    iconPath: `${_package.webRoot}/images/mpo.png`,
    description:
      '- Create and manage MPO profiles.\n' +
      '- Build profiles manually or train from labeled data.\n' +
      '- Initialize profiles from scratch or from an existing dataset.\n' +
      '- Full lifecycle support: create, edit, clone, and delete profiles.\n',
  });
}
```

### Таблица профилей

```typescript
private buildProfilesTable(): HTMLElement {
  const table = ui.table(
    MpoProfileManager.items,
    (profile) => [
      this.buildActionsButton(profile),         // кнопка ⋮ с popup-меню
      this.buildProfileLink(profile),           // кликабельное имя
      this.buildDescription(profile.description), // описание с обрезкой
    ],
    ['', 'Name', 'Description'],
  );
  table.classList.add('chem-mpo-profiles-table');
  return table;
}
```

### Кнопка ⋮ — контекстное меню

```typescript
private buildActionsButton(profile: MpoProfileInfo): HTMLElement {
  return ui.button('⋮', () => {
    ui.popupMenu()
      .item('Edit',     () => MpoProfileHandler.edit(profile))
      .item('Clone',    () => MpoProfileHandler.clone(profile))
      .item('Download', () => MpoProfileManager.download(profile))
      .separator()
      .item('Delete',   () => MpoProfileHandler.delete(profile));
  });
}
```

### Ссылка на профиль

```typescript
private buildProfileLink(profile: MpoProfileInfo): HTMLElement {
  const link = ui.link(profile.name, () => {
    this.previewedFileName = profile.fileName;   // клик → превью в Properties
  });
  link.addEventListener('dblclick', () => MpoProfileHandler.edit(profile));
  return ui.bind(profile, link);  // привязка объекта для context panel
}
```

`ui.bind(profile, link)` — при клике объект `profile` отображается в Context Panel через `MpoProfileHandler.renderProperties()`.

### Кнопки CREATE PROFILE / UPLOAD

```typescript
const createButton = ui.button(ui.h2('Create profile'), () => this.openCreateProfile());
const importButton = ui.button(ui.h2('Upload'), () => MpoProfileManager.upload());
```

**Upload** использует `DG.Utils.openFile()` с валидацией JSON:

```typescript
DG.Utils.openFile({accept: '.json', open: async (file) => {
  const parsed = JSON.parse(await file.text());
  if (!isDesirabilityProfile(parsed)) { grok.shell.warning('...'); return; }
  migrateProfile(parsed);
  await this.save(parsed, this.generateFileName(parsed.name));
}});
```

### Навигация

Используются **отдельные View** (не MultiView):

```typescript
static edit(profile: MpoProfileInfo): void {
  const view = new MpoProfileCreateView(editable, false, profile.fileName);
  grok.shell.v = grok.shell.addView(view.view);
  view.setupBreadcrumbs();
}
```

Breadcrumbs с иконкой Home:

```typescript
const breadcrumbs = ui.breadcrumbs(['Home', 'MPO Profiles', profileName]);
breadcrumbs.onPathClick.subscribe((path) => {
  if (clicked === 'MPO Profiles') {
    const listView = grok.shell.views.find(v => v.name === 'MPO Profiles');
    grok.shell.v = listView;
  }
});
```

### Реактивное обновление через custom events

```typescript
grok.events.onCustomEvent(MPO_PROFILE_CHANGED_EVENT)
  .subscribe(() => this.reloadProfiles());
grok.events.onCustomEvent(MPO_PROFILE_DELETED_EVENT)
  .subscribe((data) => { ... });
grok.events.onViewRemoving.subscribe((v) => {
  if (v.args.view.id === this.view.id) this.detach();
});
```

### MpoProfileManager — Singleton

```typescript
class MpoProfileManagerImpl {
  private profiles: MpoProfileInfo[] = [];
  async load(): Promise<MpoProfileInfo[]> { ... }
  async saveProfile(profile, fileName?): Promise<MpoSaveResult> { ... }
  confirmDelete(profile, onDeleted?): void { ... }
  upload(): void { ... }
  download(profile): void { ... }
}
export const MpoProfileManager = new MpoProfileManagerImpl();
```

### MpoProfileHandler — ObjectHandler

```typescript
export class MpoProfileHandler extends DG.ObjectHandler<MpoProfileInfo> {
  isApplicable(x: any): boolean { return isDesirabilityProfile(x); }
  renderTooltip(profile): HTMLElement { ... }
  renderProperties(profile): HTMLElement { ... }  // Context Panel
  static edit(profile): void { ... }
  static clone(profile): void { ... }
  static delete(profile): void { ... }
}
```

### Хранение данных

JSON-файлы в `System:AppData/Chem/mpo/`:

```typescript
const MPO_TEMPLATE_PATH = 'System:AppData/Chem/mpo';
await grok.dapi.files.list(MPO_TEMPLATE_PATH);
await grok.dapi.files.readAsText(`${MPO_TEMPLATE_PATH}/${file.name}`);
await grok.dapi.files.writeAsText(`${MPO_TEMPLATE_PATH}/${fileName}`, JSON.stringify(profile));
```

---

## 3. Сравнение архитектур

| Аспект | HitTriage | MPO Profiles |
|--------|-----------|-------------|
| **Сложность** | ~34 файла, глубокая иерархия | ~5 файлов, плоская структура |
| **View** | `DG.ViewBase` + `DG.MultiView` | `DG.View.fromRoot(div)` |
| **Навигация** | Вкладки MultiView + breadcrumbs | Отдельные View + breadcrumbs |
| **Данные** | JSON + CSV в файловой системе | JSON в `System:AppData/Chem/mpo/` |
| **Менеджер данных** | `_package.loadCampaigns()` | Singleton `MpoProfileManager` |
| **Context Panel** | Нет | `DG.ObjectHandler` с `renderProperties()` |
| **Обновление UI** | Ручной `init()` re-render | Custom events (`grok.events.onCustomEvent`) |
| **Шапка** | `u2.appHeader()` | `u2.appHeader()` — **одинаково** |
| **Таблица** | `ui.table()` | `ui.table()` — **одинаково** |
| **Формы** | Template-driven (`ui.input.form`) | Отдельный MpoProfileCreateView |

---

## 4. Общие паттерны UI Datagrok-приложений

### Компоненты из `datagrok-api/ui`

| Функция | Назначение |
|---------|-----------|
| `ui.div()`, `ui.divV()`, `ui.divH()` | Контейнеры (общий, вертикальный, горизонтальный) |
| `ui.h1()`, `ui.h2()` | Заголовки |
| `ui.table(data, renderer, headers)` | HTML-таблица из массива данных |
| `ui.link(text, onClick)` | Кликабельная ссылка |
| `ui.button()`, `ui.bigButton()` | Кнопки |
| `ui.iconFA(icon, onClick)` | Иконка FontAwesome |
| `ui.input.choice()`, `ui.input.string()` | Поля ввода |
| `ui.input.form(obj, props)` | Форма из DG.Property[] |
| `ui.dialog()` | Модальный диалог |
| `ui.popupMenu()` | Контекстное меню |
| `ui.breadcrumbs(paths)` | Хлебные крошки |
| `ui.bind(obj, element)` | Привязка объекта к элементу для Context Panel |
| `ui.empty(el)` | Очистка содержимого элемента |
| `ui.setUpdateIndicator(el, flag)` | Индикатор загрузки |

### Утилиты из `@datagrok-libraries/utils`

| Функция | Назначение |
|---------|-----------|
| `u2.appHeader({iconPath, description, learnMoreUrl})` | Шапка приложения с иконкой и буллетами |

### Создание View

```typescript
// Вариант 1: наследование (сложные приложения)
class MyView extends DG.ViewBase { ... }

// Вариант 2: из div (простые приложения)
const view = DG.View.fromRoot(ui.divV([...]));

// Вариант 3: MultiView (несколько вкладок)
const multiView = new DG.MultiView({
  viewFactories: { 'Tab1': () => view1, 'Tab2': () => view2 }
});
```

### Навигация

```typescript
// Breadcrumbs
const bc = ui.breadcrumbs(['App', 'Section', 'Item']);
bc.onPathClick.subscribe((path) => { ... });

// Переключение View
grok.shell.v = grok.shell.addView(newView);
grok.shell.addPreview(view);
```

### Хранение данных

```typescript
// Файловое хранилище пакета
await _package.files.readAsText('path/to/file.json');
await _package.files.writeAsText('path/to/file.json', data);
await _package.files.list('folder/', recursive);

// Или через dapi
await grok.dapi.files.readAsText('System:AppData/Package/path');
await grok.dapi.files.writeAsText('System:AppData/Package/path', data);
```

---

## 5. Model Hub — декларативный каталог на базе CardView

### Общая информация

- **Пакет:** `packages/Compute2` (регистрация) + `libraries/compute-utils/model-catalog/` (реализация)
- **Назначение:** каталог вычислительных моделей и скриптов с фильтрацией, поиском и иерархической навигацией
- **Объём:** ~3 файла основной логики + ObjectHandler

### Ключевые файлы

| Файл | Назначение |
|------|-----------|
| `packages/Compute2/src/package.ts` | Регистрация приложения, конфигурация |
| `libraries/compute-utils/model-catalog/src/model-catalog-view.ts` | Главный View (`extends DG.CustomCardView`) |
| `libraries/compute-utils/model-catalog/src/model-handler.ts` | ObjectHandler — рендеринг карточек, превью, меню |
| `libraries/compute-utils/model-catalog/src/model-catalog.ts` | Инициализация, события, tree browser, deep linking |
| `libraries/compute-utils/model-catalog/css/model-card.css` | Стили карточек и дерева |

### Регистрация приложения

В `packages/Compute2/src/package.ts`:

```typescript
@grok.decorators.app({
  browsePath: 'Compute',
  name: 'Model Hub',
})
static modelCatalog() {
  return startModelCatalog(modelCatalogOptions);
}
```

Конфигурация:

```typescript
const modelCatalogOptions = {
  _package,
  ViewClass: ModelCatalogView,
  segment: 'Modelhub',        // URI-сегмент для deep linking
  viewName: 'Model Hub',
  funcName: 'modelCatalog',
};
```

### ModelCatalogView — главный класс

Наследует `DG.CustomCardView` — платформенный класс для каталожных представлений с карточками:

```typescript
export class ModelCatalogView extends DG.CustomCardView {
  constructor(viewName: string) {
    super({dataSource: grok.dapi.functions, permanentFilter: MODEL_FILTER});

    this.root.classList.add('model-catalog-view');
    this.meta = new ModelHandler();           // кто рендерит карточки
    this.renderMode = DG.RENDER_MODE.BRIEF;   // компактный режим
    this.objectType = 'Func';                 // тип объектов
    this.showTree = true;                     // показать дерево слева
  }
}
```

**Платформа автоматически обеспечивает:**
- Загрузку данных с сервера с пагинацией
- Smart search (поиск)
- Галерею карточек
- Дерево навигации
- Фильтрацию по категориям

### Источник данных и фильтрация

**Постоянный фильтр** — серверный запрос к `grok.dapi.functions`:

```typescript
export const MODEL_FILTER = '(#model or options.role like "%model%")';
```

**Фильтры-категории** (левая панель с drill-down):

```typescript
this.categoryFilters = {
  'options.department': 'Department',
  'options.HL_process': 'HL_process',
  'options.process': 'Process',
  'tag': 'Tags',
};
```

**Быстрые текстовые фильтры:**

```typescript
this.filters = {
  'All': '',
  'Favorites': 'starredBy = @current',
  'Used by me': 'usedBy = @current',
};
```

**Иерархия для группировки (3 уровня):**

```typescript
this.hierarchy = [
  'options.department',
  'options.HL_process',
  'options.process',
];
```

### ModelHandler — рендеринг карточек

`model-handler.ts` наследует `DG.ObjectHandler` и переопределяет все методы рендеринга.

**Карточка в галерее** (`renderCard`):

```typescript
override renderCard(x: DG.Func): HTMLElement {
  const h = this.renderMarkup(x);           // иконка + имя
  h.classList.add('d4-link-label');
  return ui.bind(x, ui.h2(h));              // привязка объекта для контекста
}
```

**Разметка карточки** (`renderMarkup`):

```typescript
override renderMarkup(x: DG.Func): HTMLElement {
  const markup = ui.divH([]);
  const label = ui.span([
    this.renderIcon(x),                     // иконка модели
    ui.label(x.friendlyName),              // название
  ]);

  if (!hasMissingMandatoryGroups)
    markup.ondblclick = () => ModelHandler.openModel(x);  // дабл-клик → открыть
  else
    label.classList.add('d4-disabled');                     // или заблокировать

  markup.append(label);
  if (hasMissingMandatoryGroups)
    markup.append(warningIcon);                            // предупреждение о недостающих группах
}
```

**Иконка** (`renderIcon`) — приоритет разрешения:
1. `func.options['icon']` — HTTP URL
2. `package.getIconUrl()` — иконка пакета
3. Относительный путь от корня пакета
4. Иконка языка скрипта (`/images/entities/{language}.png`)
5. Fallback — SVG `'project'`

**Превью при клике** (`renderPreview`):

```typescript
override async renderPreview(x: DG.Func): Promise<DG.View> {
  // Если есть кастомный редактор — показать его
  if (x.options.editor === 'Compute2:TreeWizardEditor') {
    const call = x.prepare();
    const res = await DG.Func.byName(x.options.editor).prepare({call}).call();
    return res.getOutputParamValue() as DG.View;
  }

  // Иначе — стандартное превью
  v.append(ui.div([
    ui.h1(x.friendlyName),
    ui.divText(x.description),
    ui.bigButton('Launch', () => x.prepare().edit()),   // кнопка запуска
    ui.markdown(help ?? ''),                             // readme
  ]));
}
```

**Панель свойств** (`renderProperties`):

```typescript
override renderProperties(func: DG.Func) {
  const a = ui.accordion('ComputeModel');
  a.addTitle(ui.div([
    this.renderIcon(func),
    ui.label(func.friendlyName),
    ui.contextActions(func),        // стандартные действия
    ui.star(func.id),               // избранное
  ]));
  // + описание (markdown), теги, статус
}
```

**Тултип** (`renderTooltip`):

```typescript
override renderTooltip(x: DG.Func) {
  return ui.bind(x, ui.divV([
    ui.h2(ui.span([this.renderIcon(x), ui.label(x.friendlyName)])),
    ui.divV([ui.divText(x.description)]),
  ]));
}
```

### Открытие модели

**Двойной клик** на карточке:

```typescript
static openModel(x: DG.Func) {
  const fc = x.prepare();     // создать FuncCall
  fc.edit();                   // открыть редактор
}
```

**Выбор версии** (если модель поддерживает):

```typescript
static openVersionSelectDialog(func: DG.Func) {
  const versions = JSON.parse(func.options['versions']);
  const dlg = ui.dialog({title: 'Choose version'});
  const selInput = ui.input.choice('Select version', {items: versions});
  dlg.add(selInput).onOK(() => {
    func.prepare({params: {version: selInput.value}}).edit();
  });
  dlg.show();
}
```

### Deep linking по URL

`model-catalog.ts`:

```typescript
// URL: /apps/Compute/Modelhub/myModelName
async function handleInitialUri(segment: string) {
  const shortName = startPathSegments[catlogUriSegmentIdx + 1];
  const lst = await grok.dapi.functions
    .filter(`shortName = "${shortName}" and ${MODEL_FILTER}`).list();

  if (lst.length == 1) {
    const missingGroups = ModelHandler.getMissingGroups(func, userGroups);
    if (missingGroups?.length)
      grok.shell.error('User is not a part of...');
    else
      ModelHandler.openModel(func);
  }
}
```

### Tree Browser (боковая панель Browse)

Построение 3-уровневого дерева с lazy-loading:

```
Department
  └── HL_process
        └── Process → lazy load моделей при раскрытии
```

```typescript
processNode.onNodeExpanding.pipe(take(1)).subscribe(() => {
  const filteringRules = `(${MODEL_FILTER} and dept=X and hlProcess=Y and process=Z)`;
  processNode.loadSources(grok.dapi.functions.filter(filteringRules));
});
```

### Глобальные обработчики событий

```typescript
// Добавить панель «REST» с cURL/JS кодом для вызова модели через API
grok.events.onAccordionConstructed.subscribe((acc) => {
  if (acc.context?.type === 'script')
    acc.addPane('REST', () => renderRestPanel(acc.context));
});

// Привязать FuncCall к Model Hub для breadcrumbs
grok.functions.onBeforeRunAction.subscribe((fc) => {
  if (isModel(fc.func))
    ViewClass.bindModel(view, fc);
});
```

### REST-панель

Автоматическая генерация кода вызова модели через REST API:

```typescript
const curl = `curl --location --request POST '${apiUrl}/v1/func/${func.nqName}/run' \\
--header 'Authorization: ${auth}' \\
--header 'Content-Type: application/json' \\
--data-raw '${JSON.stringify(params)}'`;

const tabs = ui.tabControl({'CURL': ui.div([ui.divText(curl)]), 'JS': ...});
```

### Контроль доступа (Mandatory User Groups)

Хранится в `func.options['mandatoryUserGroups']` как JSON:

```typescript
const mandatoryUserGroups = JSON.parse(func.options['mandatoryUserGroups']);
// [{name: 'GroupName', help?: 'Description'}]
```

Если пользователь не в группе:
- Карточка заблокирована (`d4-disabled`)
- Иконка предупреждения с popover
- Ссылка «Request group membership» → POST `/api/groups/{id}/requests/{userId}`

### Метаданные модели (func.options)

| Ключ | Тип | Назначение |
|------|-----|-----------|
| `role` | string | Должен содержать «model» для попадания в каталог |
| `department` | string | Категория 1-го уровня |
| `HL_process` | string | Категория 2-го уровня |
| `process` | string | Категория 3-го уровня |
| `icon` | string | URL или относительный путь к иконке |
| `help` / `readme` | string | Путь к файлу помощи |
| `editor` | string | Кастомный редактор (напр. `Compute2:TreeWizardEditor`) |
| `mandatoryUserGroups` | JSON | Обязательные группы пользователей |
| `versions` | JSON | Доступные версии модели |
| `dev.status` | string | Статус разработки (бейдж) |

### Ribbon Menu

```typescript
async initMenu() {
  let help = this.ribbonMenu.group('Help');
  help.item('Developing Scripts', () => window.open('https://datagrok.ai/help/compute/scripting'));
  help.item('Developing Workflows', () => window.open('https://datagrok.ai/help/compute/workflows'));
  // + динамические пункты из ModelHub:getHelpItems
}
```

### CSS-стили

```css
.model-catalog-view                              /* Корневой контейнер */
  .d4-gallery-card.entity-model                  /* Карточка модели, 94px высота */
    i.grok-icon                                  /* Иконка 64x64px */
    .d4-link-label                               /* Имя + метка */
  .d4-tree-view-root                             /* Дерево иерархии */
    .d4-tree-view-group                          /* Уровень 1 (Department), 20px шрифт */
    > .d4-tree-view-group                        /* Уровень 2 (HL_process), 16px */
    > > .d4-tree-view-group                      /* Уровень 3 (Process), 14px */
```

---

## 6. Итоговое сравнение трёх подходов

| Аспект | HitTriage | MPO Profiles | Model Hub |
|--------|-----------|-------------|-----------|
| **Базовый класс** | `DG.ViewBase` + `DG.MultiView` | `DG.View.fromRoot(div)` | `DG.CustomCardView` |
| **Парадигма** | Ручная сборка DOM | Ручная сборка DOM | Декларативная конфигурация |
| **Список элементов** | `ui.table()` | `ui.table()` | Галерея карточек (встроенная) |
| **Фильтрация** | Ручная (localStorage) | Нет | Декларативная (`categoryFilters`, `filters`) |
| **Поиск** | Нет | Нет | Smart search (встроенный) |
| **Дерево** | Tree browser (ручной) | Нет | `showTree: true` + `hierarchy` |
| **Рендеринг элемента** | Inline в `ui.table()` | Методы `build*()` | `DG.ObjectHandler` (override `render*()`) |
| **Данные** | JSON/CSV файлы | JSON файлы | `grok.dapi.functions` (сервер) |
| **Шапка** | `u2.appHeader()` | `u2.appHeader()` | Ribbon menu |
| **Context Panel** | Нет | `DG.ObjectHandler` | `DG.ObjectHandler` |
| **Deep linking** | URL params (`campaignId`) | URL params (`profileId`) | URI segment (`/Modelhub/{shortName}`) |
| **Навигация** | MultiView + breadcrumbs | Отдельные View + breadcrumbs | `parentCall` hierarchy + breadcrumbs |
| **Контроль доступа** | Группы edit/view | Нет | Mandatory user groups |
| **Объём кода** | ~34 файла | ~5 файлов | ~3 файла (+handler) |

### Когда какой подход использовать

- **HitTriage (ручной DOM + MultiView):** сложные workflow с несколькими этапами, формами, шаблонами и CRUD над файлами
- **MPO Profiles (ручной DOM + View.fromRoot):** простые CRUD-приложения с небольшим числом сущностей
- **Model Hub (CustomCardView):** каталоги/галереи с фильтрацией, поиском и иерархической навигацией по серверным данным
