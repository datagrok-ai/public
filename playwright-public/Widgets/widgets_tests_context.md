# Home page Widgets — контекст проекта (handoff)

Сопровождение к `home_widgets_manual_tests.md`. Содержит весь контекст для продолжения
работы в новом чате: цель, окружение, что такое Home page Widgets, как устроен UI,
снятые вживую селекторы, подводные камни автоматизации и план переноса в Playwright.

---

## 1. Цель и текущий этап

- Пишем **мануальные тест-кейсы для Home page Widgets** (виджеты на лендинг-странице Datagrok), покрытие через UI.
- Конечная цель — **автотесты на Playwright** (по максимуму через UI, API — в крайнем случае).
- Сейчас этап: готова мануалка (`home_widgets_manual_tests.md`). Дальше — ревью, затем перенос в Playwright.
- Скоуп согласован: **framework виджетов + глубокие per-widget кейсы** (Spotlight, Community, Usage/Reports).

## 2. Окружение

- **Целевой инстанс:** `https://dev.datagrok.ai/`.
- Авторизация через `global-setup` (storageState `e2e/.auth.json`), креды в `playwright-tests/.env`
  (`DATAGROK_LOGIN=opavlenko+playwright@datagrok.ai`).
- **Важно про права (обновлено 2026-06-08):** тестовому пользователю `opavlenko+playwright`
  **добавлена роль администратора**. Поэтому виджеты `Usage` и `Reports` (`canView: Developers,Administrators`)
  теперь **видны** под основным аккаунтом → позитивные кейсы `Widgets-Usage-*` / `Widgets-Reports-*` / `Widgets-Perm-02`
  выполняются под ним.
- Для **negative**-кейса `Widgets-Perm-01` (виджеты скрыты без прав) нужен non-admin контекст —
  вторичный аккаунт `DATAGROK_SHARING_LOGIN` (`opavlenko+pwsharing`, storageState `.auth-sharing.json`).
  Перед использованием убедиться, что у него нет ролей Developers/Administrators.

## 3. Что такое Home page Widgets

- **Home Page** — особый лендинг-вью, открывается при логине или по иконке **Home** (Top Menu) / **House** в шапке Browse.
- Рендерится пакетом **PowerPack** → `welcome-view.ts` (`welcomeView()`, функция с `meta.autostartImmediate`).
- Сами виджеты — это функции с ролью `dashboard` и `returnType: widget`, разбросанные по пакетам.
  Welcome view находит их через `DG.Func.find({meta: {role: 'dashboard'}, returnType: 'widget'})`
  и сортирует по опции `order`.

**Реально зарегистрированные виджеты на dev (снято вживую 2026-06-08):**

| friendlyName | function name            | package        | order | showName | canView                  | Виден тест-юзеру |
|--------------|--------------------------|----------------|-------|----------|--------------------------|------------------|
| Spotlight    | `activityDashboardWidget`| PowerPack      | -1    | false    | —                        | ✅ да            |
| Community    | `communityWidget`        | PowerPack      | 6     | —        | —                        | ✅ да            |
| Usage        | `usageWidget`            | UsageAnalysis  | —     | —        | Developers,Administrators| ✅ да (юзер — admin)|
| Reports      | `reportsWidget`          | UsageAnalysis  | —     | —        | Developers,Administrators| ✅ да (юзер — admin)|

> Набор зависит от установленных плагинов и прав. На других инстансах (public/release)
> могут быть дополнительные виджеты (Recent projects, Learn, Why Datagrok? и т.д. — см. вики-скриншот).
> При расширении набора — перегенерировать таблицу.

## 4. Устройство UI (снято из кода + вживую)

**Иерархия DOM:**
- Корень вью: `.power-pack-welcome-view`
- Строка поиска: `.power-search-search-everywhere-input` (placeholder `Search everywhere. Try "aspirin" or "7JZK"`), хост `.d4-search-bar`
- Панель виджетов: `.power-pack-widgets-panel` → содержит:
  - `.power-pack-widgets-host` — контейнер виджетов
  - ссылку **«Customize widgets...»** (`.ui-link`)
- Хост результатов поиска: `.power-pack-search-host` (скрыт, пока строка поиска пуста)

**Один виджет:** `.power-pack-widget-host[widget-title="<friendlyName>"]`
- заголовок: `.d4-dialog-header` → `.d4-dialog-title` (текст). Если `showName=false` — заголовок скрыт классом `.d4-dialog-header-hidden`.
- контент: `.power-pack-widget-content`
- **close-иконка ("Remove"):** `i.grok-icon.grok-font-icon-close`. Скрыта до hover (`ui.tools.setHoverVisibility`), тултип «Remove» появляется на hover.
- порядок задаётся inline `style.order` из опции `order`.

**Поведение:**
1. Виджеты сортируются по `order` (Spotlight `-1` → первым, Community `6` → ниже).
2. **Hover** по виджету → показывается close-иконка. **Клик** по ней → виджет убирается из DOM **и** в user settings пишется `ignored=true`.
3. Ссылка **«Customize widgets...»** → открывает Context Panel с формой bool-тогглов (по одному на виджет, label = friendlyName). Снятие галки → виджет прячется + `ignored=true`; возврат галки → виджет показывается + `ignored=false`. Здесь же можно вернуть ранее закрытый виджет.
4. **Home-виджеты управляются только через «Customize widgets...» / close-иконку** (хранилище `widgets`).
   Раздел **Settings → Panels** к ним отношения не имеет — он управляет **context-панелями**
   (**Hidden Panels** / **Panel Order**, записи вида `type: Column, semType: Molecule`).
   (Ранее вики `browse.md` ошибочно отправляла в Settings → Panels — **исправлено** на реальное поведение 2026-06-08.)
5. **Поиск:** ввод в строку (debounce 500 мс) прячет `.power-pack-widgets-panel`, показывает `.power-pack-search-host` с результатами, URL → `search?q=<query>`. Очистка строки возвращает виджеты, URL → `search`.

**Persistence:** настройки хранятся в `grok.userSettings` под ключом `widgets`; значение —
объект `{ [functionName]: '{"ignored":bool}' }` (значение — JSON-строка). Ключи: `activityDashboardWidget`,
`communityWidget`, `usageWidget`, `reportsWidget`.

## 5. Per-widget детали (снято вживую)

**Spotlight** (`activityDashboardWidget`, без заголовка) — богатый мульти-таб виджет. Табы (`.d4-tab-header`):
- **Workspace** — закреплённые (pinned) сущности; «Select a pinned item on the left to see its controls here».
- **Spotlight** — Shared with me / Recent сущности.
- **Favorites** — избранное.
- **Notifications** (с бейджем-счётчиком непрочитанных) — список уведомлений, ссылка **«Mark all as read»**.
- **My Activity** — активность пользователя.
- **Learn** — подвкладки **VIDEO / WIKI / DEMO / TUTORIALS** с темами (Meetings, Develop, Cheminformatics, Visualize, Explore).
- Внизу виджета — **«💡 Demo of the day: <name>»** (ссылка, открывает демо).

**Community** (`communityWidget`, заголовок «Community») — список внешних ссылок на community.datagrok.ai
(Platform Releases, Plugin releases, Macromolecules/Cheminformatics updates, и т.д.). Ссылки открываются во внешней вкладке.

**Usage** (`usageWidget`, UsageAnalysis) — графики **Users** и **Errors** (line charts на canvas),
агрегации (count/avg/min/max/med/sum/stdev), селектор Time unit, блок **System** со ссылками
(Datlas, Datlas DB, Grok Connect, Jupyter, Grok Spawner), ссылка **«Open Usage Analysis»**.

**Reports** (`reportsWidget`, UsageAnalysis) — ссылка **«Open Reports»** + список недавних error-reports.
⚠️ **Снято вживую 2026-06-08:** виджет показывал 4 записи «NullError: method not found: 'clientWidth' on null».
Это недавние **auto-error-reports**, сгенерированные исходной ошибкой `NullError: method not found: 'clientWidth' on null`
(вероятно — рендер виджета в zero-width контейнере). Сам виджет работает (листает отчёты), но ошибку
стоит триажить как кандидата на баг.

## 6. Формат тест-кейсов

- Markdown, текстовое описание (как в `e2e/browse/`).
- Каждый кейс: **Цель → Тип → Предусловия → Шаги (действие + ожидаемая реакция) → Итоговая проверка → Постусловия/очистка**.
- **ID:** `Widgets-<Область>-<NN>` (совместимо с Test Track).
- **Тип:** `smoke` / `functional` / `regression` / `negative` → теги Playwright (`@smoke` и т.д.).
- Шаги атомарные: действие → `await`, ожидание → `expect`. Предусловия → `beforeEach`, постусловия → `afterEach`/teardown.

## 7. Подводные камни автоматизации (ВАЖНО)

- **Async-сохранение настроек.** `saveSettings()` зовёт `grok.userSettings.addAll(...)` — это **асинхронный** серверный вызов, который **не** await-ится. При проверке persistence нельзя сразу делать reload: сохранение `ignored` не успеет долететь до сервера и кейс будет флаки. **Снято вживую:** быстрый reload после close показал виджет снова (`ignored:false`); reload с паузой ~5 с — корректно скрыл (`ignored:true`). В Playwright: после close/toggle ждать сетевой idle или появление/исчезновение хоста + небольшую паузу перед reload.
- **canView-гейтинг.** `opavlenko+playwright` теперь admin → Usage/Reports видны (позитивные кейсы идут под ним). Для negative-кейса `Widgets-Perm-01` нужен non-admin контекст (`.auth-sharing.json`).
- **Очистка обязательна.** Кейсы, которые закрывают/прячут виджеты, **меняют серверные настройки пользователя** и обязаны восстанавливать состояние в `afterEach`/`finally` (вернуть галку в Customize → `ignored:false`). Иначе следующий прогон стартует с «дырявым» Home page.
- **Hover-видимость close-иконки.** Синтетический `mouseenter` через `dispatchEvent` НЕ включает видимость иконки (нужны реальные события). Playwright `locator.hover()` — сработает. Сам клик по close-иконке работает и без hover (можно `click({force:true})` как фолбэк).
- **Спотлайт-табы** — `.d4-tab-header` с текстом таба; Notifications-таб имеет суффикс-счётчик («Notifications1»), фильтровать по `startsWith`.
- **Внешние ссылки Community** — открываются на `community.datagrok.ai` в новой вкладке; в тесте проверять `href`, а не реальный переход (внешний ресурс).
- **assertNoErrors** — переиспользовать `watchErrors`/`expectNoErrors` из `e2e/browse/helpers.ts` (или вынести в общий helpers для Widgets).

## 8. Снятые селекторы (сводка для Playwright)

| Элемент | Селектор |
|---|---|
| Корень Home view | `.power-pack-welcome-view` |
| Панель виджетов | `.power-pack-widgets-panel` |
| Контейнер виджетов | `.power-pack-widgets-host` |
| Один виджет | `.power-pack-widget-host[widget-title="<Name>"]` |
| Заголовок виджета | `.d4-dialog-title` |
| Контент виджета | `.power-pack-widget-content` |
| Close-иконка («Remove») | `i.grok-icon.grok-font-icon-close` (внутри хоста виджета) |
| Ссылка Customize | текст `Customize widgets...` (`.ui-link`) |
| Тоггл в Customize | `.ui-input-root.ui-input-bool` (label = friendlyName) |
| Строка поиска | `.power-search-search-everywhere-input` |
| Хост результатов поиска | `.power-pack-search-host` |
| Табы Spotlight | `.power-pack-widget-host[widget-title="Spotlight"] .d4-tab-header` |
| Home (house) иконка | `[name="icon-home"]` |
| Хранилище настроек | `grok.userSettings.get('widgets')` |

## 9. План дальнейших шагов

1. Ревью мануалки (`home_widgets_manual_tests.md`).
2. Снять оставшиеся точечные селекторы при первом прогоне (Settings→Panels, под-табы Learn).
3. Перенос в Playwright: сперва `@smoke` (Display, базовый рендер), затем `@functional` (Controls, Customize, Search, Spotlight), затем `@regression`/`@negative` (persistence, canView).
4. Общий teardown: восстановить все `widgets`-настройки в `ignored:false`.

## 10. Открытые вопросы

- **[РЕШЕНО]** Аккаунт для Usage/Reports — `opavlenko+playwright` получил роль admin, кейсы идут под ним. Для negative нужен non-admin (`opavlenko+pwsharing`).
- **[РЕШЕНО]** Settings → Panels — это управление **context-панелями**, а не home-виджетами. Вики `browse.md` исправлена на реальное поведение; отдельный тест-кейс убран как бессмысленный.
- **NullError 'clientWidth' on null** в виджете Reports — завести/найти тикет на исходную ошибку (auto-error-reports засоряют список Reports).
- Подтвердить, что у `opavlenko+pwsharing` нет ролей Developers/Administrators (для надёжности `Widgets-Perm-01`).
- Стабильность набора виджетов на dev (Recent projects/Learn/Why Datagrok? из вики-скриншота на dev отсутствуют) — зафиксировать ожидаемый набор в фикстуре.
