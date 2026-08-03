# Browse selectors — reference

Каталог логических имён селекторов, используемых в `playwright-tests/e2e/browse/*.test.ts`.

**Источник истины:** `selectors.ts`. Этот файл — короткое оглавление для ревью мануалки и Test Track.

**Условные обозначения:**
- ✅ значение снято на dev и подтверждено
- ⏳ TODO — требует доразведки (см. в `selectors.ts` маркер `TODO`)

Дата разведки: **2026-06-03** (dev.datagrok.ai, user `opavlenko+playwright@datagrok.ai`).

---

## Sidebar (крайняя левая полоса)

| Имя                              | Селектор                                       | Статус |
| -------------------------------- | ---------------------------------------------- | ------ |
| `SIDEBAR`                        | `.d4-tab-header-stripe.layout-sidebar.vertical` | ✅     |
| `SIDEBAR_BROWSE_ICON`            | `.d4-tab-header[name="Browse"]`                | ✅     |
| `SIDEBAR_BROWSE_ICON_BY_DATA`    | `[data-source="tab-pane-Browse"]`              | ✅     |
| `SIDEBAR_FAVORITES_ICON`         | `.d4-tab-header[name="Favorites"]`             | ✅     |
| `SIDEBAR_TAB_SELECTED`           | `.d4-tab-header.selected`                      | ✅     |
| `SIDEBAR_VIEW_BADGE`             | бейдж-счётчик однотипных view                  | ⏳     |

---

## Browse panel — шапка и иконки

| Имя                              | Селектор                                                        | Статус |
| -------------------------------- | --------------------------------------------------------------- | ------ |
| `BROWSE_PANEL_TAB_PANE`          | `[name="tab-pane-Browse"]`                                      | ✅     |
| `BROWSE_HEADER`                  | `.panel-titlebar.disable-selection.grok-browse-header`          | ✅     |
| `BROWSE_HEADER_ICONS`            | `<BROWSE_HEADER> i.grok-icon`                                   | ✅     |
| `BROWSE_HEADER_HOME`             | `[name="icon-home"]` (`fa-home`)                                | ✅     |
| `BROWSE_HEADER_IMPORT_FILE`      | `[name="icon-folder-open"]` (`fa-folder-open`)                  | ✅     |
| `BROWSE_HEADER_IMPORT_TEXT`      | `[name="icon-file-alt"]` (`fa-file-alt`)                        | ✅     |
| `BROWSE_HEADER_REFRESH`          | `[name="icon-sync"]` (`fa-sync`)                                | ✅     |
| `BROWSE_HEADER_COLLAPSE_ALL`     | `[name="icon-chevron-double-up"]` (`fa-chevron-double-up`)      | ✅     |
| `BROWSE_HEADER_LOCATE`           | `[name="icon-dot-circle"]` (`fa-dot-circle`)                    | ✅     |
| `BROWSE_HEADER_CLOSE_PANEL`      | `.grok-font-icon-close` (внутри `BROWSE_HEADER`)                | ✅     |

> Итого в шапке Browse **6 функциональных иконок + 1 close** (= 7 кликабельных, не считая контейнер `.grok-browse-icons`).

---

## Browse tree

| Имя                              | Селектор                                                  | Статус |
| -------------------------------- | --------------------------------------------------------- | ------ |
| `TREE_NODE_GROUP_LABEL`          | `.d4-tree-view-group-label`                               | ✅     |
| `TREE_NODE_ITEM_LABEL`           | `.d4-tree-view-item-label`                                | ✅     |
| `TREE_NODE_LABEL_ANY`            | оба класса выше                                            | ✅     |
| `TREE_NODE_CONTAINER`            | `.d4-tree-view-node`                                       | ✅     |
| `TREE_EXPAND_ARROW`              | `.d4-tree-view-tri`                                        | ✅     |
| `TREE_EXPAND_ARROW_EXPANDED`     | `.d4-tree-view-tri.d4-tree-view-tri-expanded`              | ✅     |
| `TREE_NODE_DROP`                 | `.d4-tree-view-node.d4-tree-drop` (drop-зона для DnD)      | ✅     |

**Динамические локаторы:**

| Функция                    | Описание                                       | Статус |
| -------------------------- | ---------------------------------------------- | ------ |
| `treeNodeByName(p, name)`  | Узел дерева по точному имени                   | ✅     |
| `treeGroupByName(p, name)` | Узел-группа (с детьми)                          | ✅     |
| `treeItemByName(p, name)`  | Узел-лист (без детей)                          | ✅     |
| `treeNodeContainer(p, n)`  | Контейнер `.d4-tree-view-node` для узла        | ✅     |
| `treeExpandArrow(p, name)` | Стрелка expand/collapse слева                  | ✅     |
| `treeNodeChildren(p, n)`   | Дочерние узлы развёрнутой группы               | ✅     |

> Когда меню/aria уже инициализированы, есть также прямой `[name="tree-node-<Name>"]` (например `tree-node-My-files`), но он надёжно срабатывает только после первого взаимодействия с узлом.

---

## Context menu (right-click)

| Имя                              | Селектор                                                   | Статус |
| -------------------------------- | ---------------------------------------------------------- | ------ |
| `CONTEXT_MENU`                   | `.d4-menu-item-container.d4-vert-menu.d4-menu-popup`       | ✅     |
| `CONTEXT_MENU_ITEM`              | `.d4-menu-item.d4-menu-item-vert`                          | ✅     |
| `CONTEXT_MENU_ITEM_LABEL`        | `.d4-menu-item-label`                                      | ✅     |
| `contextMenuItemByName(p, name)` | по атрибуту `d4-name="..."`                                | ✅     |
| `contextMenuItem(p, text)`       | по тексту пункта                                            | ✅     |

**Точные тексты пунктов (важно: регистр и многоточия):**

| Имя                              | Текст пункта         |
| -------------------------------- | -------------------- |
| `CONTEXT_MENU_BROWSE`            | `Browse`             |
| `CONTEXT_MENU_ADD_FAVORITES`     | **`Add to favorites`** (нижний регистр `f`!) |
| `CONTEXT_MENU_SHARE`             | `Share...`           |
| `CONTEXT_MENU_RENAME`            | `Rename...`          |
| `CONTEXT_MENU_DELETE`            | `Delete...`          |
| `CONTEXT_MENU_CLONE`             | `Clone...`           |
| `CONTEXT_MENU_EDIT`              | `Edit...`            |
| `CONTEXT_MENU_TEST_CONNECTION`   | `Test connection`    |
| `CONTEXT_MENU_CLEAR_CACHE`       | `Clear cache`        |
| `CONTEXT_MENU_DOWNLOAD`          | ⏳ уточнить           |
| `CONTEXT_MENU_BROWSE_SCHEMA`     | ⏳ уточнить           |
| `CONTEXT_MENU_RUN`               | ⏳ уточнить           |

---

## Context Panel (правая панель свойств)

| Имя                              | Селектор                                                                     | Статус |
| -------------------------------- | ---------------------------------------------------------------------------- | ------ |
| `CONTEXT_PANEL`                  | `.grok-prop-panel`                                                            | ✅     |
| `CONTEXT_PANEL_BASE`             | `.panel-base.splitter-container-vertical`                                     | ✅     |
| `CONTEXT_PANEL_HEADER`           | `.panel-titlebar.disable-selection` (первый внутри `CONTEXT_PANEL_BASE`)      | ✅     |
| `CONTEXT_PANEL_BACK`             | `[aria-label="Back"]` (`fa-chevron-left`, `name="icon-chevron-left"`)         | ✅     |
| `CONTEXT_PANEL_FORWARD`          | `[aria-label="Forward"]` (`fa-chevron-right`)                                 | ✅     |
| `CONTEXT_PANEL_CLONE`            | `[aria-label="Clone and detach"]` (`fa-clone`) — кейс CtxPanel-05 удалён       | ✅     |
| `CONTEXT_PANEL_COLLAPSE_ALL`     | `[aria-label="Collapse all"]` (`fa-chevron-square-up`)                        | ✅     |
| `CONTEXT_PANEL_EXPAND_ALL`       | `[aria-label="Expand all"]` (`fa-chevron-square-down`)                        | ✅     |
| `CONTEXT_PANEL_HELP`             | `[aria-label="Show Help"]` (`grok-font-icon-help`)                            | ✅     |
| `CONTEXT_PANEL_STAR`             | `[aria-label="Favorites"]` (`fa-star`, `name="icon-star"`)                    | ✅     |
| `CONTEXT_PANEL_INNER`            | `.grok-prop-panel .grok-entity-prop-panel`                                    | ✅     |
| `CONTEXT_PANEL_ACCORDION`        | `.grok-prop-panel .d4-accordion`                                              | ✅     |
| `CONTEXT_PANEL_ACCORDION_PANE`   | `.d4-accordion-pane`                                                          | ✅     |
| `infoPaneByName(p, name)`        | `.d4-accordion-pane[name="pane-<Name>"]` (например `pane-Search`)             | ✅     |

> Звезда (`CONTEXT_PANEL_STAR`) — в шапке Context Panel. Состояние «активная/оранжевая» — снимется при первой проверке Browse-Fav-02 (вероятно через дополнительный класс или цвет).

---

## View / tabs / top ribbon

| Имя                              | Селектор                                                  | Статус |
| -------------------------------- | --------------------------------------------------------- | ------ |
| `RIBBON`                         | `.d4-ribbon`                                              | ✅     |
| `VIEW_TAB`                       | `.tab-handle.disable-selection`                            | ✅     |
| `VIEW_TAB_SELECTED`              | `.tab-handle.tab-handle-selected`                          | ✅     |
| `viewTabHandle(p, name)`         | `[name="view-handle: <name>"]`                            | ✅     |
| `VIEW_TAB_TEXT`                  | `.tab-handle-text`                                         | ✅     |
| `VIEW_TAB_CLOSE`                 | `.tab-handle-close-button[name="Close"]`                   | ✅     |
| `TOP_MENU`                       | `.grok-main-menu`                                          | ✅     |
| `TOP_MENU_ROOT`                  | `.grok-main-menu[d4-name="top"]`                           | ✅     |
| `topMenuSection(p, name)`        | `.grok-main-menu [d4-name="<name>"]`                       | ✅     |
| `topMenuItem(p, [path])`         | `[name="div-<part1>---<part2>---..."]`                     | ✅     |
| `VIEW_CLOSE`                     | `[aria-label="Close view"]`                                | ✅     |
| `VIEW_CLOSE_BY_NAME(name)`       | `[aria-label="Close <name>"]`                              | ✅     |
| `TOOLBOX_PIN_TOOLBOX`            | `[aria-label="Pin toolbox"]` (`fa-thumbtack`)              | ✅     |
| `VIEW_TAB_PIN`                   | `[aria-label*="Pin"]` для табов view                       | ⏳     |
| `VIEW_HAMBURGER` / `Split right/down` | Не найдено в `simpleMode=true`; вероятно в windows mode | ⏳     |

---

## Status bar (внизу экрана; виден когда открыт view)

| Имя                              | Селектор                                                | Статус |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `STATUS_BAR`                     | `.layout-status-bar`                                    | ✅     |
| `STATUS_BAR_TASK`                | `.layout-status-bar .d4-task-bar`                       | ✅     |
| `STATUS_BAR_VIEW_PANEL`          | `.layout-status-bar .d4-view-status-panel` (Rows/Cols/Selected) | ✅ |
| `STATUS_BAR_GLOBAL_PANEL`        | `.layout-status-bar .d4-global-status-panel`            | ✅     |
| `STATUS_BAR_AI_HELPER`           | `.layout-status-bar i.fa-user-robot`                    | ✅     |
| `STATUS_BAR_TERMINAL`            | `.layout-status-bar i.fa-terminal`                      | ✅     |
| `STATUS_BAR_MODE_MAXIMIZE`       | `.layout-status-bar i.fa-window-maximize`               | ⏳ подтвердить mapping → Tabbed/Simple/Presentation |
| `STATUS_BAR_MODE_PRESENTATION`   | `.layout-status-bar i.fa-presentation`                  | ⏳     |

> На иконках status bar нет `name` / `title` / `aria-label`. Точное соответствие иконок ↔ Tabbed/Simple/Presentation уточню при тесте Browse-Modes-01.

---

## Home Page / global search

| Имя                              | Селектор                                                | Статус |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `HOME_VIEW_HANDLE`               | `[name="view-handle: Home"]`                            | ✅     |
| `HOME_GLOBAL_SEARCH_INPUT`       | `input.power-search-search-everywhere-input` (placeholder: `Search everywhere. Try "aspirin" or "7JZK"`) | ✅ |
| `HOME_WIDGETS`                   | контейнер виджетов                                      | ⏳     |
| `HOME_SEARCH_RESULTS`            | область результатов поиска                              | ⏳     |

---

## Filtering / lists

| Имя                              | Селектор                                                | Статус |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `FILTER_PANEL`                   | панель фильтров в галерее (Apps/Files и т.п.)           | ⏳ воспроизвести в галерее |
| `FILTER_TOGGLE`                  | `[aria-label="Toggle filters"]`                         | ✅     |
| `filterProperty(p, name)`        | свойство фильтра по имени                               | ⏳     |
| `LIST_SEARCH_INPUT`              | `.d4-search-input.ui-input-editor`                      | ✅     |
| `LIST_EMPTY_STATE`               | пустое состояние                                        | ⏳     |

---

## Notifications / errors / balloons

| Имя                              | Селектор                                                | Статус |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `BALLOON_CONTAINER`              | `.d4-balloon-container`                                 | ✅     |
| `ERROR_BALLOON`                  | внутри `BALLOON_CONTAINER`, подтип error                 | ⏳ воспроизвести ошибку |
| `INFO_BALLOON`                   | внутри `BALLOON_CONTAINER`, подтип info                  | ⏳     |

---

## Toolbox (заменяет Browse-панель когда открыт view)

| Имя                              | Селектор                                                | Статус |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `TOOLBOX_TAB_PANE`               | `[name="tab-pane-Toolbox"]`                             | ✅     |
| `TOOLBOX_PANE_SEARCH`            | `[name="pane-Search"]`                                  | ✅     |
| `TOOLBOX_PANE_VIEWERS`           | `[name="pane-Viewers"]`                                 | ✅     |
| `TOOLBOX_PANE_LAYOUTS`           | `[name="pane-Layouts"]`                                 | ✅     |
| `TOOLBOX_PIN_TOOLBOX`            | `[aria-label="Pin toolbox"]`                            | ✅     |

---

## Сервисные

| Имя                              | Селектор                                                | Статус |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `APP_LOADED` (= `RIBBON`)        | `.d4-ribbon` — признак загруженного приложения          | ✅     |
| `PROGRESS_INDICATOR`             | индикатор загрузки/спиннер                              | ⏳     |
