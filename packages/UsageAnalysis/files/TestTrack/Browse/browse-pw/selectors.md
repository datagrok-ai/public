# Browse selectors — reference

Catalog of logical selector names used in `playwright-tests/e2e/browse/*.test.ts`.

**Source of truth:** `selectors.ts`. This file is a short table of contents for reviewing the manual test docs and Test Track.

**Legend:**
- ✅ value captured on dev and confirmed
- ⏳ TODO — needs further investigation (see the `TODO` marker in `selectors.ts`)

Investigation date: **2026-06-03** (dev.datagrok.ai, user `opavlenko+playwright@datagrok.ai`).

---

## Sidebar (leftmost strip)

| Name                             | Selector                                       | Status |
| -------------------------------- | ---------------------------------------------- | ------ |
| `SIDEBAR`                        | `.d4-tab-header-stripe.layout-sidebar.vertical` | ✅     |
| `SIDEBAR_BROWSE_ICON`            | `.d4-tab-header[name="Browse"]`                | ✅     |
| `SIDEBAR_BROWSE_ICON_BY_DATA`    | `[data-source="tab-pane-Browse"]`              | ✅     |
| `SIDEBAR_FAVORITES_ICON`         | `.d4-tab-header[name="Favorites"]`             | ✅     |
| `SIDEBAR_TAB_SELECTED`           | `.d4-tab-header.selected`                      | ✅     |
| `SIDEBAR_VIEW_BADGE`             | counter badge for views of the same type       | ⏳     |

---

## Browse panel — header and icons

| Name                             | Selector                                                        | Status |
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
| `BROWSE_HEADER_CLOSE_PANEL`      | `.grok-font-icon-close` (inside `BROWSE_HEADER`)                | ✅     |

> In total, the Browse header has **6 functional icons + 1 close** (= 7 clickable, not counting the `.grok-browse-icons` container).

---

## Browse tree

| Name                             | Selector                                                  | Status |
| -------------------------------- | --------------------------------------------------------- | ------ |
| `TREE_NODE_GROUP_LABEL`          | `.d4-tree-view-group-label`                               | ✅     |
| `TREE_NODE_ITEM_LABEL`           | `.d4-tree-view-item-label`                                | ✅     |
| `TREE_NODE_LABEL_ANY`            | both classes above                                        | ✅     |
| `TREE_NODE_CONTAINER`            | `.d4-tree-view-node`                                       | ✅     |
| `TREE_EXPAND_ARROW`              | `.d4-tree-view-tri`                                        | ✅     |
| `TREE_EXPAND_ARROW_EXPANDED`     | `.d4-tree-view-tri.d4-tree-view-tri-expanded`              | ✅     |
| `TREE_NODE_DROP`                 | `.d4-tree-view-node.d4-tree-drop` (drop zone for DnD)      | ✅     |

**Dynamic locators:**

| Function                   | Description                                    | Status |
| -------------------------- | ---------------------------------------------- | ------ |
| `treeNodeByName(p, name)`  | Tree node by exact name                        | ✅     |
| `treeGroupByName(p, name)` | Group node (with children)                     | ✅     |
| `treeItemByName(p, name)`  | Leaf node (without children)                   | ✅     |
| `treeNodeContainer(p, n)`  | `.d4-tree-view-node` container for a node      | ✅     |
| `treeExpandArrow(p, name)` | Expand/collapse arrow on the left              | ✅     |
| `treeNodeChildren(p, n)`   | Child nodes of an expanded group               | ✅     |

> When the menu/aria are already initialized, there is also a direct `[name="tree-node-<Name>"]` (for example `tree-node-My-files`), but it only fires reliably after the first interaction with the node.

---

## Context menu (right-click)

| Name                             | Selector                                                   | Status |
| -------------------------------- | ---------------------------------------------------------- | ------ |
| `CONTEXT_MENU`                   | `.d4-menu-item-container.d4-vert-menu.d4-menu-popup`       | ✅     |
| `CONTEXT_MENU_ITEM`              | `.d4-menu-item.d4-menu-item-vert`                          | ✅     |
| `CONTEXT_MENU_ITEM_LABEL`        | `.d4-menu-item-label`                                      | ✅     |
| `contextMenuItemByName(p, name)` | by the `d4-name="..."` attribute                           | ✅     |
| `contextMenuItem(p, text)`       | by the item text                                           | ✅     |

**Exact item texts (important: case and ellipses):**

| Name                             | Item text            |
| -------------------------------- | -------------------- |
| `CONTEXT_MENU_BROWSE`            | `Browse`             |
| `CONTEXT_MENU_ADD_FAVORITES`     | **`Add to favorites`** (lowercase `f`!) |
| `CONTEXT_MENU_SHARE`             | `Share...`           |
| `CONTEXT_MENU_RENAME`            | `Rename...`          |
| `CONTEXT_MENU_DELETE`            | `Delete...`          |
| `CONTEXT_MENU_CLONE`             | `Clone...`           |
| `CONTEXT_MENU_EDIT`              | `Edit...`            |
| `CONTEXT_MENU_TEST_CONNECTION`   | `Test connection`    |
| `CONTEXT_MENU_CLEAR_CACHE`       | `Clear cache`        |
| `CONTEXT_MENU_DOWNLOAD`          | ⏳ to be clarified    |
| `CONTEXT_MENU_BROWSE_SCHEMA`     | ⏳ to be clarified    |
| `CONTEXT_MENU_RUN`               | ⏳ to be clarified    |

---

## Context Panel (right-hand properties panel)

| Name                             | Selector                                                                     | Status |
| -------------------------------- | ---------------------------------------------------------------------------- | ------ |
| `CONTEXT_PANEL`                  | `.grok-prop-panel`                                                            | ✅     |
| `CONTEXT_PANEL_BASE`             | `.panel-base.splitter-container-vertical`                                     | ✅     |
| `CONTEXT_PANEL_HEADER`           | `.panel-titlebar.disable-selection` (first one inside `CONTEXT_PANEL_BASE`)   | ✅     |
| `CONTEXT_PANEL_BACK`             | `[aria-label="Back"]` (`fa-chevron-left`, `name="icon-chevron-left"`)         | ✅     |
| `CONTEXT_PANEL_FORWARD`          | `[aria-label="Forward"]` (`fa-chevron-right`)                                 | ✅     |
| `CONTEXT_PANEL_CLONE`            | `[aria-label="Clone and detach"]` (`fa-clone`) — case CtxPanel-05 removed      | ✅     |
| `CONTEXT_PANEL_COLLAPSE_ALL`     | `[aria-label="Collapse all"]` (`fa-chevron-square-up`)                        | ✅     |
| `CONTEXT_PANEL_EXPAND_ALL`       | `[aria-label="Expand all"]` (`fa-chevron-square-down`)                        | ✅     |
| `CONTEXT_PANEL_HELP`             | `[aria-label="Show Help"]` (`grok-font-icon-help`)                            | ✅     |
| `CONTEXT_PANEL_STAR`             | `[aria-label="Favorites"]` (`fa-star`, `name="icon-star"`)                    | ✅     |
| `CONTEXT_PANEL_INNER`            | `.grok-prop-panel .grok-entity-prop-panel`                                    | ✅     |
| `CONTEXT_PANEL_ACCORDION`        | `.grok-prop-panel .d4-accordion`                                              | ✅     |
| `CONTEXT_PANEL_ACCORDION_PANE`   | `.d4-accordion-pane`                                                          | ✅     |
| `infoPaneByName(p, name)`        | `.d4-accordion-pane[name="pane-<Name>"]` (for example `pane-Search`)          | ✅     |

> The star (`CONTEXT_PANEL_STAR`) is in the Context Panel header. The "active/orange" state will be captured during the first Browse-Fav-02 check (probably via an additional class or color).

---

## View / tabs / top ribbon

| Name                             | Selector                                                  | Status |
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
| `VIEW_TAB_PIN`                   | `[aria-label*="Pin"]` for view tabs                        | ⏳     |
| `VIEW_HAMBURGER` / `Split right/down` | Not found in `simpleMode=true`; probably in windows mode | ⏳     |

---

## Status bar (bottom of the screen; visible when a view is open)

| Name                             | Selector                                                | Status |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `STATUS_BAR`                     | `.layout-status-bar`                                    | ✅     |
| `STATUS_BAR_TASK`                | `.layout-status-bar .d4-task-bar`                       | ✅     |
| `STATUS_BAR_VIEW_PANEL`          | `.layout-status-bar .d4-view-status-panel` (Rows/Cols/Selected) | ✅ |
| `STATUS_BAR_GLOBAL_PANEL`        | `.layout-status-bar .d4-global-status-panel`            | ✅     |
| `STATUS_BAR_AI_HELPER`           | `.layout-status-bar i.fa-user-robot`                    | ✅     |
| `STATUS_BAR_TERMINAL`            | `.layout-status-bar i.fa-terminal`                      | ✅     |
| `STATUS_BAR_MODE_MAXIMIZE`       | `.layout-status-bar i.fa-window-maximize`               | ⏳ confirm mapping → Tabbed/Simple/Presentation |
| `STATUS_BAR_MODE_PRESENTATION`   | `.layout-status-bar i.fa-presentation`                  | ⏳     |

> The status bar icons have no `name` / `title` / `aria-label`. I will clarify the exact icon ↔ Tabbed/Simple/Presentation mapping during the Browse-Modes-01 test.

---

## Home Page / global search

| Name                             | Selector                                                | Status |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `HOME_VIEW_HANDLE`               | `[name="view-handle: Home"]`                            | ✅     |
| `HOME_GLOBAL_SEARCH_INPUT`       | `input.power-search-search-everywhere-input` (placeholder: `Search everywhere. Try "aspirin" or "7JZK"`) | ✅ |
| `HOME_WIDGETS`                   | widgets container                                       | ⏳     |
| `HOME_SEARCH_RESULTS`            | search results area                                     | ⏳     |

---

## Filtering / lists

| Name                             | Selector                                                | Status |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `FILTER_PANEL`                   | filter panel in the gallery (Apps/Files etc.)           | ⏳ reproduce in the gallery |
| `FILTER_TOGGLE`                  | `[aria-label="Toggle filters"]`                         | ✅     |
| `filterProperty(p, name)`        | filter property by name                                 | ⏳     |
| `LIST_SEARCH_INPUT`              | `.d4-search-input.ui-input-editor`                      | ✅     |
| `LIST_EMPTY_STATE`               | empty state                                             | ⏳     |

---

## Notifications / errors / balloons

| Name                             | Selector                                                | Status |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `BALLOON_CONTAINER`              | `.d4-balloon-container`                                 | ✅     |
| `ERROR_BALLOON`                  | inside `BALLOON_CONTAINER`, error subtype               | ⏳ reproduce an error |
| `INFO_BALLOON`                   | inside `BALLOON_CONTAINER`, info subtype                | ⏳     |

---

## Toolbox (replaces the Browse panel when a view is open)

| Name                             | Selector                                                | Status |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `TOOLBOX_TAB_PANE`               | `[name="tab-pane-Toolbox"]`                             | ✅     |
| `TOOLBOX_PANE_SEARCH`            | `[name="pane-Search"]`                                  | ✅     |
| `TOOLBOX_PANE_VIEWERS`           | `[name="pane-Viewers"]`                                 | ✅     |
| `TOOLBOX_PANE_LAYOUTS`           | `[name="pane-Layouts"]`                                 | ✅     |
| `TOOLBOX_PIN_TOOLBOX`            | `[aria-label="Pin toolbox"]`                            | ✅     |

---

## Service

| Name                             | Selector                                                | Status |
| -------------------------------- | ------------------------------------------------------- | ------ |
| `APP_LOADED` (= `RIBBON`)        | `.d4-ribbon` — indicates the application has loaded      | ✅     |
| `PROGRESS_INDICATOR`             | loading indicator/spinner                               | ⏳     |
