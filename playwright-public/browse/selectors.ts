import { Page, Locator } from '@playwright/test';

// =============================================================================
// Browse selectors — единый источник истины для тестов из e2e/browse/.
// См. selectors.md (человекочитаемое оглавление).
//
// Значения сняты на https://dev.datagrok.ai 2026-06-03.
// Маркеры:
//   * без комментария   — подтверждено в DOM
//   * // TODO ...       — требует уточнения (не воспроизведено в разведке)
// =============================================================================

// ---------- Sidebar -----------------------------------------------------------

export const SIDEBAR = '.d4-tab-header-stripe.layout-sidebar.vertical';
export const SIDEBAR_BROWSE_ICON = `${SIDEBAR} .d4-tab-header[name="Browse"]`;
export const SIDEBAR_BROWSE_ICON_BY_DATA = `${SIDEBAR} [data-source="tab-pane-Browse"]`;
export const SIDEBAR_FAVORITES_ICON = `${SIDEBAR} .d4-tab-header[name="Favorites"]`;
export const SIDEBAR_FAVORITES_ICON_BY_DATA = `${SIDEBAR} [data-source="tab-pane-Favorites"]`;
export const SIDEBAR_TAB_SELECTED = '.d4-tab-header.selected';
export const SIDEBAR_VIEW_BADGE = '.d4-tab-header .grok-icon-badge'; // TODO: подтвердить класс бейджа

// ---------- Browse panel — общая структура ------------------------------------

export const BROWSE_PANEL_TAB_PANE = '[name="tab-pane-Browse"]';
export const BROWSE_HEADER = '.panel-titlebar.disable-selection.grok-browse-header';
export const BROWSE_HEADER_ICONS_GROUP = `${BROWSE_HEADER} .grok-browse-icons`;
export const BROWSE_HEADER_ICONS = `${BROWSE_HEADER} i.grok-icon`;

// Конкретные иконки в шапке Browse — по name атрибуту (стабильнее всего).
export const BROWSE_HEADER_HOME = `${BROWSE_HEADER} [name="icon-home"]`;
export const BROWSE_HEADER_IMPORT_FILE = `${BROWSE_HEADER} [name="icon-folder-open"]`;
export const BROWSE_HEADER_IMPORT_TEXT = `${BROWSE_HEADER} [name="icon-file-alt"]`;
export const BROWSE_HEADER_REFRESH = `${BROWSE_HEADER} [name="icon-sync"]`;
export const BROWSE_HEADER_COLLAPSE_ALL = `${BROWSE_HEADER} [name="icon-chevron-double-up"]`;
export const BROWSE_HEADER_LOCATE = `${BROWSE_HEADER} [name="icon-dot-circle"]`;
export const BROWSE_HEADER_CLOSE_PANEL = `${BROWSE_HEADER} .grok-font-icon-close`;

// ---------- Browse tree — узлы и стрелки --------------------------------------

export const TREE_NODE_GROUP_LABEL = '.d4-tree-view-group-label';
export const TREE_NODE_ITEM_LABEL = '.d4-tree-view-item-label';
export const TREE_NODE_LABEL_ANY = `${TREE_NODE_GROUP_LABEL}, ${TREE_NODE_ITEM_LABEL}`;
export const TREE_NODE_CONTAINER = '.d4-tree-view-node';
export const TREE_EXPAND_ARROW = '.d4-tree-view-tri';
export const TREE_EXPAND_ARROW_EXPANDED = '.d4-tree-view-tri.d4-tree-view-tri-expanded';
export const TREE_NODE_DROP = '.d4-tree-view-node.d4-tree-drop';

// Динамические локаторы — по имени узла.
// Точный селектор по name (например [name="tree-node-My-files"]) работает только когда меню/aria
// уже навешены; для общего случая используем поиск по тексту label.

export const treeNodeByName = (page: Page, name: string): Locator =>
  page.locator(TREE_NODE_LABEL_ANY, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

export const treeGroupByName = (page: Page, name: string): Locator =>
  page.locator(TREE_NODE_GROUP_LABEL, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

export const treeItemByName = (page: Page, name: string): Locator =>
  page.locator(TREE_NODE_ITEM_LABEL, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

// Контейнер узла (для events / атрибутов).
export const treeNodeContainer = (page: Page, name: string): Locator =>
  treeNodeByName(page, name)
    .locator('xpath=ancestor::*[contains(@class,"d4-tree-view-node")][1]');

/**
 * Stable locator by tree path (uses the `name="tree-<Root>---<Child>---..."` attribute).
 * Example: treeNodeByPath(page, ['Apps', 'Tutorials']) -> [name="tree-Apps---Tutorials"]
 */
export const treeNodeByPath = (page: Page, segments: string[]): Locator =>
  page.locator(`[name="tree-${segments.map((s) => s.replace(/ /g, '-')).join('---')}"]`);

// Стрелка expand/collapse слева от узла.
export const treeExpandArrow = (page: Page, name: string): Locator =>
  treeNodeContainer(page, name).locator(TREE_EXPAND_ARROW).first();

// Контейнер с детьми развёрнутого узла.
export const treeNodeChildren = (page: Page, parentName: string): Locator =>
  treeNodeContainer(page, parentName).locator(TREE_NODE_CONTAINER);

// ---------- Context menu (right-click) ----------------------------------------

export const CONTEXT_MENU = '.d4-menu-item-container.d4-vert-menu.d4-menu-popup';
export const CONTEXT_MENU_ITEM = '.d4-menu-item.d4-menu-item-vert';
export const CONTEXT_MENU_ITEM_LABEL = '.d4-menu-item-label';

// По атрибуту d4-name="<точное имя пункта>" (например d4-name="Add to favorites").
export const contextMenuItemByName = (page: Page, d4Name: string): Locator =>
  page.locator(`${CONTEXT_MENU_ITEM}[d4-name="${d4Name}"]`).first();

// По тексту (универсально).
export const contextMenuItem = (page: Page, text: string): Locator =>
  page.locator(CONTEXT_MENU_ITEM_LABEL, { hasText: new RegExp(`^${escapeRegExp(text)}$`) }).first();

// Точные тексты пунктов меню на dev (важно: регистр и многоточия).
export const CONTEXT_MENU_BROWSE = 'Browse';
export const CONTEXT_MENU_ADD_FAVORITES = 'Add to favorites';
export const CONTEXT_MENU_SHARE = 'Share...';
export const CONTEXT_MENU_RENAME = 'Rename...';
export const CONTEXT_MENU_DELETE = 'Delete...';
export const CONTEXT_MENU_CLONE = 'Clone...';
export const CONTEXT_MENU_EDIT = 'Edit...';
export const CONTEXT_MENU_DOWNLOAD = 'Download'; // TODO: подтвердить точный текст для файлов
export const CONTEXT_MENU_BROWSE_SCHEMA = 'Browse schema'; // TODO: подтвердить для connection
export const CONTEXT_MENU_RUN = 'Run'; // TODO: подтвердить для модели/функции
export const CONTEXT_MENU_TEST_CONNECTION = 'Test connection';
export const CONTEXT_MENU_CLEAR_CACHE = 'Clear cache';

// ---------- Context Panel (правая панель свойств) -----------------------------

export const CONTEXT_PANEL = '.grok-prop-panel';
export const CONTEXT_PANEL_BASE = '.panel-base.splitter-container-vertical:has(.grok-prop-panel)';
export const CONTEXT_PANEL_HEADER = `${CONTEXT_PANEL_BASE} > .panel-titlebar.disable-selection`;

// Кнопки в шапке Context Panel — через aria-label (надёжно).
export const CONTEXT_PANEL_BACK = `${CONTEXT_PANEL_HEADER} [aria-label="Back"]`;
export const CONTEXT_PANEL_FORWARD = `${CONTEXT_PANEL_HEADER} [aria-label="Forward"]`;
export const CONTEXT_PANEL_CLONE = `${CONTEXT_PANEL_HEADER} [aria-label="Clone and detach"]`;
export const CONTEXT_PANEL_COLLAPSE_ALL = `${CONTEXT_PANEL_HEADER} [aria-label="Collapse all"]`;
export const CONTEXT_PANEL_EXPAND_ALL = `${CONTEXT_PANEL_HEADER} [aria-label="Expand all"]`;
export const CONTEXT_PANEL_HELP = `${CONTEXT_PANEL_HEADER} [aria-label="Show Help"]`;
export const CONTEXT_PANEL_STAR = `${CONTEXT_PANEL_HEADER} [aria-label="Favorites"]`;

// Содержимое.
export const CONTEXT_PANEL_INNER = `${CONTEXT_PANEL} .grok-entity-prop-panel`;
export const CONTEXT_PANEL_ACCORDION = `${CONTEXT_PANEL} .d4-accordion`;
export const CONTEXT_PANEL_ACCORDION_TITLE = '.d4-accordion-title';
export const CONTEXT_PANEL_ACCORDION_PANE = '.d4-accordion-pane';
export const infoPaneByName = (page: Page, name: string): Locator =>
  page.locator(`${CONTEXT_PANEL} ${CONTEXT_PANEL_ACCORDION_PANE}[name="pane-${name}"]`).first();

// ---------- View / tabs / top ribbon ------------------------------------------

export const RIBBON = '.d4-ribbon';
export const VIEW_TAB = '.tab-handle.disable-selection';
export const VIEW_TAB_SELECTED = '.tab-handle.tab-handle-selected';
export const viewTabHandle = (page: Page, viewName: string): Locator =>
  page.locator(`[name="view-handle: ${viewName}"]`);
export const VIEW_TAB_TEXT = '.tab-handle-text';
export const VIEW_TAB_CLOSE = '.tab-handle-close-button[name="Close"]';

// Top menu (горизонтальное основное меню над таблицей).
export const TOP_MENU = '.grok-main-menu';
export const TOP_MENU_ROOT = `${TOP_MENU}[d4-name="top"]`;
export const topMenuSection = (page: Page, name: string): Locator =>
  page.locator(`${TOP_MENU} [d4-name="${name}"]`).first();
export const topMenuItem = (page: Page, path: string[]): Locator =>
  page.locator(`[name="div-${path.join('---').replace(/ /g, '-')}"]`).first();

// View close (в шапке открытого view).
export const VIEW_CLOSE = '[aria-label="Close view"]';
export const VIEW_CLOSE_BY_NAME = (name: string) => `[aria-label="Close ${name}"]`;

// View tab pin (browsing preview): appears as a thumbtack icon inside `.tab-handle-button`
// when the view is unpinned. The aria-label is the user-facing hint text.
export const TOOLBOX_PIN = '[aria-label="Pin toolbox"]';
export const VIEW_TAB_PREVIEW_TEXT = '.tab-handle-text.grok-browse-preview';
export const VIEW_TAB_PREVIEW_TAB = `${'.tab-handle.disable-selection'}:has(${VIEW_TAB_PREVIEW_TEXT})`;
export const VIEW_TAB_PIN_BUTTON = '.tab-handle-button i.fa-thumbtack';
export const VIEW_TAB_PIN_BY_ARIA = '[aria-label="This is Browse preview. Click to keep it open"]';

// View selector (combo-popup) — shows a list of all open views and a count badge.
export const VIEW_SELECTOR = '[name="view selector"]';
export const VIEW_SELECTOR_COUNT_BADGE = `${VIEW_SELECTOR} span`;

// Hamburger / Split right / Split down — в Top Menu не найдены при simpleMode=true.
// Скорее всего доступны только в windows mode (simpleMode=false) через dock manager.
export const VIEW_HAMBURGER = ''; // TODO: воспроизвести в windows mode
export const VIEW_HAMBURGER_SPLIT_RIGHT = 'Split right'; // TODO
export const VIEW_HAMBURGER_SPLIT_DOWN = 'Split down'; // TODO

// ---------- Status bar --------------------------------------------------------

export const STATUS_BAR = '.layout-status-bar';
export const STATUS_BAR_TASK = `${STATUS_BAR} .d4-task-bar`;
export const STATUS_BAR_VIEW_PANEL = `${STATUS_BAR} .d4-view-status-panel`;
export const STATUS_BAR_GLOBAL_PANEL = `${STATUS_BAR} .d4-global-status-panel`;

// Status bar toggles (icons have NO name/aria-label — identified by FontAwesome class).
// Hover tooltips confirmed on dev:
//   fa-window-maximize → "Tabs"           (Tabs ↔ Windows mode)
//   fa-presentation    → "Presentation mode F7"
//   fa-ballot          → "Toolbox"
//   fa-sliders-h       → "Context Panel F4"
//   fa-info            → "Context Help F1"
//   fa-value-absolute  → "Variables ALT+V"
//   fa-terminal        → "Console"
//   fa-user-robot      → "AI"
export const STATUS_BAR_MODE_TABS = `${STATUS_BAR} i.fa-window-maximize`;
export const STATUS_BAR_MODE_PRESENTATION = `${STATUS_BAR} i.fa-presentation`;
export const STATUS_BAR_TOGGLE_TOOLBOX = `${STATUS_BAR} i.fa-ballot`;
export const STATUS_BAR_TOGGLE_CONTEXT_PANEL = `${STATUS_BAR} i.fa-sliders-h`;
export const STATUS_BAR_TOGGLE_HELP = `${STATUS_BAR} i.fa-info`;
export const STATUS_BAR_TOGGLE_VARIABLES = `${STATUS_BAR} i.fa-value-absolute`;
export const STATUS_BAR_TERMINAL = `${STATUS_BAR} i.fa-terminal`;
export const STATUS_BAR_AI_HELPER = `${STATUS_BAR} i.fa-user-robot`;

// ---------- Home Page / global search -----------------------------------------

export const HOME_VIEW_HANDLE = '[name="view-handle: Home"]';
export const HOME_GLOBAL_SEARCH_INPUT = 'input.power-search-search-everywhere-input';
export const HOME_WIDGETS = '.power-search-input-container, .power-search-widgets'; // TODO: уточнить
export const HOME_SEARCH_RESULTS = '.power-search-results'; // TODO

// ---------- Filtering / lists -------------------------------------------------

// Filter panel appears in gallery views (Users / Files / Apps / Models). Confirmed on dev.
export const FILTER_PANEL = '.d4-filter-panel';
export const FILTER_PANEL_SECTION = `${FILTER_PANEL} .d4-filter-panel-section`;
export const FILTER_PANEL_HEADER = `${FILTER_PANEL_SECTION} .d4-filter-panel-header`;
export const FILTER_QUICK_TAGS = `${FILTER_PANEL_SECTION} .d4-tag-group .d4-tag`;
export const FILTER_TOGGLE = '[aria-label="Toggle filters"]';
export const filterQuickTag = (page: Page, text: string): Locator =>
  page.locator(FILTER_QUICK_TAGS, { hasText: new RegExp(`^${escapeRegExp(text)}$`) }).first();
export const filterProperty = (page: Page, name: string): Locator =>
  page.locator(`${FILTER_PANEL} .ui-label, ${FILTER_PANEL} label`,
    { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();
export const LIST_SEARCH_INPUT = '.d4-search-input.ui-input-editor';
export const LIST_EMPTY_STATE = '.grok-empty-state'; // TODO confirm when there's truly no result

// ---------- Notifications / errors / balloons --------------------------------

export const BALLOON_CONTAINER = '.d4-balloon-container';
export const ERROR_BALLOON = `${BALLOON_CONTAINER} .d4-balloon-error, ${BALLOON_CONTAINER} [class*="error"]`; // TODO
export const INFO_BALLOON = `${BALLOON_CONTAINER} .d4-balloon-info, ${BALLOON_CONTAINER} [class*="info"]`; // TODO

// ---------- Toolbox -----------------------------------------------------------

export const TOOLBOX_TAB_PANE = '[name="tab-pane-Toolbox"]';
export const TOOLBOX_PANE_SEARCH = '[name="pane-Search"]';
export const TOOLBOX_PANE_VIEWERS = '[name="pane-Viewers"]';
export const TOOLBOX_PANE_LAYOUTS = '[name="pane-Layouts"]';
export const TOOLBOX_PIN_TOOLBOX = '[aria-label="Pin toolbox"]';

// ---------- Сервисные ---------------------------------------------------------

export const APP_LOADED = RIBBON; // признак, что приложение загрузилось
export const PROGRESS_INDICATOR = '.grok-progress'; // TODO

// ---------- Утилиты -----------------------------------------------------------

function escapeRegExp(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
