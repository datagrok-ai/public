import { Page, Locator } from '@playwright/test';

export const SIDEBAR = '.d4-tab-header-stripe.layout-sidebar.vertical';
export const SIDEBAR_BROWSE_ICON = `${SIDEBAR} .d4-tab-header[name="Browse"]`;
export const SIDEBAR_BROWSE_ICON_BY_DATA = `${SIDEBAR} [data-source="tab-pane-Browse"]`;
export const SIDEBAR_FAVORITES_ICON = `${SIDEBAR} .d4-tab-header[name="Favorites"]`;
export const SIDEBAR_FAVORITES_ICON_BY_DATA = `${SIDEBAR} [data-source="tab-pane-Favorites"]`;
export const SIDEBAR_TAB_SELECTED = '.d4-tab-header.selected';
export const SIDEBAR_VIEW_BADGE = '.d4-tab-header .grok-icon-badge'; 

export const BROWSE_PANEL_TAB_PANE = '[name="tab-pane-Browse"]';
export const BROWSE_HEADER = '.panel-titlebar.disable-selection.grok-browse-header';
export const BROWSE_HEADER_ICONS_GROUP = `${BROWSE_HEADER} .grok-browse-icons`;
export const BROWSE_HEADER_ICONS = `${BROWSE_HEADER} i.grok-icon`;

export const BROWSE_HEADER_HOME = `${BROWSE_HEADER} [name="icon-home"]`;
export const BROWSE_HEADER_IMPORT_FILE = `${BROWSE_HEADER} [name="icon-folder-open"]`;
export const BROWSE_HEADER_IMPORT_TEXT = `${BROWSE_HEADER} [name="icon-file-alt"]`;
export const BROWSE_HEADER_REFRESH = `${BROWSE_HEADER} [name="icon-sync"]`;
export const BROWSE_HEADER_COLLAPSE_ALL = `${BROWSE_HEADER} [name="icon-chevron-double-up"]`;
export const BROWSE_HEADER_LOCATE = `${BROWSE_HEADER} [name="icon-dot-circle"]`;
export const BROWSE_HEADER_CLOSE_PANEL = `${BROWSE_HEADER} .grok-font-icon-close`;

export const TREE_NODE_GROUP_LABEL = '.d4-tree-view-group-label';
export const TREE_NODE_ITEM_LABEL = '.d4-tree-view-item-label';
export const TREE_NODE_LABEL_ANY = `${TREE_NODE_GROUP_LABEL}, ${TREE_NODE_ITEM_LABEL}`;
export const TREE_NODE_CONTAINER = '.d4-tree-view-node';
export const TREE_EXPAND_ARROW = '.d4-tree-view-tri';
export const TREE_EXPAND_ARROW_EXPANDED = '.d4-tree-view-tri.d4-tree-view-tri-expanded';
export const TREE_NODE_DROP = '.d4-tree-view-node.d4-tree-drop';

export const treeNodeByName = (page: Page, name: string): Locator =>
  page.locator(TREE_NODE_LABEL_ANY, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

export const treeGroupByName = (page: Page, name: string): Locator =>
  page.locator(TREE_NODE_GROUP_LABEL, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

export const treeItemByName = (page: Page, name: string): Locator =>
  page.locator(TREE_NODE_ITEM_LABEL, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

export const treeNodeContainer = (page: Page, name: string): Locator =>
  treeNodeByName(page, name)
    .locator('xpath=ancestor::*[contains(@class,"d4-tree-view-node")][1]');

export const treeNodeByPath = (page: Page, segments: string[]): Locator =>
  page.locator(`[name="tree-${segments.map((s) => s.replace(/ /g, '-')).join('---')}"]`);

export const treeExpandArrow = (page: Page, name: string): Locator =>
  treeNodeContainer(page, name).locator(TREE_EXPAND_ARROW).first();

export const treeNodeChildren = (page: Page, parentName: string): Locator =>
  treeNodeContainer(page, parentName).locator(TREE_NODE_CONTAINER);

export const CONTEXT_MENU = '.d4-menu-item-container.d4-vert-menu.d4-menu-popup';
export const CONTEXT_MENU_ITEM = '.d4-menu-item.d4-menu-item-vert';
export const CONTEXT_MENU_ITEM_LABEL = '.d4-menu-item-label';

export const contextMenuItemByName = (page: Page, d4Name: string): Locator =>
  page.locator(`${CONTEXT_MENU_ITEM}[d4-name="${d4Name}"]`).first();

export const contextMenuItem = (page: Page, text: string): Locator =>
  page.locator(CONTEXT_MENU_ITEM_LABEL, { hasText: new RegExp(`^${escapeRegExp(text)}$`) }).first();

export const CONTEXT_MENU_BROWSE = 'Browse';
export const CONTEXT_MENU_ADD_FAVORITES = 'Add to favorites';
export const CONTEXT_MENU_SHARE = 'Share...';
export const CONTEXT_MENU_RENAME = 'Rename...';
export const CONTEXT_MENU_DELETE = 'Delete...';
export const CONTEXT_MENU_CLONE = 'Clone...';
export const CONTEXT_MENU_EDIT = 'Edit...';
export const CONTEXT_MENU_DOWNLOAD = 'Download'; 
export const CONTEXT_MENU_BROWSE_SCHEMA = 'Browse schema'; 
export const CONTEXT_MENU_RUN = 'Run'; 
export const CONTEXT_MENU_TEST_CONNECTION = 'Test connection';
export const CONTEXT_MENU_CLEAR_CACHE = 'Clear cache';

export const CONTEXT_PANEL = '.grok-prop-panel';
export const CONTEXT_PANEL_BASE = '.panel-base.splitter-container-vertical:has(.grok-prop-panel)';
export const CONTEXT_PANEL_HEADER = `${CONTEXT_PANEL_BASE} > .panel-titlebar.disable-selection`;

export const CONTEXT_PANEL_BACK = `${CONTEXT_PANEL_HEADER} [aria-label="Back"]`;
export const CONTEXT_PANEL_FORWARD = `${CONTEXT_PANEL_HEADER} [aria-label="Forward"]`;
export const CONTEXT_PANEL_CLONE = `${CONTEXT_PANEL_HEADER} [aria-label="Clone and detach"]`;
export const CONTEXT_PANEL_COLLAPSE_ALL = `${CONTEXT_PANEL_HEADER} [aria-label="Collapse all"]`;
export const CONTEXT_PANEL_EXPAND_ALL = `${CONTEXT_PANEL_HEADER} [aria-label="Expand all"]`;
export const CONTEXT_PANEL_HELP = `${CONTEXT_PANEL_HEADER} [aria-label="Show Help"]`;
export const CONTEXT_PANEL_STAR = `${CONTEXT_PANEL_HEADER} [aria-label="Favorites"]`;

export const CONTEXT_PANEL_INNER = `${CONTEXT_PANEL} .grok-entity-prop-panel`;
export const CONTEXT_PANEL_ACCORDION = `${CONTEXT_PANEL} .d4-accordion`;
export const CONTEXT_PANEL_ACCORDION_TITLE = '.d4-accordion-title';
export const CONTEXT_PANEL_ACCORDION_PANE = '.d4-accordion-pane';
export const infoPaneByName = (page: Page, name: string): Locator =>
  page.locator(`${CONTEXT_PANEL} ${CONTEXT_PANEL_ACCORDION_PANE}[name="pane-${name}"]`).first();

export const RIBBON = '.d4-ribbon';
export const VIEW_TAB = '.tab-handle.disable-selection';
export const VIEW_TAB_SELECTED = '.tab-handle.tab-handle-selected';
export const viewTabHandle = (page: Page, viewName: string): Locator =>
  page.locator(`[name="view-handle: ${viewName}"]`);
export const VIEW_TAB_TEXT = '.tab-handle-text';
export const VIEW_TAB_CLOSE = '.tab-handle-close-button[name="Close"]';

export const TOP_MENU = '.grok-main-menu';
export const TOP_MENU_ROOT = `${TOP_MENU}[d4-name="top"]`;
export const topMenuSection = (page: Page, name: string): Locator =>
  page.locator(`${TOP_MENU} [d4-name="${name}"]`).first();
export const topMenuItem = (page: Page, path: string[]): Locator =>
  page.locator(`[name="div-${path.join('---').replace(/ /g, '-')}"]`).first();

export const VIEW_CLOSE = '[aria-label="Close view"]';
export const VIEW_CLOSE_BY_NAME = (name: string) => `[aria-label="Close ${name}"]`;

export const TOOLBOX_PIN = '[aria-label="Pin toolbox"]';
export const VIEW_TAB_PREVIEW_TEXT = '.tab-handle-text.grok-browse-preview';
export const VIEW_TAB_PREVIEW_TAB = `${'.tab-handle.disable-selection'}:has(${VIEW_TAB_PREVIEW_TEXT})`;
export const VIEW_TAB_PIN_BUTTON = '.tab-handle-button i.fa-thumbtack';
export const VIEW_TAB_PIN_BY_ARIA = '[aria-label="This is Browse preview. Click to keep it open"]';

export const VIEW_SELECTOR = '[name="view selector"]';
export const VIEW_SELECTOR_COUNT_BADGE = `${VIEW_SELECTOR} span`;

export const VIEW_HAMBURGER = ''; 
export const VIEW_HAMBURGER_SPLIT_RIGHT = 'Split right'; 
export const VIEW_HAMBURGER_SPLIT_DOWN = 'Split down'; 

export const STATUS_BAR = '.layout-status-bar';
export const STATUS_BAR_TASK = `${STATUS_BAR} .d4-task-bar`;
export const STATUS_BAR_VIEW_PANEL = `${STATUS_BAR} .d4-view-status-panel`;
export const STATUS_BAR_GLOBAL_PANEL = `${STATUS_BAR} .d4-global-status-panel`;

export const STATUS_BAR_MODE_TABS = `${STATUS_BAR} i.fa-window-maximize`;
export const STATUS_BAR_MODE_PRESENTATION = `${STATUS_BAR} i.fa-presentation`;
export const STATUS_BAR_TOGGLE_TOOLBOX = `${STATUS_BAR} i.fa-ballot`;
export const STATUS_BAR_TOGGLE_CONTEXT_PANEL = `${STATUS_BAR} i.fa-sliders-h`;
export const STATUS_BAR_TOGGLE_HELP = `${STATUS_BAR} i.fa-info`;
export const STATUS_BAR_TOGGLE_VARIABLES = `${STATUS_BAR} i.fa-value-absolute`;
export const STATUS_BAR_TERMINAL = `${STATUS_BAR} i.fa-terminal`;
export const STATUS_BAR_AI_HELPER = `${STATUS_BAR} i.fa-user-robot`;

export const HOME_VIEW_HANDLE = '[name="view-handle: Home"]';
export const HOME_GLOBAL_SEARCH_INPUT = 'input.power-search-search-everywhere-input';
export const HOME_WIDGETS = '.power-search-input-container, .power-search-widgets'; 
export const HOME_SEARCH_RESULTS = '.power-search-results'; 

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
export const LIST_EMPTY_STATE = '.grok-empty-state'; 

export const BALLOON_CONTAINER = '.d4-balloon-container';
export const ERROR_BALLOON = `${BALLOON_CONTAINER} .d4-balloon-error, ${BALLOON_CONTAINER} [class*="error"]`; 
export const INFO_BALLOON = `${BALLOON_CONTAINER} .d4-balloon-info, ${BALLOON_CONTAINER} [class*="info"]`; 

export const TOOLBOX_TAB_PANE = '[name="tab-pane-Toolbox"]';
export const TOOLBOX_PANE_SEARCH = '[name="pane-Search"]';
export const TOOLBOX_PANE_VIEWERS = '[name="pane-Viewers"]';
export const TOOLBOX_PANE_LAYOUTS = '[name="pane-Layouts"]';
export const TOOLBOX_PIN_TOOLBOX = '[aria-label="Pin toolbox"]';

export const APP_LOADED = RIBBON; 
export const PROGRESS_INDICATOR = '.grok-progress'; 

function escapeRegExp(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
