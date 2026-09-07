import { Page, Locator } from '@playwright/test';

export {
  CONTEXT_MENU,
  CONTEXT_MENU_ITEM,
  CONTEXT_MENU_ITEM_LABEL,
  contextMenuItem,
  contextMenuItemByName,
  treeNodeByPath,
  CONTEXT_PANEL,
  CONTEXT_PANEL_STAR,
  CONTEXT_PANEL_EXPAND_ALL,
  infoPaneByName,
  RIBBON,
  BALLOON_CONTAINER,
} from '../browse/selectors';

export const GALLERY_GRID = '.grok-gallery-grid';

export const GALLERY_COUNTS = '.grok-items-view-counts';

export const galleryCardByName = (page: Page, name: string): Locator =>
  page.locator(`${GALLERY_GRID} .d4-link-label`, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

export const galleryCardByAttr = (page: Page, name: string): Locator =>
  page.locator(`${GALLERY_GRID} [name="span-${name.replace(/ /g, '-')}"]`).first();

export const NEW_USER_BUTTON = 'button[name="button-New"]';

export const ribbonButtonByText = (page: Page, text: string): Locator =>
  page.locator('.d4-ribbon .ui-btn, .d4-ribbon button', { hasText: new RegExp(escapeRegExp(text), 'i') }).first();

export const gallerySearch = (page: Page, kind: 'users' | 'groups' | 'roles'): Locator =>
  page.locator(`input[placeholder^="Search ${kind}"]`).first();

export const VIEW_TOGGLE_BRIEF = '.grok-items-view-toggle [name="icon-grip-lines"]';
export const VIEW_TOGGLE_CARD = '.grok-items-view-toggle [name="icon-grip-horizontal"]';
export const VIEW_TOGGLE_GRID = '.grok-items-view-toggle [name="icon-table"]';
export const SORT_LIST = '.grok-items-view-toggle [name="icon-sort-alt"]';
export const TOGGLE_FILTERS = '.d4-search-bar [name="icon-filter"]';
export const REFRESH_LIST = '.d4-search-bar [name="icon-sync"], .grok-gallery-search-bar [name="icon-sync"]';

export const DIALOG = '.d4-dialog';
export const DIALOG_TITLE = '.d4-dialog-title';

export const dialogInput = (page: Page, caption: string): Locator =>
  page.locator(`${DIALOG} [name="input-host-${caption.replace(/ /g, '-')}"] input, ${DIALOG} [name="input-host-${caption.replace(/ /g, '-')}"] textarea`).first();

export const dialogButton = (page: Page, text: string): Locator =>
  page.locator(`${DIALOG} .ui-btn`, { hasText: new RegExp(`^${escapeRegExp(text)}$`, 'i') }).first();

export const MEMBERS_ADD_INPUT = `${DIALOG} input[placeholder^="Search by name or email"]`;

export const manageButtonInPane = (page: Page, paneName: string): Locator =>
  page.locator(`${'.grok-prop-panel'} .d4-accordion-pane[name="pane-${paneName}"]`)
    .getByText('MANAGE', { exact: true }).first();

export function escapeRegExp(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
