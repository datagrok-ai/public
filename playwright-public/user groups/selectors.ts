import { Page, Locator } from '@playwright/test';

// =============================================================================
// Users / Groups / Roles selectors — source of truth for tests in
// e2e/user groups/. Captured live on https://dev.datagrok.ai (2026-06-11).
// Browse-level selectors (tree, context menu, context panel, dialogs) are
// reused from ../browse/selectors.
// =============================================================================

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

// ---------- Gallery (Users / Groups / Roles list view) ------------------------

/** The card grid that holds the entity cards. */
export const GALLERY_GRID = '.grok-gallery-grid';
/** "shown / total" counter in the gallery header. */
export const GALLERY_COUNTS = '.grok-items-view-counts';

/** A single entity card label by display name (e.g. a user login or group/role name). */
export const galleryCardByName = (page: Page, name: string): Locator =>
  page.locator(`${GALLERY_GRID} .d4-link-label`, { hasText: new RegExp(`^${escapeRegExp(name)}$`) }).first();

/** Same card by its stable `name="span-<Name>"` attribute (spaces -> hyphens). */
export const galleryCardByAttr = (page: Page, name: string): Locator =>
  page.locator(`${GALLERY_GRID} [name="span-${name.replace(/ /g, '-')}"]`).first();

// ---------- Gallery header controls -------------------------------------------

/** Users view "NEW" combo button (dropdown: User... / Service User... / Invite a Friend...). */
export const NEW_USER_BUTTON = 'button[name="button-New"]';

/** Ribbon button by its text — used for "NEW GROUP..." and "New Role...". */
export const ribbonButtonByText = (page: Page, text: string): Locator =>
  page.locator('.d4-ribbon .ui-btn, .d4-ribbon button', { hasText: new RegExp(escapeRegExp(text), 'i') }).first();

/** Gallery search input, scoped by its placeholder prefix ("Search users"/"groups"/"roles"). */
export const gallerySearch = (page: Page, kind: 'users' | 'groups' | 'roles'): Locator =>
  page.locator(`input[placeholder^="Search ${kind}"]`).first();

// View-mode toggles + sort + filter. Scoped to the visible gallery toolbar containers —
// the same `name=` icons also exist hidden elsewhere, so an unscoped `.first()` is flaky.
export const VIEW_TOGGLE_BRIEF = '.grok-items-view-toggle [name="icon-grip-lines"]';
export const VIEW_TOGGLE_CARD = '.grok-items-view-toggle [name="icon-grip-horizontal"]';
export const VIEW_TOGGLE_GRID = '.grok-items-view-toggle [name="icon-table"]';
export const SORT_LIST = '.grok-items-view-toggle [name="icon-sort-alt"]';
export const TOGGLE_FILTERS = '.d4-search-bar [name="icon-filter"]';
export const REFRESH_LIST = '.d4-search-bar [name="icon-sync"], .grok-gallery-search-bar [name="icon-sync"]';

// ---------- Dialogs (NEW User / Group / Role; members editor) -----------------

export const DIALOG = '.d4-dialog';
export const DIALOG_TITLE = '.d4-dialog-title';
/** Dialog input field by its caption ("Email", "Login", "Name", "Description", ...). */
export const dialogInput = (page: Page, caption: string): Locator =>
  page.locator(`${DIALOG} [name="input-host-${caption.replace(/ /g, '-')}"] input, ${DIALOG} [name="input-host-${caption.replace(/ /g, '-')}"] textarea`).first();
/** Dialog button by text (OK / CANCEL / SAVE). */
export const dialogButton = (page: Page, text: string): Locator =>
  page.locator(`${DIALOG} .ui-btn`, { hasText: new RegExp(`^${escapeRegExp(text)}$`, 'i') }).first();

/** "Search by name or email to add..." input inside the members / assignees editor. */
export const MEMBERS_ADD_INPUT = `${DIALOG} input[placeholder^="Search by name or email"]`;

// ---------- Context-panel "MANAGE" entry points -------------------------------

/** The MANAGE link inside a context-panel pane (Members for groups, Assigned to for roles). */
export const manageButtonInPane = (page: Page, paneName: string): Locator =>
  page.locator(`${'.grok-prop-panel'} .d4-accordion-pane[name="pane-${paneName}"]`)
    .getByText('MANAGE', { exact: true }).first();

// ---------- Utils -------------------------------------------------------------

export function escapeRegExp(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
