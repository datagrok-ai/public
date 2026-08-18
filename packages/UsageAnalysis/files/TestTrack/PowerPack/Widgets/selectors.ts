import { Page, Locator } from '@playwright/test';

export const WELCOME_VIEW = '.power-pack-welcome-view';
export const WIDGETS_PANEL = '.power-pack-widgets-panel';
export const WIDGETS_HOST = '.power-pack-widgets-host';
export const SEARCH_HOST = '.power-pack-search-host';

export const WIDGET_HOST = '.power-pack-widget-host';
export const WIDGET_TITLE = '.d4-dialog-title';
export const WIDGET_HEADER_HIDDEN = '.d4-dialog-header-hidden';
export const WIDGET_CONTENT = '.power-pack-widget-content';
export const WIDGET_CLOSE_ICON = 'i.grok-font-icon-close';

export const widgetByTitle = (page: Page, title: string): Locator =>
  page.locator(`${WIDGET_HOST}[widget-title="${title}"]`);

export const widgetCloseIcon = (page: Page, title: string): Locator =>
  widgetByTitle(page, title).locator(WIDGET_CLOSE_ICON);

export const widgetContent = (page: Page, title: string): Locator =>
  widgetByTitle(page, title).locator(WIDGET_CONTENT);

export const customizeLink = (page: Page): Locator =>
  page.getByText('Customize widgets...', { exact: false }).first();

const CONTEXT_PANEL = '.grok-prop-panel';
export const customizeToggle = (page: Page, label: string): Locator =>
  page.locator(`${CONTEXT_PANEL} .ui-input-root.ui-input-bool`)
    .filter({ hasText: label })
    .locator('input[type="checkbox"]');

export const SEARCH_INPUT = 'input.power-search-search-everywhere-input';

export const SPOTLIGHT_TAB = `${WIDGET_HOST}[widget-title="Spotlight"] .d4-tab-header`;
export const DEMO_OF_THE_DAY = 'Demo of the day';
export const TAB_SELECTED = /(^|\s)selected(\s|$)/;
export const MARK_ALL_AS_READ = 'Mark all as read';

export const spotlightTab = (page: Page, name: string): Locator =>
  page.locator(SPOTLIGHT_TAB).filter({ hasText: name }).first();

export const widgetLink = (page: Page, title: string, text: string): Locator =>
  widgetByTitle(page, title).locator('a.ui-link', { hasText: text }).first();

export const communityLinks = (page: Page): Locator =>
  widgetByTitle(page, 'Community').locator('a');

export const WIDGET_FUNC_KEY: Record<string, string> = {
  Spotlight: 'activityDashboardWidget',
  Community: 'communityWidget',
  Usage: 'usageWidget',
  Reports: 'reportsWidget',
};
