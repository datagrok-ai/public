/**
 * Moved into the platform. The widget-action Func vocabulary now ships inside
 * `datagrok-api` (`src/ui/domains/domains-grid.ts`, as statics of
 * `DG.DomainGrid`) — this module keeps the import path and the flat names, and
 * the Funcs stay registered ONCE per session by qualified name.
 *
 * @module actions
 */

import {DomainGrid} from 'datagrok-api/dg';

export {domainActionFunc} from 'datagrok-api/dg';

export type {IWidgetActions} from 'datagrok-api/dg';

export const WIDGET_ACTION_NAMESPACE = DomainGrid.WIDGET_ACTION_NAMESPACE;
export const WIDGET_ACTION_TAGS = DomainGrid.WIDGET_ACTION_TAGS;
export const widgetActionFunc = DomainGrid.widgetActionFunc;
export const saveFunc = DomainGrid.saveFunc;
export const discardFunc = DomainGrid.discardFunc;
export const resetFunc = DomainGrid.resetFunc;
export const addRowFunc = DomainGrid.addRowFunc;
export const deleteRowFunc = DomainGrid.deleteRowFunc;
export const refreshFunc = DomainGrid.refreshFunc;
export const actionFuncName = DomainGrid.actionFuncName;
