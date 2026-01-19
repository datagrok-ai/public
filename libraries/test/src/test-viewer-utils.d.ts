import type * as _DG from 'datagrok-api/dg';
import { Observable } from 'rxjs';
export declare function createViewer(tv: _DG.TableView, v: string, packageName?: string): Promise<_DG.Viewer>;
export declare function testViewerInternal(tv: _DG.TableView, viewerName: string, packageName: string, event: Observable<any>, actions?: (args: any, withDelay: boolean) => Promise<any>, awaitViewer?: (viewer: _DG.Viewer) => Promise<void>, layout?: _DG.ViewLayout, actionArgs?: any, actionsWithDelay?: boolean): Promise<any>;
export declare function selectFilterChangeCurrent(args: any, withDelay?: boolean): Promise<void>;
export declare function filterAsync(args: any, withDelay?: boolean): Promise<void>;
export declare function changeOptionsSaveLayout(args: any, withDelay?: boolean): Promise<{
    layout: _DG.ViewLayout;
    savedProps: any;
}>;
export declare function loadLayout(args: any, withDelay?: boolean): Promise<void>;
//# sourceMappingURL=test-viewer-utils.d.ts.map