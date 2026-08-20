/* The ONLY layer allowed to import datagrok-api (enforced by eslint). Bridges u2's
   platform-free core onto the platform's existing lifecycle and event surfaces. */
export {asDartInput, dartInputFor, DartInput} from './inputs/input-bridge.js';
export type {DartInputOptions, PropertyInputBuilder} from './inputs/input-bridge.js';
export {columnInput, tableInput, tablesInput, ColumnPicker} from './inputs/pickers.js';
export type {ColumnInputOptions, ColumnPickerOptions} from './inputs/pickers.js';
export {ObjectForm, propertyForm, objectForm, inputForProperty, PlatformInputs} from './forms/object-form.js';
export type {PropertyLike} from '../core/property-like.js';
export type {PropertySource, ObjectFormOptions, FieldOverride,
  PropertyInputOptions, InputFactory} from './forms/object-form.js';
export {propertyEditor, PropertyEditor} from './forms/property-editor.js';
export type {PropertyEditorOptions} from './forms/property-editor.js';
export {fromDartInput, PlatformInput} from './inputs/from-dart-input.js';
export type {DartInputLike} from './inputs/from-dart-input.js';
export {userInput} from './inputs/user-input.js';
export {dapiSource, dapiPager, sanitizeFilterValue} from './entities/dapi-source.js';
export type {DapiSourceLike, DapiSourceOptions, DapiPagerSourceLike, DapiPagerOptions} from './entities/dapi-source.js';
export {handlerRenderer, HandlerRenderer, chip, EntityChip, entityCard, EntityCard, entityInput}
  from './entities/entity.js';
export type {ChipOptions, EntityInputOptions} from './entities/entity.js';
export {Editors} from './forms/editors.js';
export type {EditorRule} from './forms/editors.js';
export {moleculeInput, moleculeRenderer} from './inputs/molecule.js';
export type {MoleculeRendererOptions} from './inputs/molecule.js';
export {metaInput} from './inputs/meta-input.js';
export type {MetaInputOptions} from './inputs/meta-input.js';
export {fileInput, FileInput} from './inputs/file-input.js';
export type {FileInputOptions, FileInputMode} from './inputs/file-input.js';
export {filesInput, FilesInput} from './inputs/files-input.js';
export type {FilesInputOptions} from './inputs/files-input.js';
export {rsaInput, RsaInput} from './inputs/rsa-input.js';
export type {RsaInputOptions} from './inputs/rsa-input.js';
export {columnRenderer, ColumnRenderer} from './entities/column-renderer.js';
export {columnsInput, ColumnsInput, columnsMapInput, ColumnsMapInput, aggregatedColumnsInput,
  AggregatedColumnsInput, aggregationsFor, defaultAggregation} from './inputs/columns.js';
export type {ColumnsInputOptions, ColumnsMapInputOptions, ColumnKey, ColumnAggregation,
  AggregatedColumnsInputOptions} from './inputs/columns.js';
export {host, U2Widget} from './shell/widget-host.js';
export {appView} from './shell/app-view.js';
export type {AppViewOptions} from './shell/app-view.js';
export {SpecNodeRef, SpecNodesRef, specTree, brokenCount, nodeLabel, idPath} from './designer/node-ref.js';
export type {SpecTree} from './designer/node-ref.js';
export {registerSpecNodeHandler} from './designer/handler.js';
export {Palette} from './designer/palette.js';
export {Tray, sourceNode, funcSourceNode} from './designer/tray.js';
export type {TrayOptions} from './designer/tray.js';
export {accepts, resolveDrop} from './designer/dnd.js';
export type {DropRect, DropTarget} from './designer/dnd.js';
export {makeDesignerDroppable, readDrop, dropNode, funcRef, tabularExtensions, OPEN_FILE}
  from './designer/drop.js';
export type {DropItem, DropReading, DesignerDropOptions} from './designer/drop.js';
export {designerView, SpecDesigner} from './designer/view.js';
export type {DesignerViewOptions} from './designer/view.js';
export {SAMPLES} from './designer/samples.js';
export {loadGallery, saveToGallery, listGallery, GALLERY_KEY} from './designer/gallery.js';
export {bindTree} from './designer/bind-model.js';
export type {BindTreeNode} from './designer/bind-model.js';
export {bindPicker, bindGroups, bindRows} from './designer/bind-picker.js';
export type {BindGroup, BindRows} from './designer/bind-picker.js';
export {funcPicker, funcEntries, filterFuncs, paramProps, paramValues, eventEntry, eventPick}
  from './designer/func-picker.js';
export type {FuncEntry, FuncLike, FuncPick, FuncPickerOptions} from './designer/func-picker.js';
export {sourceStatus, statusText, refreshSource} from './designer/source-status.js';
export type {SourceStatus} from './designer/source-status.js';
export {platformContext} from './shell/spec-context.js';
// side-effect only: filling `backends` is what makes the data sources work in the platform
import './shell/source-backends.js';

import * as DG from 'datagrok-api/dg';
import {Observable} from 'rxjs';
import {Scope} from '../core/scope.js';
import {signal, ReadonlySignal, rawEffect} from '../core/signals.js';
import {PlatformInputs} from './forms/object-form.js';
import {FileInput} from './inputs/file-input.js';

// the schema-driven router keeps loading without the platform, so the editors that need it are
// wired here instead of imported there
PlatformInputs.register('file', (prop, options) => new FileInput(options));

/** Mirrors an rxjs observable into a signal owned by `scope`. */
export function toSignal<T>(observable: Observable<T>, initial: T, scope: Scope): ReadonlySignal<T> {
  const s = signal(initial);
  const sub = observable.subscribe((v) => s.value = v);
  scope.own(() => sub.unsubscribe());
  return s;
}

/** Exposes a signal as an rxjs observable (emits current value on subscribe).
 * The dual event surface: new code binds signals, old code subscribes — no migration. */
export function toObservable<T>(source: ReadonlySignal<T>): Observable<T> {
  return new Observable<T>((observer) => rawEffect(() => observer.next(source.value)));
}

/** Live u2 scopes vs platform-registered widgets — the leak detector's report.
 * After a view with u2 content closes, both numbers must return to their baseline. */
export function leakReport(): {liveScopes: number, liveWidgets: number} {
  return {liveScopes: Scope.liveCount, liveWidgets: DG.Widget.getAll().length};
}
