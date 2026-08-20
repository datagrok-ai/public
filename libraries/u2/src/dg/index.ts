/* The ONLY layer allowed to import datagrok-api (enforced by eslint). Bridges u2's
   platform-free core onto the platform's existing lifecycle and event surfaces. */
export {asDartInput, dartInputFor, DartInput} from './input-bridge.js';
export type {DartInputOptions, PropertyInputBuilder} from './input-bridge.js';
export {columnInput, tableInput, tablesInput, ColumnPicker} from './pickers.js';
export type {ColumnInputOptions, ColumnPickerOptions} from './pickers.js';
export {ObjectForm, propertyForm, objectForm, inputForProperty, PlatformInputs} from './object-form.js';
export type {PropertyLike} from '../core/property-like.js';
export type {PropertySource, ObjectFormOptions, FieldOverride,
  PropertyInputOptions, InputFactory} from './object-form.js';
export {propertyEditor, PropertyEditor} from './property-editor.js';
export type {PropertyEditorOptions} from './property-editor.js';
export {fromDartInput, PlatformInput} from './from-dart-input.js';
export type {DartInputLike} from './from-dart-input.js';
export {userInput} from './user-input.js';
export {dapiSource, dapiPager, sanitizeFilterValue} from './dapi-source.js';
export type {DapiSourceLike, DapiSourceOptions, DapiPagerSourceLike, DapiPagerOptions} from './dapi-source.js';
export {handlerRenderer, HandlerRenderer, chip, EntityChip, entityCard, EntityCard, entityInput}
  from './entity.js';
export type {ChipOptions, EntityInputOptions} from './entity.js';
export {Editors} from './editors.js';
export type {EditorRule} from './editors.js';
export {moleculeInput, moleculeRenderer} from './molecule.js';
export type {MoleculeRendererOptions} from './molecule.js';
export {metaInput} from './meta-input.js';
export type {MetaInputOptions} from './meta-input.js';
export {fileInput, FileInput} from './file-input.js';
export type {FileInputOptions, FileInputMode} from './file-input.js';
export {filesInput, FilesInput} from './files-input.js';
export type {FilesInputOptions} from './files-input.js';
export {rsaInput, RsaInput} from './rsa-input.js';
export type {RsaInputOptions} from './rsa-input.js';
export {columnRenderer, ColumnRenderer} from './column-renderer.js';
export {columnsInput, ColumnsInput, columnsMapInput, ColumnsMapInput, aggregatedColumnsInput,
  AggregatedColumnsInput, aggregationsFor, defaultAggregation} from './columns.js';
export type {ColumnsInputOptions, ColumnsMapInputOptions, ColumnKey, ColumnAggregation,
  AggregatedColumnsInputOptions} from './columns.js';
export {host, U2Widget} from './widget-host.js';
export {appView} from './app-view.js';
export type {AppViewOptions} from './app-view.js';
export {SpecNodeRef, specTree, brokenCount, nodeLabel, idPath} from './designer/node-ref.js';
export type {SpecTree} from './designer/node-ref.js';
export {registerSpecNodeHandler} from './designer/handler.js';
export {designerView, SpecDesigner} from './designer/view.js';
export type {DesignerViewOptions} from './designer/view.js';

import * as DG from 'datagrok-api/dg';
import {Observable} from 'rxjs';
import {Scope} from '../core/scope.js';
import {signal, ReadonlySignal, rawEffect} from '../core/signals.js';
import {PlatformInputs} from './object-form.js';
import {FileInput} from './file-input.js';

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
