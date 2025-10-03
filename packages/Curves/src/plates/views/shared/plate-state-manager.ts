/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Subject, Observable} from 'rxjs';
import {PlateFile, TemplateState, PlateStateChangeEvent} from './types';
import {createNewPlateForTemplate, PlateTemplate, PlateType} from '../../plates-crud';
import {parsePlates} from '../../../plate/csv-plates';
import {Plate} from '../../../plate/plate';
import {MAPPING_SCOPES} from './scopes';

export class PlateStateManager {
  private templateStates = new Map<number, TemplateState>();
  private _currentTemplate: PlateTemplate;
  private _currentPlateType: PlateType;
  private _sourceDataFrame: DG.DataFrame | null = null;

  private stateChange$ = new Subject<PlateStateChangeEvent>();
  public identifierColumns: (string | null)[] = [];

  constructor(
    initialTemplate: PlateTemplate,
    initialPlateType: PlateType
  ) {
    this._currentTemplate = initialTemplate;
    this._currentPlateType = initialPlateType;
  }


  get onStateChange(): Observable<PlateStateChangeEvent> {
    return this.stateChange$.asObservable();
  }

  get currentTemplate(): PlateTemplate {
    return this._currentTemplate;
  }

  get currentPlateType(): PlateType {
    return this._currentPlateType;
  }

  get sourceDataFrame(): DG.DataFrame | null {
    return this._sourceDataFrame;
  }

  get currentState(): TemplateState | undefined {
    return this.templateStates.get(this._currentTemplate.id);
  }

  get activePlate(): PlateFile | undefined {
    const state = this.currentState;
    if (!state || state.activePlateIdx < 0 || state.activePlateIdx >= state.plates.length)
      return undefined;
    return state.plates[state.activePlateIdx];
  }

  async setTemplate(template: PlateTemplate): Promise<void> {
    this._currentTemplate = template;
    if (!this.templateStates.has(template.id)) {
      this.templateStates.set(template.id, {
        plates: [],
        activePlateIdx: -1,
      });
    }
    this.stateChange$.next({
      type: 'template-changed',
      templateId: template.id,
    });
  }

  setPlateType(plateType: PlateType): void {
    this._currentPlateType = plateType;
  }
  private autodetectIdentifierColumn(): void {
    if (!this._sourceDataFrame) {
      this.identifierColumns = [];
      return;
    }

    const commonNames = [
      'barcode', 'plate barcode', 'plate id', 'plate index', 'plate number',
      'destination plate barcode'
    ];
    const dfCols = this._sourceDataFrame.columns.names();

    for (const name of commonNames) {
      const foundCol = dfCols.find((c) => c.toLowerCase().replace(/[\s_]/g, '') === name.replace(/\s/g, ''));
      if (foundCol) {
        this.identifierColumns = [foundCol];
        return;
      }
    }
    this.identifierColumns = [null];
  }

  async loadDataFrame(df: DG.DataFrame): Promise<void> {
    this._sourceDataFrame = df;
    this.autodetectIdentifierColumn();
    await this.reprocessPlates();
  }

  async reprocessPlates(): Promise<void> {
    if (!this._sourceDataFrame) return;

    const parsedPlateResults = parsePlates(
      this._sourceDataFrame,
      this.identifierColumns
    );

    const currentState = this.templateStates.get(this._currentTemplate.id) ?? {
      plates: [],
      activePlateIdx: -1,
    };

    currentState.plates = [];
    for (const parsedResult of parsedPlateResults) {
      const dummyFile = {
        name: `Plate: ${parsedResult.plate.barcode}`,
        fullPath: '',
      } as DG.FileInfo;
      const plate = parsedResult.plate;
      const sourceCols = new Set(plate.data.columns.names());
      const templateOnlyProps = this.currentTemplate.wellProperties
        .map((p) => p.name)
        .filter((p): p is string => !!p);

      for (const prop of templateOnlyProps) {
        if (sourceCols.has(prop))
          plate.addScopedAlias(MAPPING_SCOPES.TEMPLATE, prop, prop);
      }

      currentState.plates.push({
        plate: plate,
        file: dummyFile,
        commonProperties: parsedResult.commonProperties,
      });
    }

    const newActiveIndex = currentState.plates.length > 0 ? 0 : -1;
    currentState.activePlateIdx = newActiveIndex;
    this.templateStates.set(this._currentTemplate.id, currentState);

    this.stateChange$.next({
      type: 'plate-selected',
      plateIndex: newActiveIndex,
      plate: newActiveIndex === -1 ? undefined : currentState.plates[newActiveIndex],
    });
  }
  public addIdentifierColumn(): void {
    this.identifierColumns.push(null);
    this.reprocessPlates();
    this.stateChange$.next({type: 'identifier-changed'});
  }

  public removeIdentifierColumn(index: number): void {
    if (index >= 0 && index < this.identifierColumns.length) {
      this.identifierColumns.splice(index, 1);
      this.reprocessPlates();
      this.stateChange$.next({type: 'identifier-changed'});
    }
  }

  public setIdentifierColumn(index: number, columnName: string | null): void {
    if (index >= 0 && index < this.identifierColumns.length) {
      this.identifierColumns[index] = columnName;
      this.reprocessPlates();
    }
  }

  selectPlate(index: number): void {
    const state = this.currentState;
    if (!state || index === state.activePlateIdx) return;
    if (index >= 0 && index < state.plates.length) {
      state.activePlateIdx = index;
      this.stateChange$.next({
        type: 'plate-selected',
        plateIndex: index,
        plate: state.plates[index],
      });
    } else {
      state.activePlateIdx = -1;
      this.stateChange$.next({
        type: 'plate-selected',
        plateIndex: -1,
        plate: undefined,
      });
    }
  }

  removePlate(index: number): void {
    const state = this.currentState;
    if (!state) return;

    state.plates.splice(index, 1);
    if (state.activePlateIdx >= index)
      state.activePlateIdx = Math.max(0, state.activePlateIdx - 1);

    if (state.plates.length === 0)
      this.templateStates.delete(this._currentTemplate.id);

    this.stateChange$.next({
      type: 'plate-removed',
      plateIndex: index,
    });
  }
  public notifyPlateDataChanged(): void {
    const state = this.currentState;
    if (!state) return;
    console.log('PlateStateManager: Firing "plate-data-changed" event.');
    this.stateChange$.next({
      type: 'plate-data-changed',
      plateIndex: state.activePlateIdx,
      plate: this.activePlate,
    });
  }
  public remapScopedProperty(plateIndex: number, scope: string, target: string, source: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;
    const plateFile = state.plates[plateIndex];
    plateFile.plate.addScopedAlias(scope, source, target);

    this.stateChange$.next({
      type: scope === MAPPING_SCOPES.TEMPLATE ? 'mapping-changed' : 'analysis-mapping-changed',
      plateIndex,
      analysisType: scope as any,
      plate: plateFile,
    });
  }

  public undoScopedMapping(plateIndex: number, scope: string, targetProperty: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;

    const plateFile = state.plates[plateIndex];
    plateFile.plate.removeScopedAlias(scope, targetProperty);

    this.stateChange$.next({
      type: scope === MAPPING_SCOPES.TEMPLATE ? 'mapping-changed' : 'analysis-mapping-changed',
      plateIndex,
      analysisType: scope as any,
      plate: plateFile,
    });
  }

  public getScopedMapping(plateIndex: number, scope: string): Map<string, string> {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return new Map();

    return state.plates[plateIndex].plate.getScopedAliases(scope);
  }

  async getOrCreateDefaultPlate(): Promise<Plate> {
    const active = this.activePlate;
    if (active)
      return active.plate;
    return await createNewPlateForTemplate(this._currentPlateType, this._currentTemplate);
  }

  destroy(): void {
    this.stateChange$.complete();
  }
}
