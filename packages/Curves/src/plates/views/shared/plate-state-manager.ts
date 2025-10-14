/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Subject} from 'rxjs';
import {PlateFile, TemplateState, PlateStateChangeEvent} from './types';
import {createNewPlateForTemplate, PlateTemplate, PlateType} from '../../plates-crud';
import {parsePlates} from '../../../plate/csv-plates';
import {Plate} from '../../../plate/plate';
import {MAPPING_SCOPES} from './scopes';

export class PlateStateManager {
  private templateStates = new Map<number, TemplateState>();
  private plateMappings = new Map<number, Map<string, Map<string, string>>>();
  public currentTemplate: PlateTemplate;
  public currentPlateType: PlateType;
  public sourceDataFrame: DG.DataFrame | null = null;
  public identifierColumns: (string | null)[] = [];

  public onStateChange$ = new Subject<PlateStateChangeEvent>();

  constructor(initialTemplate: PlateTemplate, initialPlateType: PlateType) {
    this.currentTemplate = initialTemplate;
    this.currentPlateType = initialPlateType;
  }

  get currentState(): TemplateState | undefined {
    return this.templateStates.get(this.currentTemplate.id);
  }

  get activePlate(): PlateFile | undefined {
    const state = this.currentState;
    if (!state || state.activePlateIdx < 0 || state.activePlateIdx >= state.plates.length)
      return undefined;
    return state.plates[state.activePlateIdx];
  }

  public setMapping(plateIndex: number, scope: string, target: string, source: string): void {
    if (!this.plateMappings.has(plateIndex))
      this.plateMappings.set(plateIndex, new Map());

    const plateMaps = this.plateMappings.get(plateIndex)!;
    if (!plateMaps.has(scope))
      plateMaps.set(scope, new Map());

    plateMaps.get(scope)!.set(target, source);

    this.onStateChange$.next({
      type: scope === MAPPING_SCOPES.TEMPLATE ? 'mapping-changed' : 'analysis-mapping-changed',
      plateIndex,
      analysisType: scope as any,
      plate: this.currentState?.plates[plateIndex],
    });
  }

  public removeMapping(plateIndex: number, scope: string, target: string): void {
    const plateMaps = this.plateMappings.get(plateIndex);
    if (plateMaps?.has(scope))
      plateMaps.get(scope)!.delete(target);

    this.onStateChange$.next({
      type: scope === MAPPING_SCOPES.TEMPLATE ? 'mapping-changed' : 'analysis-mapping-changed',
      plateIndex,
      analysisType: scope as any,
      plate: this.currentState?.plates[plateIndex],
    });
  }

  public getMappings(plateIndex: number, scope: string): Map<string, string> {
    const plateMaps = this.plateMappings.get(plateIndex);
    return plateMaps?.has(scope) ? new Map(plateMaps.get(scope)!) : new Map();
  }

  public resolveColumn(plateIndex: number, scope: string, nameOrAlias: string): DG.Column | null {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return null;

    const plate = state.plates[plateIndex].plate;
    const plateMaps = this.plateMappings.get(plateIndex);

    if (plateMaps?.has(scope)) {
      const scopeMap = plateMaps.get(scope)!;
      for (const [alias, columnName] of scopeMap) {
        if (alias === nameOrAlias)
          return plate.data.col(columnName);
      }
    }

    return plate.data.columns.contains(nameOrAlias) ? plate.data.col(nameOrAlias) : null;
  }

  async setTemplate(template: PlateTemplate): Promise<void> {
    this.currentTemplate = template;
    if (!this.templateStates.has(template.id)) {
      this.templateStates.set(template.id, {
        plates: [],
        activePlateIdx: -1,
      });
    }
    this.onStateChange$.next({
      type: 'template-changed',
      templateId: template.id,
    });
  }

  setPlateType(plateType: PlateType): void {
    this.currentPlateType = plateType;
  }

  private autodetectIdentifierColumn(): void {
    if (!this.sourceDataFrame) {
      this.identifierColumns = [];
      return;
    }
    const commonNames = ['barcode', 'plate barcode', 'plate id', 'plate index', 'plate number', 'destination plate barcode'];
    const dfCols = this.sourceDataFrame.columns.names();

    const foundCol = commonNames.map((name) =>
      dfCols.find((c) => c.toLowerCase().replace(/[\s_]/g, '') === name.replace(/\s/g, ''))
    ).find((col) => col);

    this.identifierColumns = foundCol ? [foundCol] : [null];
  }

  async loadDataFrame(df: DG.DataFrame): Promise<void> {
    this.sourceDataFrame = df;
    this.autodetectIdentifierColumn();
    await this.processPlates();
  }

  async processPlates(): Promise<void> {
    if (!this.sourceDataFrame) return;

    const parsedPlateResults = parsePlates(this.sourceDataFrame, this.identifierColumns);
    const currentState = this.currentState ?? {plates: [], activePlateIdx: -1};

    currentState.plates = parsedPlateResults.map((parsedResult) => ({
      plate: parsedResult.plate,
      file: {name: `Plate: ${parsedResult.plate.barcode}`, fullPath: ''} as DG.FileInfo,
      commonProperties: parsedResult.commonProperties,
    }));

    const newActiveIndex = currentState.plates.length > 0 ? 0 : -1;
    currentState.activePlateIdx = newActiveIndex;
    this.templateStates.set(this.currentTemplate.id, currentState);

    this.onStateChange$.next({
      type: 'plate-selected',
      plateIndex: newActiveIndex,
      plate: newActiveIndex === -1 ? undefined : currentState.plates[newActiveIndex],
    });
  }

  public addIdentifierColumn(): void {
    this.identifierColumns.push(null);
    this.processPlates();
    this.onStateChange$.next({type: 'identifier-changed'});
  }

  public removeIdentifierColumn(index: number): void {
    if (index >= 0 && index < this.identifierColumns.length) {
      this.identifierColumns.splice(index, 1);
      this.processPlates();
      this.onStateChange$.next({type: 'identifier-changed'});
    }
  }

  public setIdentifierColumn(index: number, columnName: string | null): void {
    if (index >= 0 && index < this.identifierColumns.length) {
      this.identifierColumns[index] = columnName;
      this.processPlates();
    }
  }

  selectPlate(index: number): void {
    const state = this.currentState;
    if (!state || index === state.activePlateIdx) return;

    if (index >= 0 && index < state.plates.length)
      state.activePlateIdx = index;
    else
      state.activePlateIdx = -1;


    this.onStateChange$.next({
      type: 'plate-selected',
      plateIndex: state.activePlateIdx,
      plate: state.activePlateIdx === -1 ? undefined : state.plates[state.activePlateIdx],
    });
  }

  public notifyPlateDataChanged(): void {
    if (!this.currentState) return;

    this.onStateChange$.next({
      type: 'plate-data-changed',
      plateIndex: this.currentState.activePlateIdx,
      plate: this.activePlate,
    });
  }

  async getOrCreateDefaultPlate(): Promise<Plate> {
    return this.activePlate?.plate ?? createNewPlateForTemplate(this.currentPlateType, this.currentTemplate);
  }

  destroy(): void {
    this.onStateChange$.complete();
  }
}
