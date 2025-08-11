import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Subject, Observable} from 'rxjs';
import {PlateFile, TemplateState, PlateStateChangeEvent} from './types';
import {createNewPlateForTemplate, PlateTemplate, PlateType} from '../../plates-crud';
import {parsePlates} from '../../../plate/csv-plates';
import {Plate} from '../../../plate/plate';

export class PlateStateManager {
  private templateStates = new Map<number, TemplateState>();
  private _currentTemplate: PlateTemplate;
  private _currentPlateType: PlateType;
  private _sourceDataFrame: DG.DataFrame | null = null;

  private stateChange$ = new Subject<PlateStateChangeEvent>();

  public plateIdentifierColumn: string | null = 'Destination Plate Barcode';
  public replicateIdentifierColumn: string | null = 'Technical Duplicate ID';
  public plateNumberColumn: string | null = 'Plate Number';

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
    return state ? state.plates[state.activePlateIdx] : undefined;
  }

  async setTemplate(template: PlateTemplate): Promise<void> {
    this._currentTemplate = template;

    // Initialize empty state if needed
    if (!this.templateStates.has(template.id)) {
      this.templateStates.set(template.id, {
        plates: [],
        activePlateIdx: -1
      });
    }

    this.stateChange$.next({
      type: 'template-changed',
      templateId: template.id
    });
  }

  setPlateType(plateType: PlateType): void {
    this._currentPlateType = plateType;
  }

  async loadDataFrame(df: DG.DataFrame): Promise<void> {
    this._sourceDataFrame = df;
    await this.reprocessPlates();
  }

  async reprocessPlates(): Promise<void> {
    if (!this._sourceDataFrame) return;

    const parsedPlates = parsePlates(
      this._sourceDataFrame,
      this.plateIdentifierColumn,
      this.replicateIdentifierColumn,
      this.plateNumberColumn
    );

    const currentState = this.templateStates.get(this._currentTemplate.id) ?? {
      plates: [],
      activePlateIdx: -1
    };

    currentState.plates = [];
    for (const plate of parsedPlates) {
      const dummyFile = {
        name: `Plate: ${plate.barcode}`,
        fullPath: ''
      } as DG.FileInfo;

      currentState.plates.push({
        plate,
        file: dummyFile,
        reconciliationMap: new Map<string, string>()
      });
    }

    currentState.activePlateIdx = currentState.plates.length > 0 ? 0 : -1;
    this.templateStates.set(this._currentTemplate.id, currentState);

    this.stateChange$.next({type: 'plate-added'});
  }

  selectPlate(index: number): void {
    const state = this.currentState;
    if (state && index >= 0 && index < state.plates.length) {
      state.activePlateIdx = index;
      this.stateChange$.next({
        type: 'plate-selected',
        plateIndex: index,
        plate: state.plates[index]
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
      plateIndex: index
    });
  }

  applyMapping(plateIndex: number, mappedField: string, originalField: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;

    const plateFile = state.plates[plateIndex];
    const column = plateFile.plate.data.col(originalField);
    if (!column) {
      console.warn(`Column '${originalField}' not found for mapping`);
      return;
    }

    column.name = mappedField;
    plateFile.reconciliationMap.set(mappedField, originalField);

    grok.shell.info(`Mapped column '${originalField}' to '${mappedField}'`);

    this.stateChange$.next({
      type: 'mapping-changed',
      plateIndex,
      plate: plateFile
    });
  }

  undoMapping(plateIndex: number, mappedField: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;

    const plateFile = state.plates[plateIndex];
    const originalName = plateFile.reconciliationMap.get(mappedField);

    if (originalName) {
      const column = plateFile.plate.data.col(mappedField);
      if (column) {
        column.name = originalName;
        plateFile.reconciliationMap.delete(mappedField);

        grok.shell.info(`Reverted mapping for '${mappedField}'`);

        this.stateChange$.next({
          type: 'mapping-changed',
          plateIndex,
          plate: plateFile
        });
      }
    }
  }

  syncMappings(sourceMappings: Map<string, string>, selectedIndexes: number[]): void {
    const state = this.currentState;
    if (!state) return;

    const selectedSet = new Set(selectedIndexes);
    let modified = 0;

    state.plates.forEach((plateFile, index) => {
      const shouldBeMapped = selectedSet.has(index);
      const isCurrentlyMapped = [...sourceMappings.keys()].every(
        (key) => plateFile.reconciliationMap.has(key)
      );

      if (shouldBeMapped && !isCurrentlyMapped) {
        // Apply mapping
        sourceMappings.forEach((oldName, newName) => {
          const column = plateFile.plate.data.col(oldName);
          if (column && oldName !== newName) {
            column.name = newName;
            plateFile.reconciliationMap.set(newName, oldName);
          }
        });
        modified++;
      } else if (!shouldBeMapped && isCurrentlyMapped) {
        // Undo mapping
        sourceMappings.forEach((oldName, newName) => {
          const column = plateFile.plate.data.col(newName);
          if (column && plateFile.reconciliationMap.get(newName) === oldName) {
            column.name = oldName;
            plateFile.reconciliationMap.delete(newName);
          }
        });
        modified++;
      }
    });

    if (modified > 0) {
      grok.shell.info(`Synchronized mappings for ${modified} plate(s)`);
      this.stateChange$.next({type: 'mapping-changed'});
    }
  }

  async getOrCreateDefaultPlate(): Promise<Plate> {
    const state = this.currentState;
    if (state && state.plates.length > 0)
      return state.plates[state.activePlateIdx].plate;

    return await createNewPlateForTemplate(this._currentPlateType, this._currentTemplate);
  }

  destroy(): void {
    this.stateChange$.complete();
  }
}
