/* eslint-disable max-len */
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

  public setMappings(plateIndex: number, newMappings: Map<string, string>): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;

    const plateFile = state.plates[plateIndex];
    // plateFile.reconciliationMap = new Map(newMappings);

    this.stateChange$.next({
      type: 'mapping-changed',
      plateIndex,
      plate: plateFile,
    });
  }
  // public setAnalysisMapping(plateIndex: number, analysisType: 'drc' | 'doseRatio', mappings: Map<string, string>): void {
  //   const state = this.currentState;
  //   if (!state || !state.plates[plateIndex]) return;

  //   const plateFile = state.plates[plateIndex];
  //   plateFile.analysisMappings[analysisType] = new Map(mappings);

  //   this.stateChange$.next({
  //     type: 'analysis-mapping-changed',
  //     plateIndex,
  //     analysisType,
  //     plate: plateFile,
  //   });
  // }

  // Add to PlateStateManager
  public getAnalysisMapping(plateIndex: number, analysisType: 'drc' | 'doseRatio'): Map<string, string> {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return new Map();

    const plate = state.plates[plateIndex].plate;
    const currentMappings = new Map<string, string>();

    // The plate itself is the source of truth.
    for (const layer of plate.getLayerNames()) {
      const aliases = plate.getAliases(layer);
      for (const alias of aliases)
        currentMappings.set(alias, layer);
    }
    return currentMappings;
  }


  public remapAnalysisProperty(plateIndex: number, analysisType: 'drc' | 'doseRatio', target: string, source: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;
    const plateFile = state.plates[plateIndex];

    plateFile.plate.addAlias(source, target);

    this.stateChange$.next({
      type: 'analysis-mapping-changed',
      plateIndex,
      analysisType, // Pass the analysisType so the specific view can update
      plate: plateFile,
    });
  }


  public undoAnalysisMapping(plateIndex: number, analysisType: 'drc' | 'doseRatio', targetProperty: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;

    const plateFile = state.plates[plateIndex];
    plateFile.plate.removeAlias(targetProperty);
    grok.shell.info(`Reverted mapping for '${targetProperty}' in ${analysisType} analysis`);

    this.stateChange$.next({
      type: 'analysis-mapping-changed',
      plateIndex,
      analysisType,
      plate: plateFile,
    });
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

  // public getActivePlateMappedData(): DG.DataFrame | null {
  //   const activePlate = this.activePlate;
  //   if (!activePlate || !activePlate.plate.data) return null;

  //   const sourceDf = activePlate.plate.data;
  //   const map = activePlate.reconciliationMap;

  //   console.log('[DEBUG] Mapping plate data with reconciliation map:', map);

  //   const newColumns: DG.Column[] = [];
  //   const processedSources = new Set<string>();

  //   map.forEach((sourceColName, targetColName) => {
  //     const sourceCol = sourceDf.col(sourceColName);
  //     if (sourceCol) {
  //       const newCol = sourceCol.clone();
  //       newCol.name = targetColName;
  //       newColumns.push(newCol);
  //       processedSources.add(sourceColName);
  //       console.log(`[DEBUG] Mapped column: ${sourceColName} -> ${targetColName}`);
  //     }
  //   });

  //   for (const sourceCol of sourceDf.columns) {
  //     if (!processedSources.has(sourceCol.name)) {
  //       const newCol = sourceCol.clone();
  //       newColumns.push(newCol);
  //       console.log(`[DEBUG] Added unmapped column: ${sourceCol.name}`);
  //     }
  //   }

  //   if (newColumns.length === 0) {
  //     console.log('[DEBUG] No columns to create mapped DataFrame');
  //     return null;
  //   }

  //   const mappedDf = DG.DataFrame.fromColumns(newColumns);
  //   mappedDf.name = sourceDf.name + '_mapped';

  //   console.log('[DEBUG] Final mapped columns:', mappedDf.columns.names());

  //   return mappedDf;
  // }

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

  async loadDataFrame(df: DG.DataFrame): Promise<void> {
    this._sourceDataFrame = df;
    await this.reprocessPlates();
  }

  async reprocessPlates(): Promise<void> {
    if (!this._sourceDataFrame) return;

    const parsedPlateResults = parsePlates(
      this._sourceDataFrame,
      this.plateIdentifierColumn,
      this.replicateIdentifierColumn,
      this.plateNumberColumn
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
          plate.addAlias(prop, prop);
      }

      // Auto-alias common analysis fields on load
      const commonAnalysisFields = ['Activity', 'Concentration', 'SampleID', 'Agonist_Concentration_M', 'Antagonist_Concentration_M', 'Percent_Inhibition'];
      for (const field of commonAnalysisFields) {
        if (sourceCols.has(field))
          plate.addAlias(field, field);
      }

      currentState.plates.push({
        plate: plate,
        file: dummyFile,
        // reconciliationMap: new Map(), // This is fully deprecated now
        // FIX #3: Remove the obsolete analysisMappings property
        // analysisMappings: {drc: new Map(), doseRatio: new Map()},
        commonProperties: parsedResult.commonProperties,
      });
    }

    currentState.activePlateIdx = currentState.plates.length > 0 ? 0 : -1;
    this.templateStates.set(this._currentTemplate.id, currentState);
    this.stateChange$.next({type: 'plate-added'});
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

  public remapProperty(plateIndex: number, targetProperty: string, newSourceColumn: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;
    const plateFile = state.plates[plateIndex];

    plateFile.plate.addAlias(newSourceColumn, targetProperty);

    grok.shell.info(`Mapped '${targetProperty}' to column '${newSourceColumn}'`);
    this.stateChange$.next({
      type: 'mapping-changed',
      plateIndex,
      plate: plateFile,
    });
  }

  undoMapping(plateIndex: number, targetProperty: string): void {
    const state = this.currentState;
    if (!state || !state.plates[plateIndex]) return;
    const plateFile = state.plates[plateIndex];

    plateFile.plate.removeAlias(targetProperty);
    grok.shell.info(`Reverted mapping for '${targetProperty}'`);
    this.stateChange$.next({
      type: 'mapping-changed',
      plateIndex,
      plate: plateFile,
    });
  }

  syncMappings(sourceMappings: Map<string, string>, selectedIndexes: number[]): void {
    const state = this.currentState;
    if (!state) return;
    const selectedSet = new Set(selectedIndexes);
    let modified = 0;
    state.plates.forEach((plateFile, index) => {
      if (selectedSet.has(index)) {
        for (const [targetProperty, sourceColumn] of sourceMappings.entries())
          plateFile.plate.addAlias(sourceColumn, targetProperty);

        modified++;
      }
    });

    if (modified > 0) {
      grok.shell.info(`Synchronized mappings for ${modified} plate(s)`);
      this.stateChange$.next({type: 'mapping-changed'});
    }
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
