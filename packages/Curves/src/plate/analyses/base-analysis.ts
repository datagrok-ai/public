/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../plate';
import {BaseAnalysisView} from './base-analysis-view';
import {
  AnalysisQuery,
  createAnalysisRun, getOrCreateProperty, PlateProperty,
  queryAnalysesGeneric,
  saveAnalysisResult, saveAnalysisRunParameter
} from '../../plates/plates-crud';
import {AnalysisRequiredFields} from '../../plates/views/components/analysis-mapping/analysis-mapping-panel';
import {_package} from '../../package';
import { PlateWidget } from '../plate-widget/plate-widget';

export interface IAnalysisProperty {
    name: string;
    type: DG.TYPE;
    category?: 'Parameter' | 'Output';
    choices?: string[];
    defaultValue?: any;
}

export interface IPlateAnalysis {
    readonly name: string;
    readonly friendlyName: string;
    parameters: IAnalysisProperty[];
    outputs: IAnalysisProperty[];
    getRequiredFields(): AnalysisRequiredFields[];
    createView(
        plate: Plate,
        plateWidget: PlateWidget,
        currentMappings: Map<string, string>,
        onMap: (target: string, source: string) => void,
        onUndo: (target: string) => void,
        onRerender?: () => void
    ): HTMLElement;
    queryResults(query: AnalysisQuery): Promise<DG.DataFrame>;
    formatResultsForGrid(rawResults: DG.DataFrame): DG.DataFrame;
    getSearchableProperties(): IAnalysisProperty[];
}

export abstract class AnalysisBase implements IPlateAnalysis {
    abstract readonly name: string;
    abstract readonly friendlyName: string;

    parameters: IAnalysisProperty[] = [];
    outputs: IAnalysisProperty[] = [];

    protected _parameterProperties = new Map<string, PlateProperty>();
    protected _outputProperties = new Map<string, PlateProperty>();


    async registerProperties(): Promise<void> {
      for (const param of this.parameters)
        this._parameterProperties.set(param.name, await getOrCreateProperty(param.name, param.type));
      for (const output of this.outputs)
        this._outputProperties.set(output.name, await getOrCreateProperty(output.name, output.type));
    }

    async queryResults(query: AnalysisQuery): Promise<DG.DataFrame> {
      // Use the generic query but then format results
      const rawResults = await queryAnalysesGeneric(query);
      return this.formatResultsForGrid(rawResults);
    }

    abstract formatResultsForGrid(rawResults: DG.DataFrame): DG.DataFrame;

    getSearchableProperties(): IAnalysisProperty[] {
      // By default, return all output properties that can be searched
      return this.outputs.filter((o) =>
        o.type === DG.TYPE.FLOAT ||
            o.type === DG.TYPE.INT ||
            o.type === DG.TYPE.STRING
      );
    }

    abstract getRequiredFields(): AnalysisRequiredFields[];
    abstract createView(plate: Plate, plateWidget: PlateWidget, currentMappings: Map<string, string>, onMap: (t: string, s: string) => void, onUndo: (t: string) => void, onRerender?: () => void): HTMLElement;
    protected abstract _getGroups(resultsDf: DG.DataFrame): { groupColumn: string, groups: string[] };
    protected _createStandardMappingView(
      plate: Plate,
      currentMappings: Map<string, string>,
      onMap: (t: string, s: string) => void, onUndo: (t: string) => void,
      createResultsView: (p: Plate, m: Map<string, string>) => HTMLElement | null
    ): HTMLElement {
      return new BaseAnalysisView(plate, {
        analysisName: this.friendlyName,
        requiredFields: this.getRequiredFields(),
        createResultsView: createResultsView
      }, currentMappings, onMap, onUndo).getRoot();
    }

    async saveResults(
      plate: Plate,
      resultsDf: DG.DataFrame,
      params: Record<string, any>,
      currentMappings: Map<string, string>
    ): Promise<void> {
      if (!plate.id) {
        grok.shell.warning('Please use the main "CREATE" button to save the plate before saving analysis results.');
        throw new Error('Cannot save analysis results for an unsaved plate.');
      }

      const pi = DG.TaskBarProgressIndicator.create(`Saving ${this.friendlyName} results...`);
      let runId: number | null = null;

      try {
        const {groupColumn, groups} = this._getGroups(resultsDf);
        runId = await createAnalysisRun(plate.id, this.name, groups);

        const savePromises: Promise<any>[] = [];

        for (const [name, value] of Object.entries(params)) {
          const prop = this._parameterProperties.get(name);
          if (prop) {
            savePromises.push(saveAnalysisRunParameter({
              runId: runId,
              propertyName: prop.name,
              propertyType: prop.type as DG.TYPE,
              value: value,
            }));
          }
        }

        for (const [requiredField, actualColumn] of currentMappings.entries()) {
          const mappingPropName = `Mapping: ${requiredField}`;
          savePromises.push(saveAnalysisRunParameter({
            runId: runId,
            propertyName: mappingPropName,
            propertyType: DG.TYPE.STRING,
            value: actualColumn,
          }));
        }

        for (const row of resultsDf.rows) {
          const groupKey = row.get(groupColumn);
          if (!groupKey) continue;
          const groupCombination = Array.isArray(groupKey) ? groupKey : [groupKey];

          for (const output of this.outputs) {
            const prop = this._outputProperties.get(output.name);
            if (prop && resultsDf.columns.contains(output.name)) {
              const value = row.get(output.name);
              if (value !== null && value !== undefined)
                savePromises.push(saveAnalysisResult({runId, propertyId: prop.id, propertyName: prop.name, propertyType: prop.type, value, groupCombination}));
            }
          }
        }

        await Promise.all(savePromises);

        grok.shell.info(`Saved ${this.friendlyName} results for plate ${plate.barcode}.`);
      } catch (e: any) {
        const errorMessage = `Failed to save ${this.friendlyName} results. Reason: ${e.message}`;
        _package.logger.debug(errorMessage);
        if (runId) {
          try {
            await grok.data.db.query('Curves:Plates', `DELETE FROM plates.analysis_runs WHERE id = ${runId};`);
          } catch (cleanupError: any) {
            const cleanupMessage = `CRITICAL: Failed to clean up analysis run ${runId} after a save error: ${cleanupError.message}`;
            grok.shell.error(cleanupMessage);
          }
        }
        throw new Error(errorMessage, {cause: e});
      } finally {
        pi.close();
      }
    }
}
