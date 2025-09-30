/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../plate';
import {PlateWidget} from '../plate-widget';
import {BaseAnalysisView} from '../base-analysis-view';
import {
  createAnalysisRun, getOrCreateProperty, PlateProperty,
  saveAnalysisResult, saveAnalysisRunParameter
} from '../../plates/plates-crud';
import {AnalysisRequiredFields} from '../../plates/views/components/analysis-mapping/analysis-mapping-panel';

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
        onUndo: (target: string) => void
    ): HTMLElement;
}

export abstract class AbstractPlateAnalysis implements IPlateAnalysis {
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

    abstract getRequiredFields(): AnalysisRequiredFields[];
    abstract createView(plate: Plate, plateWidget: PlateWidget, currentMappings: Map<string, string>, onMap: (t: string, s: string) => void, onUndo: (t: string) => void): HTMLElement;
    protected abstract _getGroups(resultsDf: DG.DataFrame): { groupColumn: string, groups: string[] };

    /** * An optional helper that subclasses can call within their `createView` if they
     * only need a standard mapping-to-results UI.
     */
    protected _createStandardMappingView(
      plate: Plate,
      currentMappings: Map<string, string>,
      onMap: (t: string, s: string) => void, onUndo: (t: string) => void,
      // ✅ FIX: This signature is now synchronous to match BaseAnalysisView's expectation.
      createResultsView: (p: Plate, m: Map<string, string>) => HTMLElement | null
    ): HTMLElement {
      return new BaseAnalysisView(plate, {
        analysisName: this.friendlyName,
        requiredFields: this.getRequiredFields(),
        createResultsView: createResultsView
      }, currentMappings, onMap, onUndo).getRoot();
    }

    async saveResults(plate: Plate, resultsDf: DG.DataFrame, params: Record<string, any>): Promise<void> {
      if (!plate.id) {
        grok.shell.warning('Please use the main "CREATE" button to save the plate before saving analysis results.');
        return;
      }

      const pi = DG.TaskBarProgressIndicator.create(`Saving ${this.friendlyName} results...`);
      try {
        const {groupColumn, groups} = this._getGroups(resultsDf);
        const runId = await createAnalysisRun(plate.id, this.name, groups);

        // Save parameters
        for (const [name, value] of Object.entries(params)) {
          const prop = this._parameterProperties.get(name);
          if (prop) {
            await saveAnalysisRunParameter({
              runId: runId,
              propertyName: prop.name,
              propertyType: prop.type as DG.TYPE,
              value: value,
            });
          }
        }

        // Save results per group
        for (const row of resultsDf.rows) {
          const groupKey = row.get(groupColumn);
          if (!groupKey) continue;
          const groupCombination = Array.isArray(groupKey) ? groupKey : [groupKey];

          for (const output of this.outputs) {
            const prop = this._outputProperties.get(output.name);
            if (prop && resultsDf.columns.contains(output.name)) {
              const value = row.get(output.name);
              if (value !== null && value !== undefined)
                await saveAnalysisResult({runId, propertyId: prop.id, propertyType: prop.type, value, groupCombination});
            }
          }
        }
        grok.shell.info(`Saved ${this.friendlyName} results for plate ${plate.barcode}.`);
      } catch (e: any) {
        grok.shell.error(`Failed to save ${this.friendlyName} results.`);
        console.error(e);
      } finally {
        pi.close();
      }
    }
}
