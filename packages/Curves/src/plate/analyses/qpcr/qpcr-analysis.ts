/* eslint-disable camelcase */
/* eslint-disable max-len */

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {Plate, LayerType} from '../../plate';
import {AnalysisBase, IAnalysisProperty} from '../base-analysis';
import {AnalysisRequiredFields} from '../../../plates/views/components/analysis-mapping/analysis-mapping-panel';
import {PlateWidget} from '../../plate-widget';
import {createAnalysisRun, saveAnalysisResult, saveAnalysisRunParameter} from '../../../plates/plates-crud';
import './../plate-analyses.css';

const REQUIRED_ROLES = {
  TARGET: 'Target Gene',
  REFERENCE: 'Reference Gene',
  CONTROL: 'Control',
  TREATED: 'Treated',
};

export class QpcrAnalysis extends AnalysisBase {
  readonly name: string = 'qPCR';
  readonly friendlyName: string = 'qPCR Analysis';

  parameters: IAnalysisProperty[] = [];
  outputs: IAnalysisProperty[] = [
    {name: 'Delta Delta Ct', type: DG.TYPE.FLOAT, category: 'Output'},
    {name: 'Fold Change', type: DG.TYPE.FLOAT, category: 'Output'},
  ];

  getRequiredFields(): AnalysisRequiredFields[] {
    return [{name: 'Ct', required: true, description: 'Cycle threshold value'}];
  }

  createView(
    plate: Plate,
    plateWidget: PlateWidget,
    currentMappings: Map<string, string>,
    onMap: (t: string, s: string) => void,
    onUndo: (t: string) => void,
    onRerender?: () => void

  ): HTMLElement {
    return this._createStandardMappingView(plate, currentMappings, onMap, onUndo, (p, m) => {
      if (this._areAllRolesPresent(p))
        return this._createResultsGrid(p, m);
      else
        return this._createRoleAssignmentView(p, onRerender);
    });
  }

  override formatResultsForGrid(rawResults: DG.DataFrame): DG.DataFrame {
    if (rawResults.rowCount === 0)
      return rawResults;

    const propsCol = rawResults.col('properties');
    if (!propsCol) {
      console.error('qPCR: CRITICAL - `properties` column not found!');
      return DG.DataFrame.create(0);
    }

    const finalRows: any[] = [];
    for (let i = 0; i < rawResults.rowCount; i++) {
      const row = rawResults.row(i);
      const propsJson = propsCol.get(i);
      if (!propsJson) continue;

      try {
        const properties = JSON.parse(propsJson);
        finalRows.push({
          'run_id': row.get('run_id'),
          'plate_id': row.get('plate_id'),
          'barcode': row.get('barcode'),
          'Delta Delta Ct': properties['Delta Delta Ct'],
          'Fold Change': properties['Fold Change'],
        });
      } catch (e) {
        console.error(`Failed to parse qPCR properties JSON for row ${i}:`, propsJson, e);
      }
    }

    if (finalRows.length === 0)
      return DG.DataFrame.create(0);

    const resultDf = DG.DataFrame.fromObjects(finalRows)!;

    const ddCtCol = resultDf.col('Delta Delta Ct');
    if (ddCtCol) ddCtCol.meta.format = '0.000';

    const foldChangeCol = resultDf.col('Fold Change');
    if (foldChangeCol) foldChangeCol.meta.format = '0.000';

    return resultDf;
  }

  private _createResultsGrid(plate: Plate, mappings: Map<string, string>): HTMLElement {
    try {
      const ctColumnName = mappings.get('Ct')!;
      const resultsDf = this._calculateDeltaDeltaCt(plate, ctColumnName);
      const grid = DG.Viewer.grid(resultsDf);
      grid.props.allowEdit = false;

      const saveButton = ui.button('SAVE RESULTS', async () => {
        ui.setUpdateIndicator(saveButton, true);
        try {
          await this.saveResults(plate, resultsDf, {}, mappings);
        } catch (e) {
        } finally {
          ui.setUpdateIndicator(saveButton, false);
        }
      });

      const container = ui.divV([grid.root, ui.div([saveButton], 'ui-box')], 'assay-plates--analysis-grid-container');
      return container;
    } catch (e: any) {
      console.error('qPCR Calculation Error:', e);
      return ui.divText(`Calculation Error: ${e.message}`, 'assay-plates--error-message');
    }
  }

  private _createRoleAssignmentView(plate: Plate, onRerender?: () => void): HTMLElement {
    const embeddedPlateWidget = new PlateWidget();
    embeddedPlateWidget.plate = plate;
    embeddedPlateWidget.editable = true;

    const createAssignButton = (roleToAssign: string, mutuallyExclusiveRole: string | undefined) => {
      return ui.button(`Assign ${roleToAssign}`, () => {
        const selection = embeddedPlateWidget.plate.data.selection;
        if (selection.trueCount === 0) {
          grok.shell.warning('Please select wells on the plate to assign a role.');
          return;
        }

        const selectedIndexes = selection.getSelectedIndexes();
        const roleCol = plate.data.columns.contains(roleToAssign) ? plate.data.col(roleToAssign)! : plate.data.columns.addNewBool(roleToAssign);
        plate.registerLayer(roleToAssign, LayerType.LAYOUT);
        for (const i of selectedIndexes) roleCol.set(i, true);
        if (mutuallyExclusiveRole && plate.data.columns.contains(mutuallyExclusiveRole)) {
          const exclusiveCol = plate.data.col(mutuallyExclusiveRole)!;
          for (const i of selectedIndexes) exclusiveCol.set(i, false);
        }
        selection.setAll(false, false);

        embeddedPlateWidget.refresh();

        if (onRerender)
          onRerender();
      });
    };

    const statusPanel = this._createStatusPanel(plate);
    if (onRerender) {
      plate.onDataChanged.subscribe(() => {
        const newPanel = this._createStatusPanel(plate);
        ui.empty(statusPanel).append(...Array.from(newPanel.childNodes));
      });
    }

    const assignmentPanel = ui.divV([
      ui.h3('Assign Roles'),
      ui.p('Select wells on the plate below, then assign a role.'),
      ui.buttonsInput([
        createAssignButton(REQUIRED_ROLES.TARGET, REQUIRED_ROLES.REFERENCE),
        createAssignButton(REQUIRED_ROLES.REFERENCE, REQUIRED_ROLES.TARGET),
      ]),
      ui.buttonsInput([
        createAssignButton(REQUIRED_ROLES.CONTROL, REQUIRED_ROLES.TREATED),
        createAssignButton(REQUIRED_ROLES.TREATED, REQUIRED_ROLES.CONTROL),
      ]),
      statusPanel,
    ], 'ui-box assay-plates--qpcr-assignment-panel');

    const container = ui.divH([
      embeddedPlateWidget.root,
      assignmentPanel,
    ], 'assay-plates--qpcr-role-view-container');

    return container;
  }

  private _createStatusPanel(plate: Plate): HTMLElement {
    const statusItems = Object.values(REQUIRED_ROLES).map((role) => {
      const isAssigned = plate.data.columns.contains(role);
      const icon = ui.iconFA(isAssigned ? 'check-circle' : 'times-circle', null, isAssigned ? 'Assigned' : 'Not Assigned');

      icon.classList.add(
        'assay-plates--qpcr-status-icon',
        isAssigned ? 'assay-plates--qpcr-status-icon--assigned' : 'assay-plates--qpcr-status-icon--unassigned'
      );

      return ui.divH([icon, ui.label(role)], 'assay-plates--qpcr-status-item');
    });
    return ui.divV([ui.h3('Status'), ...statusItems], 'assay-plates--qpcr-status-panel');
  }

  async saveResults(
    plate: Plate,
    resultsDf: DG.DataFrame,
    params: Record<string, any>,
    currentMappings: Map<string, string>
  ): Promise<void> {
    if (!plate.id) {
      grok.shell.error('Plate must be saved before saving analysis results.');
      return;
    }
    const pi = DG.TaskBarProgressIndicator.create('Saving qPCR results...');
    try {
      const runId = await createAnalysisRun(plate.id, this.name, Object.values(REQUIRED_ROLES));

      for (const [requiredField, actualColumn] of currentMappings.entries()) {
        const mappingPropName = `Mapping: ${requiredField}`;
        await saveAnalysisRunParameter({
          runId: runId,
          propertyName: mappingPropName,
          propertyType: DG.TYPE.STRING,
          value: actualColumn,
        });
      }

      const groupCombination: string[] = [];
      const ddCt = parseFloat(resultsDf.cell(3, 'Target Gene').value);
      const foldChange = parseFloat(resultsDf.cell(4, 'Target Gene').value);

      if (!isNaN(ddCt))
        await this.saveAnalysisProperty(runId, 'Delta Delta Ct', ddCt, groupCombination);
      if (!isNaN(foldChange))
        await this.saveAnalysisProperty(runId, 'Fold Change', foldChange, groupCombination);
      grok.shell.info('Saved qPCR analysis results.');
    } catch (e: any) {
      grok.shell.error(`Failed to save qPCR results: ${e.message}`);
      console.error(e);
    } finally {
      pi.close();
    }
  }

  private async saveAnalysisProperty(runId: number, propName: string, value: any, group: string[]): Promise<void> {
    const prop = this._outputProperties.get(propName);
    if (prop && value !== null && value !== undefined)
      await saveAnalysisResult({runId, propertyId: prop.id, propertyName: prop.name, propertyType: prop.type, value, groupCombination: group});
  }


  private _areAllRolesPresent(plate: Plate): boolean {
    const df = plate.data;
    const roleNames = Object.values(REQUIRED_ROLES);

    if (!roleNames.every((role) => df.columns.contains(role)))
      return false;

    const hasCombination = (role1Name: string, role2Name: string): boolean => {
      const role1Col = df.col(role1Name)!;
      const role2Col = df.col(role2Name)!;
      for (let i = 0; i < df.rowCount; i++) {
        if (role1Col.get(i) && role2Col.get(i))
          return true;
      }
      return false;
    };

    return hasCombination(REQUIRED_ROLES.TARGET, REQUIRED_ROLES.CONTROL) &&
           hasCombination(REQUIRED_ROLES.TARGET, REQUIRED_ROLES.TREATED) &&
           hasCombination(REQUIRED_ROLES.REFERENCE, REQUIRED_ROLES.CONTROL) &&
           hasCombination(REQUIRED_ROLES.REFERENCE, REQUIRED_ROLES.TREATED);
  }


  private _calculateDeltaDeltaCt(plate: Plate, ctColumnName: string): DG.DataFrame {
    const df = plate.data;
    const ctCol = df.col(ctColumnName)!;
    const targetCol = df.col(REQUIRED_ROLES.TARGET)!;
    const refCol = df.col(REQUIRED_ROLES.REFERENCE)!;
    const controlCol = df.col(REQUIRED_ROLES.CONTROL)!;
    const treatedCol = df.col(REQUIRED_ROLES.TREATED)!;

    const getCtValues = (role1: DG.Column<boolean>, role2: DG.Column<boolean>): number[] => {
      const values: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        if (role1.get(i) && role2.get(i)) {
          const ct = ctCol.get(i);
          if (ct !== null && Number.isFinite(ct))
            values.push(ct);
        }
      }
      if (values.length === 0)
        throw new Error(`No wells found for combination: ${role1.name} + ${role2.name}`);
      return values;
    };
    const average = (arr: number[]) => arr.reduce((a, b) => a + b, 0) / arr.length;
    const targetControlCt = getCtValues(targetCol, controlCol);
    const targetTreatedCt = getCtValues(targetCol, treatedCol);
    const refControlCt = getCtValues(refCol, controlCol);
    const refTreatedCt = getCtValues(refCol, treatedCol);
    const avgCt = {
      targetControl: average(targetControlCt),
      targetTreated: average(targetTreatedCt),
      refControl: average(refControlCt),
      refTreated: average(refTreatedCt),
    };
    const deltaCt_Control = avgCt.targetControl - avgCt.refControl;
    const deltaCt_Treated = avgCt.targetTreated - avgCt.refTreated;
    const deltaDeltaCt = deltaCt_Treated - deltaCt_Control;
    const foldChange = Math.pow(2, -deltaDeltaCt);

    return DG.DataFrame.fromObjects([
      {'Category': 'Avg Ct (Control)', 'Target Gene': avgCt.targetControl.toFixed(2), 'Reference Gene': avgCt.refControl.toFixed(2), 'Notes': `N=${targetControlCt.length} / N=${refControlCt.length}`},
      {'Category': 'Avg Ct (Treated)', 'Target Gene': avgCt.targetTreated.toFixed(2), 'Reference Gene': avgCt.refTreated.toFixed(2), 'Notes': `N=${targetTreatedCt.length} / N=${refTreatedCt.length}`},
      {'Category': 'ΔCt', 'Target Gene': deltaCt_Control.toFixed(3), 'Reference Gene': deltaCt_Treated.toFixed(3), 'Notes': 'ΔCt = AvgCt(Target) - AvgCt(Reference)'},
      {'Category': 'ΔΔCt', 'Target Gene': deltaDeltaCt.toFixed(3), 'Reference Gene': '', 'Notes': 'ΔΔCt = ΔCt(Treated) - ΔCt(Control)'},
      {'Category': 'Fold Change', 'Target Gene': foldChange.toFixed(3), 'Reference Gene': '', 'Notes': 'Calculated as 2^(-ΔΔCt)'},
    ])!;
  }

  protected override _getGroups(resultsDf: DG.DataFrame): { groupColumn: string; groups: string[]; } {
    // qPCR analysis results are aggregated for the whole plate, so there are no specific groups.
    return {groupColumn: '', groups: []};
  }
}
