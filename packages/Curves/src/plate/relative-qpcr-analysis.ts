// src/plate/relative-qpcr-analysis.ts

/* eslint-disable camelcase */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {Plate, LayerType} from './plate';
import {BaseAnalysisView} from './base-analysis-view';
import {PlateWidget} from './plate-widget';

const REQUIRED_ROLES = {
  TARGET: 'Target Gene',
  REFERENCE: 'Reference Gene',
  CONTROL: 'Control',
  TREATED: 'Treated',
};
function areAllGroupsPresent(plate: Plate): boolean {
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
export class PlateQpcrAnalysis {
  static createAnalysisView(
    plate: Plate,
    currentMappings: Map<string, string>,
    onMap: (target: string, source: string) => void,
    onUndo: (target: string) => void,
    onDataChanged: () => void
  ): HTMLElement {
    console.log('%c--- Rendering qPCR Analysis View ---', 'color: #2083f5; font-weight: bold;');
    const allPrerequisitesMet = areAllGroupsPresent(plate) && currentMappings.has('Ct');

    if (allPrerequisitesMet) {
      const resultsView = this.createQpcrResultsUI(plate, currentMappings, onDataChanged);
      if (resultsView)
        return resultsView;
    }

    const analysisView = new BaseAnalysisView(
      plate, {
        analysisName: 'qPCR Analysis',
        requiredFields: [
          {name: 'Ct', required: true, description: 'Cycle threshold value'},
        ],
        createResultsView: (p, m) => this.createQpcrResultsUI(p, m, onDataChanged),
      },
      currentMappings,
      onMap,
      onUndo
    );

    return analysisView.getRoot();
  }

  private static createQpcrResultsUI(plate: Plate, mappings: Map<string, string>, onDataChanged: () => void): HTMLElement | null {
    if (areAllGroupsPresent(plate)) {
      try {
        const ctColumnName = mappings.get('Ct')!;
        const resultsDf = this.calculateDeltaDeltaCt(plate, ctColumnName);
        const grid = DG.Viewer.grid(resultsDf);
        grid.props.allowEdit = false;
        return grid.root;
      } catch (e: any) {
        console.error('qPCR Calculation Error:', e);
        return ui.divText(`Calculation Error: ${e.message}`, 'error-message');
      }
    }

    const embeddedPlateWidget = new PlateWidget();
    embeddedPlateWidget.plate = plate;
    embeddedPlateWidget.editable = true;
    const createAssignButton = (roleToAssign: string, mutuallyExclusiveRole: string | undefined, onDataChangedCallback: () => void) => {
      return ui.button(`Assign ${roleToAssign}`, () => {
        const selection = plate.data.selection;
        if (selection.trueCount === 0) {
          grok.shell.warning('Please select wells on the plate to assign a role.');
          return;
        }
        const selectedIndexes = selection.getSelectedIndexes();
        let roleCol: DG.Column<boolean>;
        if (plate.data.columns.contains(roleToAssign))
          roleCol = plate.data.col(roleToAssign)!;
        else
          roleCol = plate.data.columns.addNewBool(roleToAssign);

        plate.registerLayer(roleToAssign, LayerType.LAYOUT);
        for (const i of selectedIndexes)
          roleCol.set(i, true);

        if (mutuallyExclusiveRole && plate.data.columns.contains(mutuallyExclusiveRole)) {
          const exclusiveCol = plate.data.col(mutuallyExclusiveRole)!;
          for (const i of selectedIndexes)
            exclusiveCol.set(i, false);
        }
        selection.setAll(false, false);
        embeddedPlateWidget.refresh();
        onDataChangedCallback();
      });
    }; ;


    const statusItems = Object.values(REQUIRED_ROLES).map((role) => {
      const isAssigned = plate.data.columns.contains(role);
      const icon = ui.iconFA(isAssigned ? 'check-circle' : 'times-circle');
      icon.style.color = isAssigned ? 'var(--green-2)' : 'var(--red-3)';
      icon.style.paddingRight = '8px';
      const label = ui.label(role);
      return ui.divH([icon, label], {style: {padding: '2px 0'}});
    });
    const statusPanel = ui.divV([ui.h3('Status'), ...statusItems], 'ui-box');
    statusPanel.style.marginTop = '20px';
    statusPanel.style.padding = '10px';

    const assignmentPanel = ui.divV([
      ui.h3('Assign Roles'),
      ui.p('Select wells on the plate, then assign a role.'),
      ui.divV([
        createAssignButton(REQUIRED_ROLES.TARGET, REQUIRED_ROLES.REFERENCE, onDataChanged),
        createAssignButton(REQUIRED_ROLES.REFERENCE, REQUIRED_ROLES.TARGET, onDataChanged),
      ], {style: {paddingBottom: '10px'}}),
      ui.divV([
        createAssignButton(REQUIRED_ROLES.CONTROL, REQUIRED_ROLES.TREATED, onDataChanged),
        createAssignButton(REQUIRED_ROLES.TREATED, REQUIRED_ROLES.CONTROL, onDataChanged),
      ]),
    ], {style: {maxWidth: '300px', paddingLeft: '20px'}});

    const container = ui.divH([
      embeddedPlateWidget.root,
      assignmentPanel,
    ], {style: {width: '100%', height: '100%'}});

    embeddedPlateWidget.root.style.flexGrow = '1';

    return container;
  }

  private static calculateDeltaDeltaCt(plate: Plate, ctColumnName: string): DG.DataFrame {
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
      console.log(`Found ${values.length} wells for [${role1.name} + ${role2.name}]. Values: [${values.join(', ')}]`);
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
    console.log('Average Ct values:', avgCt);

    const deltaCt_Control = avgCt.targetControl - avgCt.refControl;
    const deltaCt_Treated = avgCt.targetTreated - avgCt.refTreated;
    console.log(`ΔCt (DeltaCt): Control=${deltaCt_Control.toFixed(3)}, Treated=${deltaCt_Treated.toFixed(3)}`);

    const deltaDeltaCt = deltaCt_Treated - deltaCt_Control;
    console.log(`ΔΔCt (DeltaDeltaCt): ${deltaDeltaCt.toFixed(3)}`);

    const foldChange = Math.pow(2, -deltaDeltaCt);
    console.log(`Fold Change (2^-ΔΔCt): ${foldChange.toFixed(3)}`);


    return DG.DataFrame.fromObjects([
      {
        'Category': 'Avg Ct (Control)',
        'Target Gene': avgCt.targetControl.toFixed(2),
        'Reference Gene': avgCt.refControl.toFixed(2),
        'Notes': `N=${targetControlCt.length} / N=${refControlCt.length}`
      },
      {
        'Category': 'Avg Ct (Treated)',
        'Target Gene': avgCt.targetTreated.toFixed(2),
        'Reference Gene': avgCt.refTreated.toFixed(2),
        'Notes': `N=${targetTreatedCt.length} / N=${refTreatedCt.length}`
      },
      {
        'Category': 'ΔCt',
        'Target Gene': deltaCt_Control.toFixed(3),
        'Reference Gene': deltaCt_Treated.toFixed(3),
        'Notes': 'ΔCt = AvgCt(Target) - AvgCt(Reference)'
      },
      {
        'Category': 'ΔΔCt',
        'Target Gene': deltaDeltaCt.toFixed(3),
        'Reference Gene': null,
        'Notes': 'ΔΔCt = ΔCt(Treated) - ΔCt(Control)'
      },
      {
        'Category': 'Fold Change',
        'Target Gene': foldChange.toFixed(3),
        'Reference Gene': null,
        'Notes': 'Calculated as 2^(-ΔΔCt)'
      },
    ])!;
  }
}
