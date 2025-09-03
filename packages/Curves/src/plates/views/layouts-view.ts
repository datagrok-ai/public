import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateWidget} from '../../plate/plate-widget';
import {Plate} from '../../plate/plate';
import {LayoutPattern, LayoutStateManager} from './shared/layout-state-manager';

/** Creates an interactive control for applying a layout pattern with a specific role. */
function createPatternControl(
  pattern: LayoutPattern,
  onApply: (pattern: LayoutPattern, role: string) => void
): HTMLElement {
  const roleInput = ui.input.string('', {value: 'Sample'});
  const applyBtn = ui.button('Apply', () => onApply(pattern, roleInput.value), 'Apply this pattern');

  // Use standard Datagrok classes for alignment and spacing
  const control = ui.divH([
    ui.divText(pattern.name, 'pattern-name'),
    roleInput.root,
    applyBtn,
  ], 'd4-flex-row d4-justify-space-between d4-align-center');

  control.style.padding = '4px 0'; // Add some vertical padding
  roleInput.root.style.marginLeft = 'auto'; // Push input to the right
  roleInput.root.style.marginRight = '8px';
  roleInput.root.style.maxWidth = '150px';

  return control;
}

/** Main view for the interactive plate layout designer. */
export function layoutsView(): DG.View {
  const view = DG.View.create();
  view.name = 'Layouts';
  view.box = true;

  // --- Define Available Patterns ---
  const patterns: LayoutPattern[] = [
    {
      name: 'Checkerboard',
      apply: (plate: Plate, role: string) => {
        let roleCol = plate.data.col('Role');
        if (roleCol === null)
          roleCol = plate.data.columns.addNewString('Role');

        for (let r = 0; r < plate.rows; r++) {
          for (let c = 0; c < plate.cols; c++) {
            if ((r + c) % 2 === 0)
              roleCol.set(plate._idx(r, c), role);
          }
        }
      },
    },
    {
      name: 'Vertical Stripes',
      apply: (plate: Plate, role: string) => {
        let roleCol = plate.data.col('Role');
        if (roleCol === null)
          roleCol = plate.data.columns.addNewString('Role');

        for (let c = 0; c < plate.cols; c += 2) {
          for (let r = 0; r < plate.rows; r++)
            roleCol.set(plate._idx(r, c), role);
        }
      },
    },
    {
      name: 'Control Columns',
      apply: (plate: Plate, role: string) => {
        let roleCol = plate.data.col('Role');
        if (roleCol === null)
          roleCol = plate.data.columns.addNewString('Role');

        // Apply to first and last column
        for (let r = 0; r < plate.rows; r++) {
          roleCol.set(plate._idx(r, 0), role);
          roleCol.set(plate._idx(r, plate.cols - 1), role);
        }
      },
    },
  ];

  // --- Initialize State and Components ---
  const initialPlate = new Plate(8, 12);
  const stateManager = new LayoutStateManager(initialPlate);
  const plateWidget = new PlateWidget();
  plateWidget.plate = stateManager.plate;

  // --- UI Panels ---
  const patternCatalogue = ui.divV(
    patterns.map((p) => createPatternControl(p, stateManager.applyLayout.bind(stateManager))),
    'ui-form' // Use a standard form class for better styling
  );

  const historyPanel = ui.divV([], 'ui-form');

  // --- Main Layout ---
  // Use a splitter for a resizable layout
  const mainLayout = ui.splitH([
    // Left side contains the plate and history
    ui.splitV([
      plateWidget.root,
      ui.divV([ui.h2('History'), historyPanel], {style: {padding: '10px'}}),
    ]),
    // Right side contains the patterns
    ui.divV([ui.h2('Patterns'), patternCatalogue], {style: {padding: '10px'}}),
  ]);
  // CORRECTED: Cast to 'any' to access the special 'ratio' property.
  (mainLayout as any).ratio = 0.65; // Give the plate widget more initial space

  view.root.appendChild(mainLayout);

  // --- State Change Handling ---
  const renderHistory = () => {
    ui.empty(historyPanel);
    stateManager.history.forEach((action) => {
      const deleteIcon = ui.iconFA('times', () => stateManager.undoAction(action.id), 'Undo this action');
      const historyItem = ui.divH([
        ui.divText(`${action.patternName} (${action.role})`),
        deleteIcon,
      ], 'd4-flex-row d4-justify-space-between d4-align-center');
      historyItem.style.padding = '4px 0';
      historyPanel.appendChild(historyItem);
    });
  };

  view.subs.push(stateManager.onStateChange.subscribe(() => {
    plateWidget.plate = stateManager.plate;
    plateWidget.refresh();
    renderHistory();
  }));

  renderHistory(); // Initial render

  return view;
}


