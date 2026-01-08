import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {DesirabilityProfile} from './mpo';
import {MpoDesirabilityLineEditor} from './mpo-line-editor';
import {Subject} from 'rxjs';

import '../../css/styles.css';

export class MpoProfileEditor {
  root = ui.div([]);
  dataFrame?: DG.DataFrame;
  onChanged = new Subject();
  profile?: DesirabilityProfile;

  constructor(dataFrame?: DG.DataFrame) {
    this.dataFrame = dataFrame;
    this.setProfile();
  }

  getProfile(): DesirabilityProfile | undefined {
    return this.profile;
  }

  setProfile(profile?: DesirabilityProfile) {
    this.profile = profile;
    ui.empty(this.root);
    if (!profile) {
      this.root.append(ui.divText('No profile specified.'));
      return;
    }

    // Create header row for the table
    const header = ui.divH([
      ui.divText('Property', 'statistics-mpo-header-property'),
      ui.divText('Weight', 'statistics-mpo-header-weight'),
      ui.divText('Desirability', 'statistics-mpo-header-desirability'), // Let editor take space
    ], 'statistics-mpo-header');

    const propertyRows = Object.entries(profile.properties).map(([propertyName, prop]) => {
      const lineEditor = new MpoDesirabilityLineEditor(prop, 300, 80);
      lineEditor.onChanged.subscribe((_) => this.onChanged.next());

      // Input for weight - updates the *copy* of the template data
      const weightInput = ui.input.float('', { // No label needed here
        value: prop.weight, min: 0, max: 1,
        onValueChanged: (newValue) => { // Changed parameter name for clarity
          if (profile && profile.properties[propertyName]) {
            // Update the weight in the temporary template object
            let clampedWeight = newValue ?? 0;
            clampedWeight = Math.max(0, Math.min(1, clampedWeight)); // Clamp 0-1
            profile.properties[propertyName].weight = clampedWeight;
          }
        },
      });
      weightInput.root.classList.add('statistics-mpo-weight-input');

      const matchedColumnName = this.dataFrame ?
        this.dataFrame.columns.names().find((name) => name.toLowerCase() == propertyName.toLowerCase()) :
        null;

      const columnInput = ui.input.choice('', {
        nullable: true,
        items: this.dataFrame?.columns?.names() ?? [''],
        value: matchedColumnName ?? ''});

      const rowDiv = ui.divH([
        ui.divV([
          ui.divText(propertyName, 'statistics-mpo-property-name'),
          this.dataFrame ? columnInput.root : null,
        ]),
        weightInput.root,
        lineEditor.root, // Add the Konva container div
      ]);
      rowDiv.classList.add('statistics-mpo-row');

      return rowDiv;
    }).filter((el) => el !== null); // Filter out skipped properties

    if (propertyRows.length > 0) {
      this.root.append(header);
      this.root.append(ui.divV(propertyRows as HTMLElement[])); // Cast needed after filter
    } else
      this.root.append(ui.divText('No matching properties found in the table for this template.'));
  }
}
