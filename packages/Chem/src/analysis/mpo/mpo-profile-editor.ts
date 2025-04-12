import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {DesirabilityProfile} from "./mpo";
import {MpoDesirabilityLineEditor} from "./mpo-line-editor";
import {Observable, Subject} from "rxjs";


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
      ui.divText('Property', { style: { fontWeight: 'bold', width: '150px' } }),
      ui.divText('Weight', { style: { fontWeight: 'bold', width: '60px' } }),
      ui.divText('Desirability', { style: { fontWeight: 'bold', flexGrow: '1' } }) // Let editor take space
    ], { style: {
            marginTop: '10px',
            paddingBottom: '5px',
            borderBottom: '1px solid #ccc'
          }
    });

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
        }
      });

      weightInput.root.style.width = '60px';
      weightInput.root.style.marginTop = '21px';

      const matchedColumnName = this.dataFrame
        ? this.dataFrame.columns.names().find(name => name.toLowerCase() == propertyName.toLowerCase())
        : null;

      const columnInput = ui.input.choice('', {
        nullable: true,
        items: this.dataFrame?.columns?.names() ?? [''],
        value: matchedColumnName ?? ''});

      const rowDiv = ui.divH([
        ui.divV([
          ui.divText(propertyName, { style: { width: '150px', paddingTop: '5px', marginLeft: '4px'} }),
          this.dataFrame ? columnInput.root : null,
        ]),
        weightInput.root,
        lineEditor.root // Add the Konva container div
      ]);
      rowDiv.style.alignItems = 'center'; // Vertically align items in the row
      rowDiv.style.marginBottom = '5px'; // Space between rows
      rowDiv.style.minHeight = '70px'; // Ensure consistent row height even if editor takes time

      return rowDiv;
    }).filter(el => el !== null); // Filter out skipped properties

    if (propertyRows.length > 0) {
      this.root.append(header);
      this.root.append(ui.divV(propertyRows as HTMLElement[])); // Cast needed after filter
    }
    else
      this.root.append(ui.divText('No matching properties found in the table for this template.'));
  }
}
