// base-analysis-view.ts
import * as ui from 'datagrok-api/ui';
import {Plate} from '../plate';
import {AnalysisMappingPanel} from '../../plates/views/components/analysis-mapping/analysis-mapping-panel';
import './plate-analyses.css'; // Import the new stylesheet

export interface AnalysisFieldMapping {
  name: string;
  required: boolean;
  description?: string;
}

export interface AnalysisViewConfig {
  analysisName: string;
  requiredFields: AnalysisFieldMapping[];
  createResultsView: (plate: Plate, mappings: Map<string, string>) => HTMLElement | null;
}

export class BaseAnalysisView {
  private container: HTMLElement;
  private contentHost: HTMLElement;

  constructor(
    private plate: Plate,
    private config: AnalysisViewConfig,
    private currentMappings: Map<string, string>,
    private onMap: (target: string, source: string) => void,
    private onUndo: (target: string) => void
  ) {
    this.container = ui.divV([], 'assay-plates--analysis-view-container');
    this.contentHost = ui.divV([], 'assay-plates--analysis-content-host');

    this.container.appendChild(this.contentHost);
    this.render();
  }


  private render(): void {
    ui.empty(this.contentHost);

    const requiredFieldNames = this.config.requiredFields
      .filter((f) => f.required)
      .map((f) => f.name);

    const allRequiredMapped = requiredFieldNames.every((field) =>
      this.currentMappings.has(field)
    );

    if (allRequiredMapped)
      this.renderResults();
    else
      this.renderMappingPanel();
  }

  private renderMappingPanel(): void {
    const mappingPanel = new AnalysisMappingPanel({
      analysisName: this.config.analysisName,
      requiredFields: this.config.requiredFields,
      sourceColumns: this.plate.data.columns.names(),
      currentMappings: this.currentMappings,
      onMap: (target: string, source: string) => {
        this.currentMappings.set(target, source);
        this.onMap(target, source);
        this.render();
      },
      onUndo: (target: string) => {
        this.currentMappings.delete(target);
        this.onUndo(target);
        this.render();
      }
    });

    const panelRoot = mappingPanel.getRoot();
    panelRoot.classList.add('assay-plates--analysis-mapping-panel-root');

    this.contentHost.appendChild(panelRoot);
  }

  private renderResults(): void {
    const resultsView = this.config.createResultsView(this.plate, this.currentMappings);

    if (resultsView) {
      resultsView.classList.add('assay-plates--results-view');
      this.contentHost.appendChild(resultsView);
    } else {
      const errorMsg = ui.divText(
        `Unable to create ${this.config.analysisName.toLowerCase()} with current mapping.`,
        'assay-plates--warning-message'
      );
      this.contentHost.appendChild(errorMsg);
    }
  }

  /**
   * Update the view when mappings change
   */
  public updateMappings(newMappings: Map<string, string>): void {
    this.currentMappings = newMappings;
    this.render();
  }

  /**
   * Update the source plate data
   */
  public updatePlate(newPlate: Plate): void {
    this.plate = newPlate;
    this.render();
  }

  /**
   * Get the root container element
   */
  public getRoot(): HTMLElement {
    return this.container;
  }
}
