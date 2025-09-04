// base-analysis-view.ts
import * as ui from 'datagrok-api/ui';
import {Plate} from './plate';
import {AnalysisMappingPanel} from '../plates/views/components/analysis-mapping/analysis-mapping-panel';
import {PlateWidget} from './plate-widget';
import './base-analysis-view.css';

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

/**
 * Standardized base class for all analysis views.
 * Handles the mapping panel <-> results view pattern consistently.
 */
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
    this.container = ui.divV([], 'analysis-view-container');
    this.contentHost = ui.divV([], 'analysis-content-host');

    this.setupLayout();
    this.render();
  }

  private setupLayout(): void {
    // Consistent container styling that works for both mapping and results
    this.container.style.cssText = `
      width: 100%;
      height: 100%;
      display: flex;
      flex-direction: column;
      min-height: 400px;
    `;

    this.contentHost.style.cssText = `
      width: 100%;
      flex: 1;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    `;

    this.container.appendChild(this.contentHost);
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
      onMap: this.onMap,
      onUndo: this.onUndo
    });

    const panelRoot = mappingPanel.getRoot();
    // Ensure mapping panel fills the content host properly
    panelRoot.style.cssText = `
      width: 100%;
      flex: 1;
      padding: 16px;
      box-sizing: border-box;
    `;

    this.contentHost.appendChild(panelRoot);
  }

  private renderResults(): void {
    const resultsView = this.config.createResultsView(this.plate, this.currentMappings);

    if (resultsView) {
      // Ensure results view fills the content host properly
      resultsView.style.cssText = `
        width: 100%;
        height: 100%;
        flex: 1;
        display: flex;
        flex-direction: column;
      `;
      this.contentHost.appendChild(resultsView);
    } else {
      const errorMsg = ui.divText(
        `Unable to create ${this.config.analysisName.toLowerCase()} with current mapping.`,
        'warning-message'
      );
      errorMsg.style.cssText = `
        width: 100%;
        flex: 1;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        color: var(--orange-2);
        font-style: italic;
      `;
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
