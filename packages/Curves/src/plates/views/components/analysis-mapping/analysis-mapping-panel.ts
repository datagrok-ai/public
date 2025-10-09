import * as ui from 'datagrok-api/ui';
import {renderMappingEditor, TargetProperty} from '../mapping-editor/mapping-editor';

export interface AnalysisRequiredFields {
  name: string;
  required: boolean;
  description?: string;
}

export interface AnalysisMappingConfig {
  analysisName: string;
  requiredFields: AnalysisRequiredFields[];
  sourceColumns: string[];
  currentMappings: Map<string, string>;

  onMap: (target: string, source: string) => void;
  onUndo: (target: string) => void;
}

export class AnalysisMappingPanel {
  private root: HTMLElement;
  private mappingHost: HTMLElement;

  constructor(private config: AnalysisMappingConfig) {
    this.root = ui.divV([], 'analysis-mapping-panel');
    this.mappingHost = ui.divV([]);
    this.buildPanel();
  }

  private buildPanel(): void {
    const header = ui.divH([
      ui.h3(`${this.config.analysisName} Required Fields`),
    ], 'analysis-mapping-header');
    header.style.marginBottom = '12px';

    const infoText = ui.divText(
      `Map your data columns to the required fields for ${this.config.analysisName} analysis.`,
      'analysis-mapping-info'
    );
    infoText.style.fontSize = '12px';
    infoText.style.color = 'var(--grey-3)';
    infoText.style.marginBottom = '16px';

    this.root.appendChild(header);
    this.root.appendChild(infoText);
    this.root.appendChild(this.mappingHost);

    this.renderMapping();
  }

  private renderMapping(): void {
    const targetProperties: TargetProperty[] = this.config.requiredFields.map((field) => ({
      name: field.name,
      required: field.required,
    }));

    renderMappingEditor(this.mappingHost, {
      targetProperties,
      sourceColumns: this.config.sourceColumns,
      mappings: this.config.currentMappings,
      onMap: (target, source) => this.config.onMap(target, source), // Make sure order matches
      onUndo: (target) => this.config.onUndo(target),
    });
  }
  public updateConfig(newConfig: Partial<AnalysisMappingConfig>): void {
    Object.assign(this.config, newConfig);
    this.renderMapping();
  }

  public getRoot(): HTMLElement {
    return this.root;
  }

  // Helper method to check if all required fields are mapped
  public hasAllRequiredFieldsMapped(): boolean {
    const requiredFields = this.config.requiredFields.filter((f) => f.required);
    return requiredFields.every((field) =>
      this.config.currentMappings.has(field.name)
    );
  }

  // Get validation status for UI feedback
  public getValidationStatus(): { isValid: boolean; missingFields: string[] } {
    const requiredFields = this.config.requiredFields.filter((f) => f.required);
    const missingFields = requiredFields
      .filter((field) => !this.config.currentMappings.has(field.name))
      .map((field) => field.name);

    return {
      isValid: missingFields.length === 0,
      missingFields
    };
  }
}
