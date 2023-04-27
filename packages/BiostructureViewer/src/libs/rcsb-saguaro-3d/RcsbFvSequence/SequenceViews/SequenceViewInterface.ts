import {SaguaroPluginModelMapType} from "../../RcsbFvStructure/StructurePlugins/SaguaroPluginInterface";

export interface SequenceViewInterface {
    structureSelectionCallback(): void;
    structureHoverCallback(): void;
    representationChangeCallback(): void;
    modelChangeCallback(modelMap:SaguaroPluginModelMapType): void;
    updateDimensions(): void;
    additionalContent(): JSX.Element | null;
}
