import {
    PresetStructureRepresentations,
    StructureRepresentationPresetProvider
} from "molstar/lib/mol-plugin-state/builder/structure/representation-preset";
import {TrajectoryHierarchyPresetProvider} from "molstar/lib/mol-plugin-state/builder/structure/hierarchy-preset";
import {StateObjectSelector} from "molstar/lib/mol-state";
import {PluginStateObject} from "molstar/lib/mol-plugin-state/objects";
import {StateObject} from "molstar/lib/mol-state/object";
import {StateTransformer} from "molstar/lib/mol-state/transformer";

type StructureObject = StateObjectSelector<PluginStateObject.Molecule.Structure, StateTransformer<StateObject<any, StateObject.Type<any>>, StateObject<any, StateObject.Type<any>>, any>>

export const RcsbRepresentationPreset: TrajectoryHierarchyPresetProvider = TrajectoryHierarchyPresetProvider({
    id: "rcsb-saguaro-3d",
    display: {
        name: 'Feature View 3D'
    },
    params: () => ({
    }),
    async apply(trajectory, params, plugin) {
        const builder = plugin.builders.structure;
        const model = await builder.createModel(trajectory, {modelIndex: 0});
        const modelProperties = await builder.insertModelProperties(model);
        const structure: StructureObject = await builder.createStructure(modelProperties);
        const structureProperties: StructureObject = await builder.insertStructureProperties(structure);
        const unitcell: StateObjectSelector | undefined = await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true });
        const representation: StructureRepresentationPresetProvider.Result | undefined = await plugin.builders.structure.representation.applyPreset(structureProperties, PresetStructureRepresentations.auto);
        components:
        for (const c of plugin.managers.structure.hierarchy.currentComponentGroups) {
            for (const comp of c) {
                if(comp.cell.obj?.label === "Water") {
                    plugin.managers.structure.component.toggleVisibility(c);
                    break components;
                }
            }
        }
        components:
        for (const c of plugin.managers.structure.hierarchy.currentComponentGroups) {
            for (const comp of c) {
                if(comp.cell.obj?.label === "Polymer") {
                    plugin.managers.structure.component.updateRepresentationsTheme([comp], { color: 'chain-id' });
                    break components;
                }
            }
        }
        return {
            model,
            modelProperties,
            unitcell,
            structure,
            structureProperties,
            representation
        };
    }
});
