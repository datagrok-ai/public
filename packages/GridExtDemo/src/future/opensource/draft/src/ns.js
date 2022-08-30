export const CHEM = function() {throw new Error("The function cannot be called. Use the static properties.")};
import {MolService} from "./chem/MolService";
CHEM.MolService = MolService;

CHEM.RDKit = function() {throw new Error("The function cannot be called. Use the static properties.")};
import {RDKitEngine} from "./chem/rdkit/RDKitEngine";
CHEM.RDKit.RDKitEngine = RDKitEngine;

CHEM.OCL = function() {throw new Error("The function cannot be called. Use the static properties.")};
import {OCLEngine} from "./chem/ocl/OCLEngine";
CHEM.OCL.OCLEngine = OCLEngine;

export const ENTITY = function() {throw new Error("The function cannot be called. Use the static properties.")};
import {SemEntity} from "./entity/SemEntity";
import {SemType} from "./entity/SemType";
import {SemSorter} from "./entity/SemSorter";
ENTITY.SemEntity = SemEntity;
ENTITY.SemType = SemType;
ENTITY.SemSorter = SemSorter;

ENTITY.UI = function() {};
import {SemEntityRenderer} from "./entity/ui/SemEntityRenderer";
import {FlowTextRenderer} from "./entity/ui/FlowTextRenderer";
import {KeyValueRenderer} from "./entity/ui/KeyValueRenderer";
import  {AbstractVertLayoutTextRenderer} from "./entity/ui/AbstractVertLayoutTextRenderer";
ENTITY.UI.SemEntityRenderer = SemEntityRenderer;
ENTITY.UI.FlowTextRenderer = FlowTextRenderer;
ENTITY.UI.KeyValueRenderer = KeyValueRenderer;
ENTITY.UI.AbstractVertLayoutTextRenderer = AbstractVertLayoutTextRenderer;