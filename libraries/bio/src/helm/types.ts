import type {
  HelmType, PolymerType, MonomerType, WebEditorRGroups, MonomerSetType,
} from '@datagrok-libraries/js-draw-lite/src/types/org';

import type {Point} from '@datagrok-libraries/js-draw-lite/src/Point';
import type {Atom} from '@datagrok-libraries/js-draw-lite/src/Atom';
import type {IBio, IJsAtom} from '@datagrok-libraries/js-draw-lite/src/types/jsdraw2';
import type {Bond} from '@datagrok-libraries/js-draw-lite/src/Bond';
import type {Mol} from '@datagrok-libraries/js-draw-lite/src/Mol';
import type {IMolHandler} from '@datagrok-libraries/js-draw-lite/src/types/mol-handler';
import type {Editor} from '@datagrok-libraries/js-draw-lite/src/JSDraw.Editor';

import type {
  IHelmBio, HelmAtom, HelmBond, HelmMol, HelmEditor, HelmString,
  IOrgHelmWebEditor, IOrgHelmMonomers,
  GetMonomerFunc, GetMonomerResType, MonomersFuncs,
  IHelmEditorOptions, IHelmDrawOptions, IMonomerColors, IWebEditorMonomer,
} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';
import type {HweWindow} from '@datagrok-libraries/helm-web-editor/src/types';
import type {MonomerExplorer} from '@datagrok-libraries/helm-web-editor/helm/MonomerExplorer';
import type {App} from '@datagrok-libraries/helm-web-editor/helm/App';


import type {TabDescType} from '@datagrok-libraries/js-draw-lite/form/Tab';
import type {JSDraw2ModuleType, ScilModuleType} from '@datagrok-libraries/js-draw-lite/src/types';
import type {OrgType} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';
import type {Monomers} from '@datagrok-libraries/helm-web-editor/helm/Monomers';
import type {DojoType, DojoxType} from '@datagrok-libraries/js-draw-lite/src/types/dojo';

export {HelmType, PolymerType, MonomerType, WebEditorRGroups};
export {Point, IBio, Atom, IJsAtom, Bond, Mol, Editor, MonomerExplorer};

export {IHelmBio, HelmAtom, HelmBond, HelmMol, HelmEditor, HelmString};

export {
  IWebEditorMonomer, IMonomerColors, IOrgHelmWebEditor, IOrgHelmMonomers,
  IMolHandler, IHelmEditorOptions, IHelmDrawOptions, GetMonomerFunc, GetMonomerResType,
  App, Monomers,
};

export {DojoType, DojoxType};
export {HweWindow, ScilModuleType, JSDraw2ModuleType, OrgType};


export interface ISeqMonomer {
  position: number;
  biotype: HelmType;
  symbol: string;
}

export interface IHelmWebEditor {
  get editor(): HelmEditor;
  get host(): HTMLDivElement;

  resizeEditor(width: number, height: number): void;
}

export {MonomerSetType, MonomersFuncs, TabDescType};
