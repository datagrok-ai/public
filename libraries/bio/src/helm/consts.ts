// HELM type vocabulary is owned by `@datagrok-libraries/hwe` (the standalone
// HELM editor) and re-exported here as the single source of truth. The string
// values match the canonical Pistoia HELM Web Editor, except `NUCLEOTIDE` which
// uses the corrected spelling `HELM_NUCLEOTIDE` (was `HELM_NUCLETIDE`).
import {HelmTypes, MonomerTypes, PolymerTypes} from '@datagrok-libraries/hwe';
import {HelmTabKeys, MonomerNumberingTypes} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';

export {HelmTypes, MonomerTypes, PolymerTypes, MonomerNumberingTypes, HelmTabKeys};
