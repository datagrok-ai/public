// HELM type vocabulary is owned by `@datagrok-libraries/hwe` (the standalone
// HELM editor) and re-exported here as the single source of truth. The string
// values match the canonical Pistoia HELM Web Editor, except `NUCLEOTIDE` which
// uses the corrected spelling `HELM_NUCLEOTIDE` (was `HELM_NUCLETIDE`).
import {HelmTypes, MonomerTypes, PolymerTypes} from '@datagrok-libraries/hwe';

export {HelmTypes, MonomerTypes, PolymerTypes};

// hwe migration: this small enum is not part of the hwe value vocabulary.
// It is defined here locally (replicating the legacy HELM Web Editor values)
// so `@datagrok-libraries/bio` no longer depends on
// `@datagrok-libraries/helm-web-editor`.
export enum MonomerNumberingTypes {
  default = 0,
  continuous = 1,
}
