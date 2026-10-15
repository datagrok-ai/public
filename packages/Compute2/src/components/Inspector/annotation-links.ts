import {AnnotationLinkKind} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';

/** Missing-value checks are never shown; other links generated from annotations only with the toggle on. */
export const isLinkShown = (annotation: AnnotationLinkKind | undefined, showAnnotationLinks: boolean) =>
  annotation !== 'required' && (showAnnotationLinks || annotation == null);

