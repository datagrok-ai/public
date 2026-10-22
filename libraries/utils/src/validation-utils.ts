export const isLinkedRestriction = (restriction?: string) =>
  restriction === 'disabled' || restriction === 'restricted';

export const hasValidationIssue = (
  validation?: {errors?: unknown[], warnings?: unknown[]},
  consistency?: {inconsistent?: boolean},
) => !!validation?.errors?.length || !!validation?.warnings?.length || !!consistency?.inconsistent;
