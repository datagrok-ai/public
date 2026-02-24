export const MPO_SCORE_CHANGED_EVENT = 'grok-mpo-score-changed';
export const MPO_PROFILE_CHANGED_EVENT = 'chem-mpo-profile-changed';
export const MPO_PROFILE_DELETED_EVENT = 'chem-mpo-profile-deleted';

export function getNextAvailable(
  base: string,
  existing: Set<string>,
  format: (base: string, num?: number) => string,
): string {
  const first = format(base);
  if (!existing.has(first))
    return first;

  let num = 2;
  while (existing.has(format(base, num)))
    num++;

  return format(base, num);
}

export function generateMpoFileName(profileName: string, existingFileNames: Set<string>): string {
  const base = profileName.trim().replace(/[\s/\\]+/g, '-').replace(/^-|-$/g, '') || 'profile';
  return getNextAvailable(base, existingFileNames, (b, n) => n ? `${b}-${n}.json` : `${b}.json`);
}
