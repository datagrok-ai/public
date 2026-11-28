/** A row's group comes from its animal: a row left without one (a day typed into a blinded table, where the group
 *  column is hidden) takes the group its animal has in its other rows. An animal with no group anywhere stays empty,
 *  and so does one whose rows disagree: no guess. */
export function groupsByAnimal<T>(animal: (string | null)[], group: (T | null)[]): (T | null)[] {
  const empty = (g: T | null) => g === null || (g as unknown) === '';
  const known = new Map<string, T | null>();
  animal.forEach((a, i) => {
    if (a === null || empty(group[i]))
      return;
    known.set(a, known.has(a) && known.get(a) !== group[i] ? null : group[i]);
  });
  return group.map((g, i) => !empty(g) || animal[i] === null ? g : known.get(animal[i]!) ?? g);
}
