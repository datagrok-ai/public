/** The ligands a viewer holds: the selected rows', then the current row's and the mouse-over row's,
 * which may repeat a row. */
export function ligandMapItems<T>(map: {selected: T[], current: T | null, hovered: T | null}): T[] {
  return [...map.selected, ...(map.current ? [map.current] : []), ...(map.hovered ? [map.hovered] : [])];
}
