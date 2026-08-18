/** How to show an object of type `T` on u2 surfaces — the platform-free twin of the
 * platform's ObjectHandler render surface (`u2/dg`'s `handlerRenderer` adapts one to the
 * other). Collection controls take a renderer wherever they take per-item callbacks; an
 * explicit callback always wins over the renderer's member. */
export interface ObjectRenderer<T> {
  caption(x: T): string;
  icon?(x: T): HTMLElement;
  /** A row in a popup or list. */
  listItem?(x: T): HTMLElement;
  /** Inline content — chips, breadcrumbs, selected values. */
  markup?(x: T): HTMLElement;
  tooltip?(x: T): HTMLElement | string;
  card?(x: T): HTMLElement;
}
