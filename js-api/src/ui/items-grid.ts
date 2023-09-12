import {Property} from "../entities";
import {div, empty, input, span} from "../../ui";

export class ItemsGrid {
  root: HTMLDivElement = div([], 'd4-items-grid');
  items: any[] = [];
  props: Property[] = [];

  constructor(items: any[], props: Property[]) {
    this.items = items;
    this.props = props;
    this.render();
  }

  render(): void {
    this.root.style.gridTemplateColumns = `repeat(${this.props.length}, 1fr)`;
    empty(this.root);

    for (const p of this.props)
      this.root.appendChild(span([p.name], 'd4-items-grid-column-header'));

    for (const item of this.items)
      for (const p of this.props) {
        this.root.appendChild(input.forProperty(p, item).input);
      }
  }
}