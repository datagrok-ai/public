import * as d3 from "d3";
import {TreeData, Branch} from './tree-data-builder';

interface Rectangle {
    x0: number;
    y0: number;
    x1: number;
    y1: number;
}

export class SunburstRenderer {

    private readonly format = d3.format(",d");

    constructor(private readonly colors: string[],
                private readonly clickHandler: (rowIds: number[]) => void) {
    }

    private static shade(d: TreeData) {
        return 0.1 + 0.5 / Math.pow(2, d.depth - 1);
    }

    private static sameBranch(target: TreeData, d: TreeData) {
        do {
            if (d === target) {
                return true;
            }
            target = target.parent!;
        } while (target != null);
        return false;
    }

    public render(htmlElement: HTMLElement, data: TreeData, width: number, height: number) {
        htmlElement.appendChild(this.createSvg(data, width, height));
    }

    private partitionLayout(data: TreeData, radius: number) {
        return d3.partition<Branch>()
            .size([2 * Math.PI, radius])
            (data.sort((a, b) => b.value! - a.value!))
    }

    private arc(radius: number) {
        
        return d3.arc<Rectangle>()
            .startAngle(d => d.x0)
            .endAngle(d => d.x1)
            .padAngle(d => Math.min((d.x1 - d.x0) / 2, 0.005))
            .padRadius(radius / 2)
            .innerRadius(d => d.y0)
            .outerRadius(d => d.y1 - 1);
    }

    private createSvg(data: TreeData, width: number, height: number) {
        function autoBox(this: SVGGraphicsElement): string {
            document.body.appendChild(this);
            const {x, y, width, height} = this.getBBox();
            document.body.removeChild(this);
            return [x, y, width, height].join(',');
        }
        
        const center = Math.min(width, height) / 2;
        const radius = center * 0.9;
        
        const root = this.partitionLayout(data, radius);

        const svg = d3.create("svg");
        
        const segment = svg.append("g")
            .selectAll("path")
            .data(root.descendants().filter(d => d.depth))
            .join("path")
            .attr("fill", this.defaultSegmentFill(root))
            .attr("fill-opacity", SunburstRenderer.shade)
            .attr("d", this.arc(radius));

        segment.on("click", this.onClick)
            .on("mouseover", target => {
                segment.attr("fill-opacity", d => {
                    return SunburstRenderer.sameBranch(target, d) ? 0.8 : SunburstRenderer.shade(d);
                });
            })
            .on("mouseleave", target => {
                segment.attr("fill-opacity", SunburstRenderer.shade);
            })
            .append("title")
            .text(d => `${d.ancestors().map(d => d.data.value).reverse().filter((v, i) => !!i).join("/")}\n${this.format(d.value || 0)}`);

        svg.append("g")
            .attr("pointer-events", "none")
            .attr("text-anchor", "middle")
            .attr("font-size", 10)
            .attr("font-family", "sans-serif")
            .selectAll("text")
            .data(root.descendants().filter(d => d.depth && (d.y0 + d.y1) / 2 * (d.x1 - d.x0) > 10))
            .join("text")
            .attr("transform", function (d) {
                const x = (d.x0 + d.x1) / 2 * 180 / Math.PI;
                const y = (d.y0 + d.y1) / 2;
                return `rotate(${x - 90}) translate(${y},0) rotate(${x < 180 ? 0 : 180})`;
            })
            .attr("dy", "0.35em")
            .text((d) => d.data.value);
        
        return svg.attr("viewBox", `-${center} -${center} ${width} ${height}`).node()!;
    }
    
    private onClick = (d: TreeData) => {
        const rowIds = d.descendants().flatMap(x => x.data.leafIds);
        this.clickHandler(rowIds);
    }
    
    private defaultSegmentFill(root: TreeData) {
        const color = d3.scaleOrdinal(this.colors.slice(0, (root.children?.length || 0) + 1));
        
        return (d: TreeData) => {
            let v: typeof d | null = d;
            while (!!v && v.depth > 1) v = v.parent;
            return color(v?.data.value || 0);
        }
    }
}
