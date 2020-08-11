import * as d3 from "d3";
import { TreeData } from './tree-data-builder';

interface Rectangle {
    x0: number;
    y0: number;
    x1: number;
    y1: number;
}

export interface D3SunburstParams {
    htmlElement: HTMLElement;
    data: TreeData;
    radius: number;
    clickHandler: (categoryIdsFromTargetToTheRoot: string[], targetNodeDepth: number) => void;
    colors: string[]; // "rgb(123, 45, 6)"
}

// https://observablehq.com/@d3/sunburst
export function d3sunburst(params: D3SunburstParams) {
    const {htmlElement, data, radius, clickHandler, colors} = params;

    function autoBox(this: SVGGraphicsElement): string {
        document.body.appendChild(this);
        const { x, y, width, height } = this.getBBox();
        document.body.removeChild(this);
        return [x, y, width, height].join(',');
    }

    const partition = (data: TreeData) => d3.partition()
        .size([2 * Math.PI, radius])
        (d3.hierarchy(data)
            .sum(d => d.data.value)
            .sort((a, b) => b.value! - a.value!))

    // console.error(d3.quantize(d3.interpolateRainbow, data.children!.length + 1));
    const color = d3.scaleOrdinal(colors.slice(0, data.children!.length + 1));
    const shade = (d: TreeData) => {
        return 0.1 + 0.5 / Math.pow(2, d.depth - 1);
    }

    const format = d3.format(",d");

    const arc = d3.arc<Rectangle>()
        .startAngle(d => d.x0)
        .endAngle(d => d.x1)
        .padAngle(d => Math.min((d.x1 - d.x0) / 2, 0.005))
        .padRadius(radius / 2)
        .innerRadius(d => d.y0)
        .outerRadius(d => d.y1 - 1)

    const sameBranch = (target: TreeData, d: TreeData) => {
        do {
            if (d === target) {
                return true;
            }
            target = target.parent!;
        } while (target != null);
        return false;
    }

    const chart = (data: TreeData): SVGSVGElement => {
        const root = partition(data);

        const svg = d3.create("svg");

        const segment = svg.append("g")
            .selectAll("path")
            .data<TreeData>(root.descendants().filter((d: any) => d.depth) as any)
            .join("path")
            .attr("fill", d => {
                while (d.depth > 1) d = d.parent!;
                return color(d.data.id);
            })
            .attr("fill-opacity", shade)
            .attr("d", arc as any);

            segment.on("click", d => {
                clickHandler(d.ancestors().map(n => n.data.id), d.depth);
            })
            .on("mouseover", target => {
                segment.attr("fill-opacity", d => {
                    return sameBranch(target, d) ? 0.8 : shade(d);
                });
            })
            .on("mouseleave", target => {
                segment.attr("fill-opacity", shade);
            })
            .append("title")
            .text(d => `${d.ancestors().map(d => d.data.id).reverse().filter((v, i) => !!i).join("/")}\n${format(d.data.value)}`);

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
            .text((d: any) => d.data.id);

        return svg.attr("viewBox", autoBox).node()!;
    };

    const result = chart(data);
    htmlElement.appendChild(result!)
}
