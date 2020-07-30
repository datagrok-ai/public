import * as d3 from "d3";
import { TreeData } from './tree-data-builder';

interface Rectangle {
    x0: number;
    y0: number;
    x1: number;
    y1: number;
}
//
// export class SunburstChart {
//     private htmlElement: HTMLElement;
//     private width: number;
//     private height: number;
//     private radius: number;
//     private data: TreeData;
//
//     private color = d3.scaleOrdinal(d3.quantize(d3.interpolateRainbow, data.children!.length + 1));
//     private format = d3.format(",d");
//     private arc = d3.arc()
//         .startAngle(d => d.x0)
//         .endAngle(d => d.x1)
//         .padAngle(d => Math.min((d.x1 - d.x0) / 2, 0.005))
//         .padRadius(radius / 2)
//         .innerRadius(d => d.y0)
//         .outerRadius(d => d.y1 - 1);
//
//     constructor(private containerElementId: string) {
//         this.htmlElement = document.getElementById(containerElementId)!;
//         this. width = this.htmlElement.offsetWidth * 0.95;
//         this. height = 800
//         this. radius = Math.min(this.width, this.height) / 2;
//     }
//
//     autoBox(this: SVGGraphicsElement) {
//         document.body.appendChild(this);
//         const { x, y, width, height } = this.getBBox();
//         document.body.removeChild(this);
//         return [x, y, width, height];
//     }
//
//     partition(data: TreeItem[]) {
//         return d3.partition()
//         .size([2 * Math.PI, this.radius])
//         (d3.hierarchy(data)
//             .sum(d => d.value)
//             .sort((a, b) => b.value! - a.value!))
//     }
//
//     chart(data: TreeItem[]) {
//         const root = this.partition(data);
//
//         const svg = d3.create("svg");
//
//         svg.append("g")
//             .attr("fill-opacity", 0.6)
//             .selectAll("path")
//             .data(root.descendants().filter(d => d.depth))
//             .join("path")
//             .attr("fill", d => {
//                 while (d.depth > 1) d = d.parent;
//                 return color(d.data.id);
//             })
//             .attr("d", arc)
//             .append("title")
//             .text(d => `${d.ancestors().map(d => d.data.id).reverse().filter((v, i) => !!i).join("/")}\n${format(d.value)}`);
//
//         svg.append("g")
//             .attr("pointer-events", "none")
//             .attr("text-anchor", "middle")
//             .attr("font-size", 10)
//             .attr("font-family", "sans-serif")
//             .selectAll("text")
//             .data(root.descendants().filter(d => d.depth && (d.y0 + d.y1) / 2 * (d.x1 - d.x0) > 10))
//             .join("text")
//             .attr("transform", function (d) {
//                 const x = (d.x0 + d.x1) / 2 * 180 / Math.PI;
//                 const y = (d.y0 + d.y1) / 2;
//                 return `rotate(${x - 90}) translate(${y},0) rotate(${x < 180 ? 0 : 180})`;
//             })
//             .attr("dy", "0.35em")
//             .text(d => d.data.id);
//
//         return svg.attr("viewBox", autoBox).node();
//
//     }
//
//     draw(data: TreeItem[]) {
//         const result = chart(data);
//         htmlElement
//         htmlElement.appendChild(result)
//     }
// }

// https://observablehq.com/@d3/sunburst
export function d3sunburst(containerElementId: string, data: TreeData) {
    console.error(data);

    const htmlElement = document.getElementById(containerElementId)!;
    const width = htmlElement.offsetWidth * 0.95;
    const height = 800
    const radius = Math.min(width, height) / 2;

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

    const color = d3.scaleOrdinal(d3.quantize(d3.interpolateRainbow, data.children!.length + 1))

    const format = d3.format(",d")

    const arc = d3.arc<Rectangle>()
        .startAngle(d => d.x0)
        .endAngle(d => d.x1)
        .padAngle(d => Math.min((d.x1 - d.x0) / 2, 0.005))
        .padRadius(radius / 2)
        .innerRadius(d => d.y0)
        .outerRadius(d => d.y1 - 1)

    const chart = (data: TreeData): SVGSVGElement => {
        const root = partition(data);

        const svg = d3.create("svg");

        svg.append("g")
            .attr("fill-opacity", 0.6)
            .selectAll("path")
            .data<TreeData>(root.descendants().filter((d: any) => d.depth) as any)
            .join("path")
            .attr("fill", d => {
                while (d.depth > 1) d = d.parent!;
                return color(d.data.id);
            })
            .attr("d", arc as any)
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
