import * as d3 from "d3";
import { Branch, TreeData } from './tree-data-builder';
import { ColorMode, OpacityMode, SunburstRendererColor } from './sunburst-renderer-color';
import { HierarchyRectangularNode } from 'd3';

interface Rectangle {
    x0: number;
    y0: number;
    x1: number;
    y1: number;
    data: Branch;
}

export type ClickHandler = (rowIds: string[]) => void;

export class SunburstRenderer {

    private readonly format = d3.format("~r");
    private readonly colorPicker: SunburstRendererColor;

    constructor(colors: string[],
                private readonly clickHandler: ClickHandler) {
        this.colorPicker = new SunburstRendererColor(colors);
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

    public render(htmlElement: HTMLElement, data: TreeData, width: number, height: number, colorMode: ColorMode) {
        const center = Math.min(width, height) / 2;
        const radius = center * 0.9;
        const root = this.partitionLayout(data, radius);

        this.colorPicker.setup(root);
        this.colorPicker.colorMode = colorMode;
        this.colorPicker.opacityMode = OpacityMode.GRADIENT;

        const svg = this.createSvg(root, radius);
        const html = svg.attr("viewBox", `-${center} -${center} ${width} ${height}`).node()!;
        htmlElement.innerHTML = '';
        htmlElement.appendChild(html);
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

    private arcSelection(radius: number) {
        return d3.arc<Rectangle>()
            .startAngle(d => d.x0)
            .endAngle(d => d.x1)
            .padAngle(d => Math.min((d.x1 - d.x0) / 2, 0.005))
            .padRadius(radius / 2)
            .innerRadius(d => d.y0)
            .outerRadius(d => {
                const leavesTotal = d.data.statsOverall?.count || 0;
                const leavesSelected = d.data.statsSelected?.count || 0;
                const ratio = leavesTotal == 0 ? 0 : 1 - leavesSelected / leavesTotal;
                return d.y1 - 1 - (d.y1 - 1 - d.y0) * ratio;
            });
    }

    private createSvg(root: HierarchyRectangularNode<Branch>, radius: number) {
        const elements = root.descendants().filter(d => d.depth);
        const svg = d3.create("svg");

        // Selection (partial) segments
        const segmentFiltered = svg.append("g")
            .selectAll("path")
            .data(elements)
            .join("path")
            .attr("d", this.arcSelection(radius))
            .attr("fill-opacity", d => this.colorPicker.getOpacity(d))
            .attr("fill", d => this.colorPicker.getColor(d));

        // Sunburst segments
        const segment = svg.append("g")
            .selectAll("path")
            .data(elements)
            .join("path")
            .attr("d", this.arc(radius))
            .attr("fill-opacity", d => this.colorPicker.getOpacity(d))
            .attr("fill", d => this.colorPicker.getColor(d));

        segment.on("click", this.onClick)
            .on("mouseover", target => {
                segment.attr("fill-opacity", d => {
                    return SunburstRenderer.sameBranch(target, d) ? 0.8 : this.colorPicker.getOpacity(d);
                });
            })
            .on("mouseleave", target => {
               segment.attr("fill-opacity", d => this.colorPicker.getOpacity(d));
            })
            .append("title")
            .text(this.getTooltipText);

        const fontSize = radius / 20;
        svg.append("g")
            .attr("pointer-events", "none")
            .attr("text-anchor", "middle")
            .attr("font-size", fontSize)
            .attr("font-family", "sans-serif")
            .selectAll("text")
            .data(root.descendants().filter(d => d.depth && (d.y0 + d.y1) / 2 * (d.x1 - d.x0) > fontSize))
            .join("text")
            .attr("transform", function (d) {
                const x = (d.x0 + d.x1) / 2 * 180 / Math.PI;
                const y = (d.y0 + d.y1) / 2;
                return `rotate(${x - 90}) translate(${y},0) rotate(${x < 180 ? 0 : 180})`;
            })
            .attr("dy", "0.35em")
            .attr("fill", "black")
            .attr("fill-opacity", 1)
            // .attr("stroke", "black")
            // .attr("stroke-width", ".5px")
            // .attr("stroke-linecap", "round")
            // .attr("stroke-linejoin", "round")
            // .attr("stroke-opacity", .3)
            // .attr("stroke-alignment", "outer")
            .text((d) => d.data.category);

        return svg
    }

    private getTooltipText = (d: TreeData): string => {
        const selectedCount = d.data.statsSelected?.count || 0;
        const selectedSum = d.data.statsSelected?.sum || 0;
        const overallCount = d.data.statsOverall?.count || 0;
        const overallSum = d.data.statsOverall?.sum || 0;
        const itemPath = this.getCategories(d).join("/");
        const selectedCountStr = d.data.statsSelected?.count !== undefined
            ?  this.format(selectedCount) + ' / '
            : '';

        let tooltipText = `${itemPath}\n` +
            `count: ${selectedCountStr}${this.format(overallCount)}\n`
        if (d.data.statsOverall?.sum === undefined) {
            return tooltipText;
        }

        // If value column is selected
        const selectedSumStr = d.data.statsSelected?.count !== undefined
            ? this.format(selectedSum) + ' / '
            : '';
        const selectedAvgStr = d.data.statsSelected?.count !== undefined
            ? this.format(selectedSum / selectedCount) + ' / '
            : '';
        tooltipText += `sum:    ${selectedSumStr}${this.format(overallSum)} \n` +
            `avg:     ${selectedAvgStr}${this.format(overallSum / overallCount)}`;
        return tooltipText;
    }

    private getCategories = (d: TreeData): string[] => {
        return d.ancestors().map(d => d.data.category).reverse().filter((v, i) => !!i);
    }

    private onClick = (d: TreeData) => {
        this.clickHandler(this.getCategories(d));
    }
}
