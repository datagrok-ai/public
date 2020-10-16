import * as d3 from "d3";
import {TreeData, Branch} from './tree-data-builder';

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

    constructor(private readonly colors: string[],
                private readonly clickHandler: ClickHandler) {
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

    public render(htmlElement: HTMLElement, data: TreeData, width: number, height: number, colorMode: number) {
        htmlElement.innerHTML = '';
        htmlElement.appendChild(this.createSvg(data, width, height, colorMode));
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

    private createSvg(data: TreeData, width: number, height: number, colorMode: number) {
        const center = Math.min(width, height) / 2;
        const radius = center * 0.9;

        const root = this.partitionLayout(data, radius);
        const elements = root.descendants().filter(d => d.depth);

        const svg = d3.create("svg");

        // Selection (partial) segments
        const segmentFiltered = svg.append("g")
            .selectAll("path")
            .data(elements)
            .join("path")
            .attr("d", this.arcSelection(radius))
        ;
        if (!colorMode) {
            segmentFiltered
                .attr("fill-opacity", SunburstRenderer.opacityMode1)
                .attr("fill", this.colorMode1(root))
        } else {
            segmentFiltered
                .attr("fill-opacity", SunburstRenderer.opacityMode2)
                .attr("fill", this.colorMode2(root))
        }

        // Sunburst segments
        const segment = svg.append("g")
            .selectAll("path")
            .data(elements)
            .join("path")
            .attr("d", this.arc(radius))
        ;
        if (!colorMode) {
            segment
                .attr("fill-opacity", SunburstRenderer.opacityMode1)
                .attr("fill", this.colorMode1(root))
        } else {
            segment
                .attr("fill-opacity", SunburstRenderer.opacityMode2)
                .attr("fill", this.colorMode2(root))
        }

        segment.on("click", this.onClick)
            .on("mouseover", target => {
                segment.attr("fill-opacity", d => {
                    return SunburstRenderer.sameBranch(target, d) ? 0.8 : (!colorMode ? SunburstRenderer.opacityMode1(d) : SunburstRenderer.opacityMode2(d));
                });
            })
            .on("mouseleave", target => {
               !colorMode ? segment.attr("fill-opacity", SunburstRenderer.opacityMode1) : segment.attr("fill-opacity", SunburstRenderer.opacityMode2);
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

        return svg.attr("viewBox", `-${center} -${center} ${width} ${height}`).node()!;
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

    private colorMode1(root: TreeData) {
        const color = d3.scaleOrdinal(this.colors.slice(0, (root.children?.length || 0) + 1));

        return (d: TreeData) => {
            let v: typeof d | null = d;
            while (!!v && v.depth > 1) v = v.parent;
            return color(v?.data?.category || '');
        }
    }

    private colorMode2(root: TreeData) {
        return (d: TreeData) => {
            const nodeColor = d.data?.properties?.color || 'rgb(230, 247, 255)';
            return nodeColor; // d3.interpolate(nodeColor, "white")(.5);
        }
    }

    private colorMode3(root: TreeData) {
        return (d: TreeData) => {
            const nodeColor = d.data?.properties?.color;
            const parentColor = d.parent?.data?.properties?.color;
            if (nodeColor) {
                return parentColor ? d3.interpolateRgb.gamma(2.2)(nodeColor, parentColor)(.3) : nodeColor;
            }
            return 'rgb(230, 247, 255)';
        }
    }

    private static opacityMode2(d: TreeData) {
        return 0.1 + 0.5 / Math.pow(2, d.depth - 1);
    }

    private static opacityMode1(d: TreeData) {
        return 0.3 + 0.5 / Math.pow(2, d.depth - 1);
    }

}
