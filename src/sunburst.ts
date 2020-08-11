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
    clickHandler: (columnName: string, columnIndex: number) => void;
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

    const format = d3.format(",d");

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
            .selectAll("path")
            .data<TreeData>(root.descendants().filter((d: any) => d.depth) as any)
            .join("path")
            .attr("fill", d => {
                while (d.depth > 1) d = d.parent!;
                return color(d.data.id);
            })
            .attr("fill-opacity", d => {
                return 0.3 + 0.4 / Math.pow(2, d.depth - 1);
            })
            .attr("d", arc as any)
            .on("click", d => {
                clickHandler(d.data.id, d.depth - 1);
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

        // var path = vis.selectAll("path")
        //     .data(nodes)
        //     .enter().append("path")
        //     .attr("display", function (d) {
        //         return d.depth ? null : "none";
        //     })
        //     .attr("d", arc)
        //     .style("fill", function (d) {
        //         return getFillValue(d);
        //     }) // 4
        //     .on("mouseover", mouseover); // 5
        //
        // d3.select("#chart-container").on("mouseleave", mouseleave); // 5

        // svg.on("click", function () {
        //     var coords = d3.mouse(this);
        //
        //     console.error({arguments, coords});
        //     // // Normally we go from data to pixels, but here we're doing pixels to data
        //     // var newData = {
        //     //     x: Math.round(xScale.invert(coords[0])),  // Takes the pixel number to convert to number
        //     //     y: Math.round(yScale.invert(coords[1]))
        //     // };
        //     //
        //     // dataset.push(newData);   // Push data to our array
        //     //
        //     // svg.selectAll("circle")  // For new circle, go through the update process
        //     //     .data(dataset)
        //     //     .enter()
        //     //     .append("circle")
        //     //     .attr(circleAttrs)  // Get attributes from circleAttrs var
        //     //     .on("mouseover", handleMouseOver)
        //     //     .on("mouseout", handleMouseOut);
        // })

        return svg.attr("viewBox", autoBox).node()!;
    };

    // function mouseover(d) {
    //     d3.select("#amount")
    //         .text(d.value.toLocaleString('fr-FR', { minimumFractionDigits: 0, style: 'currency', currency: 'EUR' }));
    //
    //     d3.selectAll("path") // On grise tous les segments
    //         .style("opacity", 0.3);
    //
    //     vis.selectAll("path") // Ensuite on met en valeur uniquement ceux qui sont ancêtres de la sélection
    //         .filter(function (node) {
    //             return (sequenceArray.indexOf(node) >= 0);
    //         })
    //         .style("opacity", 1);
    // }
    //
    // function mouseleave(d) {
    //     // On désactive la fonction mouseover le temps de la transition
    //     d3.selectAll("path").on("mouseover", null);
    //
    //     // Transition pour revenir à l'état d'origine et on remet le mouseover
    //     d3.selectAll("path")
    //         .transition()
    //         .duration(1000)
    //         .style("opacity", 1)
    //         .on("end", function () {
    //             d3.select(this).on("mouseover", mouseover);
    //         });
    // }

    const result = chart(data);
    htmlElement.appendChild(result!)
}
