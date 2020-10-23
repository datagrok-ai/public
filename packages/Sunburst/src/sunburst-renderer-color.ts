import { TreeData } from './tree-data-builder';
import * as d3 from 'd3';
import { ScaleOrdinal } from 'd3';

export enum ColorMode {
    BRANCH,
    NODE
}

export enum OpacityMode {
    FLAT,
    GRADIENT
}

export class SunburstRendererColor {
    private static FALLBACK_COLOR = 'rgb(230,247,255)';

    private _colorMode: ColorMode = ColorMode.NODE;
    private _opacityMode: OpacityMode = OpacityMode.GRADIENT;
    private rootColors?: ScaleOrdinal<string, string>;

    constructor(private readonly colors: string[]) {
    }

    set colorMode(value: ColorMode) {
        this._colorMode = value;
    }

    set opacityMode(value: OpacityMode) {
        this._opacityMode = value;
    }

    public getColor(node: TreeData): string {
        switch (this._colorMode) {
            case ColorMode.BRANCH:
                return this.getColorBranch(node);
            default:
            case ColorMode.NODE:
                return this.getColorNode(node);
        }
    }

    public getOpacity(node: TreeData): number {
        switch (this._opacityMode) {
            case OpacityMode.FLAT:
                return SunburstRendererColor.getOpacityFlat(node);
            default:
            case OpacityMode.GRADIENT:
                return SunburstRendererColor.getOpacityFade1(node);
        }
    }

    public setup(root: TreeData): void {
        switch (this._colorMode) {
            case ColorMode.BRANCH:
                this.rootColors = d3.scaleOrdinal(this.colors.slice(0, (root.children?.length || 0) + 1));
                break
            default:
        }
    }

    private getColorBranch(node: TreeData): string {
        let v: typeof node | null = node;
        while (!!v && v.depth > 1) v = v.parent;
        return this.rootColors!(v?.data?.category || '');
    }

    private getColorNode(node: TreeData): string {
        return node.data?.properties?.color || SunburstRendererColor.FALLBACK_COLOR;
    }

    private getColorNodeBlendWithParent(node: TreeData): string {
        const nodeColor = node.data?.properties?.color;
        const parentColor = node.parent?.data?.properties?.color;
        if (nodeColor) {
            return parentColor ? d3.interpolateRgb.gamma(2.2)(nodeColor, parentColor)(.3) : nodeColor;
        }
        return SunburstRendererColor.FALLBACK_COLOR;
    }

    // Fade outer layers / deeper levels
    private static getOpacityFade1(d: TreeData): number {
        return 0.3 + 0.4 / Math.pow(2, d.depth - 1);
    }

    // Fade outer layers / deeper levels (alternative)
    private static getOpacityFade2(d: TreeData): number {
        return 0.1 + 0.5 / Math.pow(2, d.depth - 1);
    }

    // Flat
    private static getOpacityFlat(_d: TreeData): number {
        return .6;
    }

}
