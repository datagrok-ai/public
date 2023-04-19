
export interface ResidueSelectionInterface {
    modelId: string;
    labelAsymId: string;
    seqIds: Set<number>;
}
export interface RegionSelectionInterface{
    begin:number;
    end:number;
    isEmpty?: boolean;
    source:'structure'|'sequence';
}

export interface ChainSelectionInterface {
    modelId: string;
    labelAsymId: string;
    regions: Array<RegionSelectionInterface>;
}

export class RcsbFvSelectorManager {

    private lastSelection: ChainSelectionInterface | null = null;
    private selection: Array<ChainSelectionInterface> = new Array<ChainSelectionInterface>();
    private hover: Array<ChainSelectionInterface> = new Array<ChainSelectionInterface>();

    public setSelectionFromRegion(modelId: string, labelAsymId: string, region: RegionSelectionInterface, mode:'select'|'hover'): void {
        this.clearSelection(mode);
        this.addSelectionFromRegion(modelId, labelAsymId, region, mode);
    }

    public addSelectionFromRegion(modelId: string, labelAsymId: string, region: RegionSelectionInterface, mode:'select'|'hover'): void {
        if(mode === 'select'){
            this.selection.push({modelId:modelId, labelAsymId:labelAsymId, regions:[region]});
        }else{
            this.hover.push({modelId:modelId, labelAsymId:labelAsymId, regions:[region]});
        }
    }

    public setSelectionFromMultipleRegions(regions: {modelId: string, labelAsymId: string, region: RegionSelectionInterface}[], mode:'select'|'hover'): void {
        this.clearSelection(mode);
        this.addSelectionFromMultipleRegions(regions, mode);
    }

    public addSelectionFromMultipleRegions(regions: {modelId: string, labelAsymId: string, region: RegionSelectionInterface}[], mode:'select'|'hover'): void {
        regions.forEach(r=>{
            this.addSelectionFromRegion(r.modelId, r.labelAsymId, r.region, mode);
        });
    }

    public setSelectionFromResidueSelection(res: Array<ResidueSelectionInterface>, mode:'select'|'hover', source: 'structure'|'sequence'): void {
        if(mode==='select'){
            this.selection = selectionFromResidueSelection(res, mode, source);
        }else{
            this.hover = selectionFromResidueSelection(res, mode, source);
        }
    }

    public getSelection(mode:'select'|'hover'): Array<ChainSelectionInterface> {
        if(mode === 'select')
            return this.selection;
        else
            return this.hover;
    }

    public getLastSelection(mode:'select'|'hover'): ChainSelectionInterface | null{
       return this.lastSelection;
    }

    public setLastSelection(mode:'select'|'hover', selection: ChainSelectionInterface | null): void {
        this.lastSelection = selection;
    }

    public getSelectionWithCondition(modelId: string, labelAsymId: string, mode:'select'|'hover'): ChainSelectionInterface | undefined{
        const sel: Array<ChainSelectionInterface> = mode === 'select' ?
            this.selection.filter(d=>(d.modelId===modelId && d.labelAsymId === labelAsymId)) :
            this.hover.filter(d=>(d.modelId===modelId && d.labelAsymId === labelAsymId));
        if(sel.length > 0)
            return {
                modelId: sel[0].modelId,
                labelAsymId: sel[0].labelAsymId,
                regions: ([] as RegionSelectionInterface[]).concat.apply([], sel.map(s => s.regions))
            };
    }

    public clearSelection(mode:'select'|'hover', selection?:{modelId?:string, labelAsymId?: string}): void {
        if(!selection)
            if(mode === 'select')
                this.selection = new Array<ChainSelectionInterface>();
            else
                this.hover = new Array<ChainSelectionInterface>();
        else
            if(selection.labelAsymId || selection.modelId)
                if(mode === 'select')
                    this.selection = this.selection.filter(r=>selectionFilter(r, selection));
                else
                    this.hover = this.hover.filter(r=>selectionFilter(r, selection));
    }

    public selectionSource(mode:'select'|'hover', region:{modelId:string;labelAsymId:string;begin:number;end:number;}): 'structure'|'sequence'|undefined{
        if(mode === 'select')
            return this.selection
                .filter(r=>(r.modelId === region.modelId && r.labelAsymId === region.labelAsymId))[0]?.regions
                .filter(r=>(r.begin === region. begin && r.end === region.end))[0]?.source;
        else
            return this.hover
                .filter(r=>(r.modelId === region.modelId && r.labelAsymId === region.labelAsymId))[0]?.regions
                .filter(r=>(r.begin === region. begin && r.end === region.end))[0]?.source;
    }
}

function selectionFromResidueSelection(res: Array<ResidueSelectionInterface>, mode:'select'|'hover', source: 'structure'|'sequence'): Array<ChainSelectionInterface> {
    const selMap: Map<string,Map<string,Set<number>>> = new Map<string, Map<string, Set<number>>>();
    res.forEach(r=>{
        if(!selMap.has(r.modelId))
            selMap.set(r.modelId,new Map<string, Set<number>>());
        if(!selMap.get(r.modelId)!.has(r.labelAsymId))
            selMap.get(r.modelId)!.set(r.labelAsymId, new Set<number>());
        r.seqIds.forEach(s=>{
            selMap.get(r.modelId)!.get(r.labelAsymId)!.add(s);
        })
    });
    const selection = new Array<ChainSelectionInterface>();
    selMap.forEach((labelMap, modelId)=>{
        labelMap.forEach((seqSet,labelId)=>{
            selection.push({modelId:modelId, labelAsymId: labelId, regions:buildIntervals(seqSet, source)});
        });
    });
    return selection;
}

function buildIntervals(ids: Set<number>, source: 'structure'|'sequence'): Array<RegionSelectionInterface>{
    const out: Array<RegionSelectionInterface> = new Array<RegionSelectionInterface>();
    const sorted: Array<number> = Array.from(ids).sort((a,b)=>(a-b));
    let begin: number = sorted.shift()!;
    let end: number = begin;
    for(const n of sorted){
        if(n==(end+1)){
            end = n;
        }else{
            out.push({begin:begin,end:end,source:source});
            begin = n;
            end = n;
        }
    }
    out.push({begin:begin,end:end,source:source});
    return out;
}

function selectionFilter(r:ChainSelectionInterface, selection:{modelId?:string, labelAsymId?: string}): boolean{
    if(selection.modelId && selection.labelAsymId)
        return (r.modelId != selection.modelId || r.labelAsymId != selection.labelAsymId);
    else if(selection.modelId)
        return (r.modelId != selection.modelId);
    else if(selection.labelAsymId)
        return (r.labelAsymId != selection.labelAsymId);
    else
        return false;
}
