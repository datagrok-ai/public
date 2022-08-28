import {ClassMap} from "../lang/ClassMap";

export class GridRenderersConfig
{
    constructor()
    {
        this.m_mapCellRenderers = new ClassMap();
        this.m_mapColHeaderRenderers = new ClassMap();
        this.m_rendererRowHeader = null;
        this.m_rendererColRowHeader = null;
        this.m_rendererColHeaderDefault = null;
    }

    getDefaultColHeaderRenderer() {return this.m_rendererColHeaderDefault; }
    setDeffaultColHeaderRenderer(renderer)
    {
        this.m_rendererColHeaderDefault = renderer;
    }



    getColRowHeaderRenderer() {return this.m_rendererColRowHeader; }
    setColRowHeaderRenderer(renderer)
    {
        this.m_rendererColRowHeader = renderer;
    }

    getRowHeaderRenderer() {return this.m_rendererRowHeader; }
    setRowHeaderRenderer(renderer)
    {
        this.m_rendererRowHeader = renderer;
    }

    getCellRenderer(classSemType)
    {
        return this.m_mapCellRenderers.get(classSemType);
    }

    setCellRenderer(classSemType, renderer)
    {
        this.m_mapCellRenderers.set(classSemType, renderer);
    }

    getColHeaderRenderer(classSemType)
    {
        return this.m_mapColHeaderRenderers.get(classSemType);
    }

    setColHeaderRenderer(classSemType, renderer)
    {
        this.m_mapColHeaderRenderers.set(classSemType, renderer);
    }
}