import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {FiltersUtils} from "./FiltersUtils";
import {ErrorUtils} from "../../utils/ErrorUtils";

export class FilterPanel extends DG.Filter
{
    constructor(arFilterClasses)
    {
        super();
        this.root = ui.divV(null);//, 'filter-panel');
        this.root.classList.add("d4-flex-col");
        this.root.classList.add("d4-filter");

        //const c = document.getElementsByClassName("panel-base splitter-container-horizontal");
        //if(c.length > 0)
          //  c.item(0).style.width = "250px";

        this.m_arFilterClasses = arFilterClasses;
        this.m_arFilters = null;
        this.m_dframe = null;
        this.m_bDetaching = false;

        this.m_fnOnRowsFiltering = null;
    }

    addFilter(col)
    {
        ErrorUtils.verifyClass(col, DG.Column);

        const strColName = col.name;
        let filter = this.getFilter(strColName);
        if(filter !== null)
            throw new Error("Filter already exists for column '" + strColName + "'");

        filter = FiltersUtils.createFilter(col, this.m_arFilterClasses);
        if(filter !== null)
        {
            this.m_arFilters.push(filter);
            filter.m_parent = this;
            filter.attach(col.dataFrame);//will call render()
            const viewerFilter = FiltersUtils.getFiltersViewer();
            viewerFilter.root.appendChild(filter.root);
            filter.m_nId = FilterPanel.FILER_INSTANCE_COUNT;
            ++FilterPanel.FILER_INSTANCE_COUNT;
            return filter;
        }

        return null;
    }

    removeFilter(strColName)
    {
        ErrorUtils.verifyType(strColName, String);

        const filter = this.getFilter(strColName);
        const nFilter = this.m_arFilters.indexOf(filter);
        if (nFilter < 0)
            throw new Error("Never should happen.");

        this.m_arFilters.splice(nFilter, 1);
        filter.detach();
        const viewerFilter = FiltersUtils.getFiltersViewer();
        viewerFilter.root.removeChild(filter.root);

        if(filter.hasFiltered())
          this.requestFilter();
    }

    getFilter(strColName)
    {
        ErrorUtils.verifyType(strColName, String);
        const nFilterCount = this.m_arFilters.length;

        let ctxCol = null;
        let filter = null;
        for(let nFilter=0; nFilter<nFilterCount; ++nFilter)
        {
           filter = this.m_arFilters[nFilter];
           ctxCol = filter.getColumnCtx();

           if(strColName === ctxCol.getName())
               return filter;
        }

        return null;
    }


    attach(dframe)
    {
        ErrorUtils.verifyClass(dframe, DG.DataFrame);

        /*
        const viewTable = grok.shell.v;//getTableView(dframe.name);
        const lstViewers = viewTable.viewers;
        for(var viewer of lstViewers)
        {
             if(viewer.type === DG.VIEWER.FILTERS)
              viewer.d.m_panelFilters = this;
        } */

       this.m_dframe = dframe;
       this.m_arFilters = [];

      this.m_fnOnRowsFiltering = this.m_dframe.onRowsFiltering.subscribe((_) => this.applyFilter());

      this.render();
    }

    detach()
    {
        this.m_bDetaching = true;
        this.requestFilter();
        //this.m_bDetaching = false;

        let filter = null;
        while(this.m_arFilters.length > 0)
        {
            filter = this.m_arFilters.pop();
        }

        this?.m_fnOnRowsFiltering?.unsubscribe();

        this.m_arFilters = null;
        this.m_dframe = null;
    }

    requestFilter()
    {
        this.m_dframe.rows.requestFilter();
    }

    applyFilter()
    {
        if(this.m_bDetaching)
          return;

        let nFilter = -1;
        let filter = null;

        const nRecordCount = this.m_dframe.rowCount;
        const bitset = this.m_dframe.filter;
        const nFilterCount = this.m_arFilters.length;

        for (nFilter = 0; nFilter < nFilterCount; ++nFilter)
        {
            filter = this.m_arFilters[nFilter];
            filter.onBeforeApplyFilter();
        }

         let b = false;
        for(let nRecord=0; nRecord<nRecordCount; ++nRecord)
        {
            for (nFilter = 0; nFilter < nFilterCount; ++nFilter)
            {
                filter = this.m_arFilters[nFilter];
                b = filter.isFiltered(nRecord);
                if(b)
                {
                    bitset.set(nRecord, false, false);
                    break;
                }
            }
        }

        for (nFilter = 0; nFilter < nFilterCount; ++nFilter)
        {
            filter = this.m_arFilters[nFilter];
            filter.onAfterApplyFilter();
        }

        this.m_dframe.filter.fireChanged();
    }

    hasFiltered()
    {
        let filter = null;
        let b = false;
        const nFilterCount = this.m_arFilters.length;
        for(var nFilter = 0; nFilter < nFilterCount; ++nFilter)
        {
            filter = this.m_arFilters[nFilter];
            b = filter.isFiltered(nRecord);
            if(b)
              return true;

        }
        return b;
    }

    render()
    {
       $(this.root).empty();

       const viewerFilter = FiltersUtils.getFiltersViewer(this.m_dframe);

       const nFilterCount = this.m_arFilters.length;
       const ar = new Array(nFilterCount);
       for(var nFilter=0; nFilter<nFilterCount; ++nFilter)
       {
         ar[nFilter] = this.m_arFilters[nFilter].root;
         viewerFilter.root.appendChild(this.m_arFilters[nFilter].root);
           this.m_arFilters[nFilter].render();
       }

        //const eSplit = ui.splitV(ar);
        //this.root.appendChild(eSplit);
    }
}

FilterPanel.FILER_INSTANCE_COUNT = 0;
