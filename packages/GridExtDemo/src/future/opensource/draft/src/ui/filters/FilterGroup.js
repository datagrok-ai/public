import "../../../styles/styles.css";
import {SingleColumnFilter} from "./SingleColumnFilter";

export class FilterGroup extends SingleColumnFilter
{
  constructor(ctxCol)
  {
      super(ctxCol);

      const eDivRoot = this.root;
      eDivRoot.style.paddingTop = "15px";

      const arFilters = this.createFilters(ctxCol);

      const nFilterCount = arFilters.length;
      for(var n=0; n<nFilterCount; ++n)
      {
          arFilters[n].m_parent = this;
          arFilters[n].root.classList.value = "";
      }

      this.m_arFilters = arFilters;
  }

  isCompatible(ctxCol)
  {
       throw new NotImplementedError();
  }

  createFilters(ctxCol)
  {
      throw new NotImplementedError();
  }

  getChildFilterCount()
  {
      return this.m_arFilters.length;
  }

  getChildFilter(nChild)
  {
      return this.m_arFilters[nChild];
  }


     /*
    isCloseButtonEnabled()
    {
        return true;
    }

  /*

  close()
  {

  }


  closeChildFilter(nChild)
  {

  } */


  isFiltered(nRecord)
  {
      const nFilterCount = this.m_arFilters.length;
      for(var n=0; n<nFilterCount; ++n)
      {
          if(this.m_arFilters[n].isFiltered(nRecord))
              return true;
      }

      return false;
  }

  hasFiltered()
  {
      const nFilterCount = this.m_arFilters.length;
      for(var n=0; n<nFilterCount; ++n)
      {
          if(this.m_arFilters[n].hasFiltered())
              return true;
      }

      return false;
  }


  onBeforeApplyFilter()
  {
      super.onBeforeApplyFilter();

      const nFilterCount = this.m_arFilters.length;
      for(var n=0; n<nFilterCount; ++n)
      {
          this.m_arFilters[n].onBeforeApplyFilter();
      }
  }

  onAfterApplyFilter()
  {
      super.onAfterApplyFilter();

      const nFilterCount = this.m_arFilters.length;
      for(var n=0; n<nFilterCount; ++n)
      {
          this.m_arFilters[n].onAfterApplyFilter();
      }
  }

  resetFilter(source)
  {
      let filter = null;
      const nFilterCount = this.getChildFilterCount();
      for(var n=0; n<nFilterCount; ++n)
      {
          filter = this.getChildFilter(n);
          filter.resetFilter(this);
      }

      super.resetFilter(source);
  }

  render()
  {
      super.render();

      this.root.childNodes[0].style.padding = "0px";

      const nFilterCount = this.m_arFilters.length;
      //const ar = new Array(nFilterCount);
       for(var n=0; n<nFilterCount; ++n)
      {
          //this.m_arFilters[n].root.style.marginLeft = "0px";
          this.m_arFilters[n].render();
          this.addUIItem(this.m_arFilters[n].root);
      }
  }
}
