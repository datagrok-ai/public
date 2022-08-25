import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {IncompatibleFilterError} from "./IncompatibleFilterError";
import {FiltersUtils} from "./FiltersUtils";
import ImageReset from "../../../images/reset.png";
import ImageClose from "../../../images/close.png";
import {NotImplementedError} from "../../lang/NotImplementedError";
import {ErrorUtils} from "../../utils/ErrorUtils";
import {ColumnCtx} from "./ColumnCtx";

export class SingleColumnFilter extends DG.Filter
{
    constructor(ctxCol)
    {
        super();

        ErrorUtils.verifyClass(ctxCol, ColumnCtx);

        if(!this.isCompatible(ctxCol))
         throw new IncompatibleFilterError("Filter '" + this.constructor.name + "' is not compatible with the column " + ctxCol.name);

        this.m_column = ctxCol;
        const eDivRoot = document.createElement("div");
        //eFieldSet.classList.add('single-col-query-device');
        eDivRoot.classList.add("d4-flex-col");
        eDivRoot.classList.add("d4-filter");
        eDivRoot.classList.add('y-axis-layout');

        ///////////eDivRoot.style.marginBottom = "10px";
        //eFieldSet.style.marginRight = "20px";
        eDivRoot.style.paddingTop = "10px";
        ////////////eDivRoot.style.paddingBottom = "10px";
        //eDivRoot.style.height = "35px";
        eDivRoot.style.width = "100%";
        this.root = eDivRoot;

        this.m_nId = -1;
        this.m_parent = null;
        this.m_eBtnReset = null;
        this.m_eBtnClose = null;
        this.m_bFloating = false;

        this.m_fnOnRowsFiltering = null;

    }

    getId() {return this.m_nId;}

    getColumnCtx() {return this.m_column;}

    isCompatible(ctxCol)
    {
        throw new NotImplementedError();
    }

    attach(dframe)
    {
        ErrorUtils.verifyClass(dframe, DG.DataFrame);

        if(this.m_parent === null)
          this.m_fnOnRowsFiltering = dframe.onRowsFiltering.subscribe((_) => this.applyFilter());

        this.render();
    }

    detach()
    {
        if(this.m_parent === null) {
            this?.m_fnOnRowsFiltering?.unsubscribe();
            this.m_fnOnRowsFiltering = null;
        }

    }

    isFiltered(nRecord)
    {
        throw new NotImplementedError();
    }

    hasFiltered()
    {
        throw new NotImplementedError();
    }

    onBeforeApplyFilter() {}
    onAfterApplyFilter() {}


    requestFilter()
    {
        //filterThis.getColumnCtx().dataFrame.rows.requestFilter();
        this.m_parent.requestFilter();

        const b = this.hasFiltered();
        this.m_eBtnReset.style.visibility = b ? "visible" : "hidden";
    }


    close()
    {
        const ctxCol = this.getColumnCtx();
        this.m_parent.removeFilter(ctxCol.getName());
    }

    isFloating()
    {
        return this.m_bFloating;
    }

    float(nX, nY)
    {
        ErrorUtils.verifyType(nX, Number);
        ErrorUtils.verifyType(nY, Number);

        if(SingleColumnFilter.m_filterFloating !== null)
         SingleColumnFilter.m_filterFloating.dock();

        const ctxCol = this.getColumnCtx();
        const dframe = ctxCol.getDataFrame();
        let viewerFilters = FiltersUtils.getFiltersViewer();
        if(viewerFilters === null)
        {
            viewerFilters = FiltersUtils.openFilterViewer(dframe);
            FiltersUtils.showFilters(false);
        }

        this.m_bFloating = true;

        if(this.isCloseButtonEnabled())
         this.m_eBtnClose.style.visibility = "hidden";

        const dial = ui.dialog("Filters");

        const eDiv = this.root;
        eDiv.style.display = "block";
        eDiv.classList.value = "";
        dial.add(eDiv);

        const filterThis = this;

        dial.onCancel(function (e) {
            if(!dial.dart.m_bVisible)
                return;

        dial.dart.m_bVisible = false;
        const b = filterThis.m_bFloating;
        filterThis.dock();
            // viewerFilters.root.appendChild(eDiv);
        });

        dial.onClose.subscribe(function() {
            if(!dial.dart.m_bVisible)
                return;

            dial.dart.m_bVisible = false;
            filterThis.dock();
            //viewerFilters.root.appendChild(eDiv);
        });


        const nHRoot = this.root.clientHeight;
        dial.show({modal: false, x: nX, y: nY, width: 250, height: nHRoot + 100});
        dial.dart.m_bVisible = true;


        SingleColumnFilter.m_dlgFloating = dial;
        SingleColumnFilter.m_filterFloating = this;
    }

    dock()
    {
       if(!this.m_bFloating)
            throw new Error("The filter is not floating.");

        const ctxCol = this.getColumnCtx();
        const dframe = ctxCol.getDataFrame();
        const viewerFilters = FiltersUtils.getFiltersViewer(dframe);
        if(viewerFilters === null)
            throw new Error("Filters viewer cannot be null.");

        const eDiv = this.root;
        viewerFilters.root.appendChild(eDiv);


        if(this.isCloseButtonEnabled())
            this.m_eBtnClose.style.visibility = "visible";

        if(SingleColumnFilter.m_dlgFloating !== null && SingleColumnFilter.m_dlgFloating.dart.m_bVisible) {
            SingleColumnFilter.m_dlgFloating.dart.m_bVisible = false;
            SingleColumnFilter.m_dlgFloating.close();

        }
        this.m_bFloating = false;
        SingleColumnFilter.m_dlgFloating = null;
        SingleColumnFilter.m_filterFloating = null;
    }


    applyFilter() {}

    resetFilter(source)
    {
        if(source === this)
            this.requestFilter();
        else
            this.m_eBtnReset.style.visibility = "hidden";
    }

    isCloseButtonEnabled()
    {
     return this.m_parent === null || !(this.m_parent instanceof SingleColumnFilter);
    }

    render()
    {
      $(this.root).empty();

      let img = new Image();
      img.src = ImageReset;

      const filterThis = this;
      const btnReset = ui.button("R", function() {
            filterThis.resetFilter(filterThis);
      }, "Reset");


      btnReset.style.position = "relative";

      btnReset.style.width = "15px";
      btnReset.style.height = "15px";
      //btnReset.style.left = "97%";
      btnReset.style.float = "right";
      btnReset.style.top = "-25px";
      btnReset.innerHTML = img.outerHTML;
      btnReset.style.visibility = "hidden";
      this.m_eBtnReset = btnReset;
      //btnReset.classList.add("fa-sync");
       // btnReset.classList.add("grok-icon");
      //  btnReset.classList.add("fal");

      const eFieldSet = document.createElement("fieldset");
      const eLegend = document.createElement("legend");
      eLegend.style.fontWeight = "bold";
      eLegend.innerHTML = this.getColumnCtx().getName();



      //Close Button
      if(this.isCloseButtonEnabled())
      {
          img = new Image();
          img.src = ImageClose;

            const btnClose = ui.button("C", function () {
                  filterThis.close();
              }, "Close");


            btnClose.style.position = "relative";

            btnClose.style.width = "15px";
            btnClose.style.height = "15px";
            //btnReset.style.left = "97%";
            btnClose.style.float = "right";
            btnClose.style.top = "-25px";
            btnClose.innerHTML = img.outerHTML;

            this.m_eBtnClose = btnClose;
            eFieldSet.appendChild(btnClose);
      }

      eFieldSet.appendChild(btnReset);
      eFieldSet.appendChild(eLegend);
      //eFieldSet.style.height="100%";
      //eFieldSet.style.overflow = "auto";

      this.root.appendChild(eFieldSet);
    }

    addUIItem(eDOMElement)
    {
        const eFieldSet = this.root.childNodes[0];
        eFieldSet.appendChild(eDOMElement);
    }
}

SingleColumnFilter.INSTAMCE_COUNT = 0;
SingleColumnFilter.m_dlgFloating = null;
SingleColumnFilter.m_filterFloating = null;
