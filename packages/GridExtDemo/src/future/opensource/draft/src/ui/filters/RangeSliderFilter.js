import * as ui from 'datagrok-api/ui';
import {SingleColumnFilter} from "./SingleColumnFilter";
import {ColumnCtx} from "./ColumnCtx";
import {TextUtils} from "../../utils/TextUtils";
import {MathUtils} from "../../utils/MathUtils";
import "../../../styles/styles.css";

export class RangeSliderFilter extends SingleColumnFilter
{
    constructor(ctxCol)
    {
        super(ctxCol);

        this.m_bDate = ctxCol.getType() === DG.COLUMN_TYPE.DATE_TIME;

        this.m_eLabelL = null;
        this.m_eLabelR = null;
        this.m_eCheckNulls = null;

        this.m_nFirstNonNullIndex = ColumnCtx.getFirstNonNullSortedValueIndex(ctxCol);
        this.m_nRowMin = this.m_nFirstNonNullIndex;
        this.m_nRowMax = ctxCol.getLength() -1;

        const arOrder = ctxCol.getSortedOrder();
        const f = ctxCol.getValue(arOrder[this.m_nRowMin]);
        this.m_fMinimum = ctxCol.getValue(arOrder[this.m_nRowMin]);

        let strTypeOf = typeof this.m_fMinimum;
        if(this.m_bDate && (typeof this.m_fMinimum !== "number")) {
            let fTmp = this.m_fMinimum.a;
            if(fTmp === undefined)
                fTmp = this.m_fMinimum.af.a;
            this.m_fMinimum = fTmp;
        }
        strTypeOf = typeof this.m_fMinimum;
        if(strTypeOf !== "number")
            throw new Error("Minimum value must be a number, bu it is " + strTypeOf);

        this.m_fMaximum = ctxCol.getValue(arOrder[this.m_nRowMax]);
        if(this.m_bDate && (typeof this.m_fMaximum !== "number")) {
            let fTmp = this.m_fMaximum.a;
            if(fTmp === undefined)
                fTmp = this.m_fMaximum.af.a;

            this.m_fMaximum = fTmp;
        }

        strTypeOf = typeof this.m_fMaximum;
        if(strTypeOf !== "number")
            throw new Error("Maximum value must be a number, bu it is " + strTypeOf);

        if(this.m_fMinimum > this.m_fMaximum)
            throw new Error("Minimum value (" + this.m_fMinimum + ") cannot exceed maximum value (" + this.m_fMaximum + ")");

        this.m_fMin = this.m_fMinimum;
        this.m_fMax = this.m_fMaximum;

        this.m_bFilterNulls = false;

        this.m_bAPIIsCalling = false;
        this.m_bGUIIsCalling = false;
    }

    isCompatible(ctxCol)
    {
        const b = ctxCol.getType() === DG.COLUMN_TYPE.FLOAT || ctxCol.getType() === DG.COLUMN_TYPE.INT ||
                  ctxCol.getType() === DG.COLUMN_TYPE.DATE_TIME;

        if(!b)
            return false;

        const nFirstNonNullIndex = ColumnCtx.getFirstNonNullSortedValueIndex(ctxCol);
        return  nFirstNonNullIndex >= 0;
    }


    isFiltered(nRecord)
    {
        const ctxCol = this.getColumnCtx();
        const ob = ctxCol.getValue(nRecord);

        if(MathUtils.isNullValue(ob))
          return this.m_bFilterNulls;

        const fVal = this.m_bDate && (typeof ob !== "number") ? ob.a : ob;

        const fMin = this.m_fMin;
        const fMax = this.m_fMax;

        let b = fVal < fMin;
        if(b)
            return true;

        if(fVal > fMax)
          return true;

        return false;
    }

    hasFiltered()
    {
      if(this.m_nRowMin !== this.m_nFirstNonNullIndex)
          return true;

      const ctxCol = this.getColumnCtx();
      if(this.m_nRowMax !== ctxCol.getLength() -1)
        return true;

      return this.getNullsFiltered();
    }


    getNullsFiltered() {return this.m_bFilterNulls;}
    setNullsFiltered(bFiltered)
    {
        this.m_bFilterNulls = bFiltered;
        //my changes this.requestFilter();

        this.m_bAPIIsCalling = true;
        this.m_eCheckNulls.checked = !bFiltered;
        this.m_bAPIIsCalling = false;
    }


    getMin() {return this.m_fMin;}
    getMax() {return this.m_fMax;}
    getMinimum() {return this.m_fMinimum;}
    getMaximum() {return this.m_fMaximum;}

    resetFilter(source)
    {
     this.setValues(this.getMinimum(), this.getMaximum());
     this.setNullsFiltered(false);

     super.resetFilter(source);
    }


    setValues(fMin, fMax)
    {
     const fValueLoOld = this.m_fMin;
     const fValueHiOld = this.m_fMax;

    if(fMin > fMax)
    {
     if(fMin === this.getMin())
       fMax = this.getMin();
     else if(fMax === this.getMax())
      fMin = this.getMax();
     else throw new Error("Min value " + fMin + " cannot exceed the max value " + fMax);
    }

    if(fMin < this.m_fMinimum)
        fMin = this.m_fMinimum;

    if(fMax < fMin)
    fMax = fMin;

    if(fMax > this.m_fMaximum)
     fMax = this.m_fMaximum;

    if(fMin > fMax)
     fMin = fMax;

        //console.log(fMin + " " + fMax);

    let strTextToSet = RangeSliderFilter.formatValue(fMin, this.m_bDate);
    let edit = this.m_eLabelL;


    if(!edit.innerHTML  !== strTextToSet)
      edit.innerHTML = strTextToSet;

    strTextToSet = RangeSliderFilter.formatValue(fMax, this.m_bDate);
    edit = this.m_eLabelR;
    if(!edit.innerHTML !== strTextToSet)
      edit.innerHTML = strTextToSet;

    //Determine directions
    const nDirMin = fMin < this.m_fMin ? -1 : fMin > this.m_fMin ? 1 : 0;
    const nDirMax = fMax < this.m_fMax ? -1 : fMax > this.m_fMax ? 1 : 0;

    if(nDirMin === 0 && nDirMax === 0)
     return;

    this.m_fMin = fMin;
    this.m_fMax = fMax;

    const ctxCol = this.getColumnCtx();
    if(ctxCol.getLength() === 0)
       return;


    if(nDirMin > 0 && nDirMax > 0)
    {
        this.m_nRowMax = RangeSliderFilter.processDoubleMaxValuePlusDir(fMax, this.m_nRowMax, ctxCol, this.m_bDate);
        this.m_nRowMin = RangeSliderFilter.processDoubleMinValuePlusDir(fMin, this.m_nRowMin, ctxCol, this.m_bDate);
    }
    else if(nDirMin < 0 && nDirMax < 0)
    {
        this.m_nRowMin = RangeSliderFilter.processDoubleMinValueMinusDir(fMin, this.m_nRowMin, this.m_nFirstNonNullIndex, ctxCol, this.m_bDate);
        this.m_nRowMax = RangeSliderFilter.processDoubleMaxValueMinusDir(fMax, this.m_nRowMax, this.m_nFirstNonNullIndex, ctxCol, this.m_bDate);
    }
    else
    {
     if(nDirMin > 0)
     {
        try{this.m_nRowMin = RangeSliderFilter.processDoubleMinValuePlusDir(fMin, this.m_nRowMin, ctxCol, this.m_bDate);}
        catch(ex)
        {
            //Afx.TRACE("L+ m_nLowRow= " + m_nLowRow + " " + "m_nHighRow= " + m_nHighRow);
            throw new Error("fValueLoOld " + fValueLoOld + " fValueHiOld " + fValueHiOld +
                "fValueLo " + fMin + " fValueHi " + fMax);
        }


             /*
        if(this.m_fMin === this.m_fMax) {
            let nRec = -1;
            let fV = null;
            let b = false;
            const arOrder = ctxCol.getSortedOrder();
            for(var nR=this.m_nRowMin; nR<=this.m_nRowMax; ++nR)
            {
                nRec = arOrder[nR];
                fV = ctxCol.getValue(nRec);
                b = this.m_fMin <= fV && fV <= this.m_fMax;
                console.log(this.m_fMin + " " + nR + " " + nRec + " " + fV + b);
            }
        }      */
     }
     else if(nDirMin < 0)
     {
        this.m_nRowMin = RangeSliderFilter.processDoubleMinValueMinusDir(fMin, this.m_nRowMin, this.m_nFirstNonNullIndex, ctxCol, this.m_bDate);
        //Afx.TRACE("L- m_nLowRow= " + m_nLowRow + " " + "m_nHighRow= " + m_nHighRow);
     }

    if(nDirMax > 0)
     this.m_nRowMax = RangeSliderFilter.processDoubleMaxValuePlusDir(fMax, this.m_nRowMax, ctxCol, this.m_bDate);

    else if(nDirMax < 0)
     this.m_nRowMax = RangeSliderFilter.processDoubleMaxValueMinusDir(fMax, this.m_nRowMax, this.m_nFirstNonNullIndex, ctxCol, this.m_bDate);
}
//Afx.TRACE("m_nLowRow= " + m_nLowRow + " " + "m_nHighRow= " + m_nHighRow);

/////////my changes changeFilter(iarToShow, iarToHide); moved inside the process.. functions
//updateChangedFilter(null, null);

       //my changes this.requestFilter();

        if(this.m_bGUIIsCalling)
            return;

//update gui
const slider = this.m_slider;
const fK =(slider.maxRange - slider.minRange)/(this.m_fMaximum - this.m_fMinimum);
const nMin  =slider.minRange + Math.floor((fMin - this.m_fMinimum)*fK);
const nMax  =slider.maxRange - Math.floor((this.m_fMaximum- fMax)*fK);

this.m_bAPIIsCalling = true;
this.m_slider.setValues(slider.minRange, slider.maxRange, nMin, nMax);
this.m_bAPIIsCalling = false;
}




    getShowNullValues() {return this.m_bShowNulls;}

    parseInput(strText)
    {
        let fValue = NaN;
        if(this.m_bDate)
         fValue = Date.parse(strText)
        else
         fValue = Number(strText);

        return fValue;
    }

    render()
    {
        super.render();

        const ctxCol = this.getColumnCtx();
        const nRecordCount = ctxCol.getLength();
        this.m_slider = ui.rangeSlider(0, 500, 0, 500);

        const filterThis = this;

        this.m_slider.onValuesChanged.subscribe(function (e) {

            if(filterThis.m_bAPIIsCalling)
                return;

            const slider = filterThis.m_slider;
            const fK = (filterThis.m_fMaximum - filterThis.m_fMinimum)/(slider.maxRange - slider.minRange);
            let fLow  = slider.min === slider.maxRange ? filterThis.m_fMaximum :
                filterThis.m_fMinimum + (slider.min - slider.minRange)*fK;
            const fHigh = slider.max === slider.minRange ? filterThis.m_fMinimum :
                filterThis.m_fMaximum - (slider.maxRange - slider.max)*fK;

            if(fLow > fHigh)
                fLow = fHigh;

            //Afx.TRACE("range_slider: " + fLow + " " + fHigh);
            filterThis.m_bGUIIsCalling = true;
            filterThis.setValues(fLow, fHigh);
            filterThis.m_bGUIIsCalling = false;
            filterThis.requestFilter();

        });



        const eDivHeader = document.createElement("div");
        eDivHeader.style.width = "100%";
        eDivHeader.style.whiteSpace = "nowrap";

        const eLblL = document.createElement("div");
        eLblL.setAttribute("contenteditable", "true");


        const nRowMin = ctxCol.getSortedOrder()[this.m_nRowMin];
        let fValue = ctxCol.getValue(nRowMin);
        if(this.m_bDate && (typeof fValue !== "number"))
            fValue = fValue.a;
        const strValMin = RangeSliderFilter.formatValue(fValue, this.m_bDate);
        eLblL.innerHTML = strValMin;
        eLblL.style.width = "40%";
        eLblL.style.overflow = "hidden";
        //eLblL.style.whiteSpace = "nowrap";
        eLblL.style.display = "inline-block";

        eDivHeader.appendChild(eLblL);

        eLblL.addEventListener('keydown', function(e) {
            const strCode = e.code;
            if (strCode === "Enter")
                e.preventDefault();
        });

        eLblL.addEventListener('keyup', function(e){

            e.preventDefault();

            const strText =  eLblL.innerHTML;
            const fVal = filterThis.parseInput(strText);
            const bValid = !isNaN(fVal);
            eLblL.style.color = bValid ? "black" : "red";

            const strCode = e.code;
            if(strCode === "Enter")
            {
               if(!bValid)
                   return;

                eLblL.blur();
                filterThis.setValues(fVal, filterThis.getMax());
                filterThis.requestFilter();
            }
         });





        this.m_eLabelL = eLblL;

        const eLblR = document.createElement("div");
        eLblR.setAttribute("contenteditable", "true");

        const nRowMax = ctxCol.getSortedOrder()[this.m_nRowMax];
        fValue = ctxCol.getValue(nRowMax);
        if(this.m_bDate  && (typeof fValue !== "number"))
              fValue = fValue.a;

        const strValMax = RangeSliderFilter.formatValue(fValue, this.m_bDate);

        eLblR.innerHTML = strValMax;
        eLblR.style.width = "40%";
        eLblR.style.display = "inline-block";
        eLblR.style.overflow = "hidden";
        eDivHeader.appendChild(eLblR);


        eLblR.addEventListener('keydown', function(e) {
            const strCode = e.code;
            if (strCode === "Enter")
                e.preventDefault();
        });

        eLblR.addEventListener('keyup', function(e){

            e.preventDefault();

            const strText =  eLblR.innerHTML;
            const fVal = filterThis.parseInput(strText);
            const bValid = !isNaN(fVal);
            eLblR.style.color = bValid ? "black" : "red";

            const strCode = e.code;
            if(strCode === "Enter")
            {
                if(!bValid)
                    return;

                eLblR.blur();
                filterThis.setValues(filterThis.getMin(), fVal);
                filterThis.requestFilter();
            }
        });




        this.m_eLabelR = eLblR;

        const eCBNull = document.createElement("INPUT");
        eCBNull.style.width = "20px";
       eCBNull.checked = !this.m_bFilterNulls;
        eCBNull.setAttribute("type", "checkbox");

        eCBNull.addEventListener('change', (event) => {

            if(filterThis.m_bAPIIsCalling)
                return;

            const eCBTmp =  event.currentTarget;
            const bSelected = eCBTmp.checked;
            filterThis.setNullsFiltered(!bSelected);
            filterThis.requestFilter();
        });


        eDivHeader.appendChild(eCBNull);
        this.m_eCheckNulls = eCBNull;


        this.addUIItem(eDivHeader);
        this.addUIItem(this.m_slider.root);


        const lstChildren = this.m_slider.root.childNodes;
        for(var n=0; n<lstChildren.length; ++n)
        {
            if(lstChildren[n].tagName === "svg")
            {
                lstChildren[n].setAttribute("width", "100%");
                lstChildren[n].setAttribute("height", "15px");
                lstChildren[n].style.overflow = "hidden";
            }
        }

        //this.slider.root.style.width = "75px";
        //this.root.style.height = "55px";
        this.root.style.overflow = "auto";
    }

}

RangeSliderFilter.formatValue = function (obValue, bDate)
{
    if(bDate)
    {
        const nTime = obValue;//.a;
        return TextUtils.formatDate(nTime);
    }

    return obValue.toFixed(1);
}

RangeSliderFilter.processDoubleMinValuePlusDir = function(fLow, nLowIndex, ctxCol, bDate)
{
    const arOrder = ctxCol.getSortedOrder();
    let fValue = NaN;
    let nRow = -1;
    for(;nLowIndex<arOrder.length; ++nLowIndex) //we don't know new high index
    {
        nRow = arOrder[nLowIndex];
        fValue = ctxCol.getValue(nRow);//tdm.getDoubleValue(nRow, nColumn);
        if(bDate && (typeof fValue !== "number"))
            fValue = fValue.a;

        //console.log(nLowIndex + " " + fValue);

        if (isNaN(fValue) || MathUtils.isNullValue(fValue))
        {
         let b = isNaN(fValue);
         b = MathUtils.isNullValue(fValue);
         throw new Error("Double value cannot be null (" + fValue + ") at  index (" + nLowIndex + ")");
    }
        if(fLow <= fValue)
            break;

        ///////if(iarToRem != null)
        //////iarToRem.append(nRow);
        //if(device != null)
          //  device.changeFilter(nRow, true, false);//my changes
    }
    //m_nLowIndex may be equal to arRows.length
    return nLowIndex;
}

RangeSliderFilter.processDoubleMinValueMinusDir = function(fLow, nLowIndex, nFirstNonNullIndex, ctxCol, bDate)
{
    const arOrder = ctxCol.getSortedOrder();
    let fValue = NaN;
    let nRow = -1;
    for(--nLowIndex; nLowIndex>=nFirstNonNullIndex; --nLowIndex)
    {
        nRow = arOrder[nLowIndex];
        fValue =  ctxCol.getValue(nRow);
        if(bDate && (typeof fValue !== "number"))
            fValue = fValue.a;

        if(MathUtils.isNullValue(fValue))
            throw new Error("Double value cannot be null (" + fValue + ") at  index (" + nLowIndex + ")");

        if(fLow > fValue)
            break;

        //if(device != null)
          //  device.changeFilter(nRow, false, false);
    }

    return ++nLowIndex;
}



RangeSliderFilter.processDoubleMaxValuePlusDir = function(fHigh, nHighIndex, ctxCol, bDate)
{
    const arOrder = ctxCol.getSortedOrder();
    let fValue = NaN;
    let nRow = -1;
    for(++nHighIndex; nHighIndex<arOrder.length; ++nHighIndex)
    {
        nRow = arOrder[nHighIndex];
        fValue = ctxCol.getValue(nRow);
        if(bDate && (typeof fValue !== "number"))
            fValue = fValue.a;

        if(isNaN(fValue) || MathUtils.isNullValue(fValue))
            throw new Error("Double value cannot be null (" + fValue + ") at  index (" + nLowIndex + ")");

        if(fHigh < fValue)
            break;

        //if(device != null)
          //  device.changeFilter(nRow, false, false);
    }

    return --nHighIndex;
}

RangeSliderFilter.processDoubleMaxValueMinusDir = function(fHigh, nHighIndex, nFirstNonNullndex, ctxCol, bDate)
{
    const arOrder = ctxCol.getSortedOrder();
    let fValue = NaN;
    let nRow = -1;
    for(;nHighIndex>=nFirstNonNullndex; --nHighIndex)
    {
        nRow = arOrder[nHighIndex];
        fValue = ctxCol.getValue(nRow);
        if(bDate && (typeof fValue !== "number"))
            fValue = fValue.a;

        if(isNaN(fValue) || MathUtils.isNullValue(fValue))
            throw new Error("Double value cannot be null (" + fValue + ") at  index (" + nLowIndex + ")");

        if(fHigh >= fValue)
            break;

        //if(device != null)
            //device.changeFilter(nRow, true, false);
    }
    return nHighIndex;
}
