
import {AllRangeCondition} from "./RangeCondition";
import {RangeSliderFilter} from "./RangeSliderFilter";

export class CombinedRangeSlider extends RangeSliderFilter
{

constructor(ctxCol, arConditions)
{
    super(ctxCol);

    if(arConditions === null || arConditions === undefined)
        throw new Error("Array containing conditions cannot be null.");

    this.m_arConditions = [];
    this.m_arMinMaxs = [NaN, NaN];

    const fMinimum = this.getMinimum();
    const fMaximum = this.getMaximum();

   // const cbb = new JComboBox();
    //cbb.setRenderer(new ThisComboRenderer(cbb.getRenderer()));

    const conditionAll = new AllRangeCondition();
    conditionAll.m_bIncludeNulls = true;

    for(var nC = 0; nC < arConditions.length; ++nC)
    {
        if(arConditions[nC] == null)
            throw new Error("Condition at index (" + nC + ") cannot be null.");

        if(arConditions[nC].fillRangeBounds(fMinimum, fMaximum, this.m_arMinMaxs))
        {

            if(!Number.isFinite(this.m_arMinMaxs[0]) || isNaN(this.m_arMinMaxs[0]) ||
               !Number.isFinite(this.m_arMinMaxs[1]) || isNaN(this.m_arMinMaxs[1]) ||
                this.m_arMinMaxs[0] > this.m_arMinMaxs[1] ||
                this.m_arMinMaxs[0] < fMinimum || this.m_arMinMaxs[0] > fMaximum ||
                this.m_arMinMaxs[1] < fMinimum || this.m_arMinMaxs[1] > fMaximum)
                throw new Error("Condition returned invalid parameters");

            if(this.m_arMinMaxs[0] === fMinimum && this.m_arMinMaxs[1] === fMaximum && nC < arConditions.length -1)
                continue;

            this.m_arConditions.push(arConditions[nC]);
            //cbb.addItem(arConditions[nC]);
        }
    }

    this.m_arConditions.push(conditionAll);
    this.m_nSelCondIndex = this.m_arConditions.length -1;//-1

    this.m_eSelectCond = null;
    //cbb.addItem(conditionAll);

    //setComboBox(cbb);
  }


    getConditionCount()
    {
        return this.m_arConditions.length;
    }

    getSelectedConditionIndex()
    {
     return this.m_nSelCondIndex;
    }

    setSelectedConditionIndex(nIndex)
    {
     if(nIndex >=  this.m_arConditions.length)
      throw new Error("Index is out of bounds " + nIndex);

     const nIndexOld = this.m_nSelCondIndex;
     if(nIndexOld === nIndex)
      return;

    this.m_nSelCondIndex = nIndex;

    if(nIndex >=0 && !this.m_arConditions[nIndex].fillRangeBounds(this.getMinimum(), this.getMaximum(), this.m_arMinMaxs))
     throw new Error("Never should happen.");

    this.setValues(this.m_arMinMaxs[0], this.m_arMinMaxs[1]);
     /////////////setShowNullValues(m_listConditions.get(nIndex).m_bIncludeNulls, true);

    this.m_bAPIIsCalling = true;
    this.m_eSelectCond.selectedIndex = nIndex;
    this.m_bAPIIsCalling = false;
}

setValues( fValueLo, fValueHi)
{
    super.setValues(fValueLo, fValueHi);

    if(this.m_bAPIIsCalling)
        return;

    this.m_bAPIIsCalling = true;
    //const cbb = getComboBox();
    if(this.m_eSelectCond.selectedIndex !== -1)//cbb.getSelectedIndex() !== -1)
    {
        this.m_eSelectCond.selectedIndex = -1;
        //cbb.setSelectedIndex(-1);
        this.m_nSelCondIndex = -1;
    }
    this.m_bAPIIsCalling = false;
}


 resetFilter(source)
 {
  this.setSelectedConditionIndex(-1);

  super.resetFilter(source);
 }

 render()
 {
    super.render();

     const nCondCount = this.m_arConditions.length;

     const eDivCond = document.createElement("div");
     eDivCond.style.width = "100%";
     eDivCond.style.whiteSpace = "nowrap";
     eDivCond.style.display = "inline-block";
     eDivCond.style.overflowX = "auto";
     eDivCond.style.overflowY = "hidden";


     let strTitle = "";
     let eOption = null;
     const eSelect = document.createElement("select");

     for(var nCond=0; nCond<nCondCount; ++nCond)
     {
         strTitle = this.m_arConditions[nCond].toString();
         eOption = document.createElement("option");
         eOption.value = strTitle;
         eOption.text = strTitle;
         eSelect.appendChild(eOption);
     }

     eSelect.selectedIndex = this.m_nSelCondIndex;

     const filterThis = this;
     eSelect.addEventListener('change', (event) => {

         if(filterThis.m_bAPIIsCalling)
             return;

         const eSelectTmp =  event.currentTarget;
         const nCondTmp = eSelectTmp.selectedIndex;
         filterThis.setSelectedConditionIndex(nCondTmp);
         filterThis.requestFilter();
     });

     eDivCond.appendChild(eSelect);

     this.m_eSelectCond = eSelect;

     this.addUIItem(eDivCond);
     //let nHH = this.root.clientHeight;
     //this.root.style.height = "80px";
 }
}
