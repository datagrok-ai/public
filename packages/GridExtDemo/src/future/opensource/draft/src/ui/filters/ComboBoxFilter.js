import {DescreteFlterItem, DescreteItemsFilter} from "./DescreteItemsFilter";

export class ComboBoxFilter extends DescreteItemsFilter
{
    constructor(ctxCol)
    {
        super(ctxCol);
        this.m_nSelItemIdx = 0;
        this.m_eSelect = null;
        this.m_bAPIIsCalling = false;
    }

    isCompatible(ctxCol)
    {
        const arCategories = ctxCol.getCategories();

        const b = arCategories.length < 50;
        return b;
    }

    createItems(ctxCol)
    {
        const arItems = super.createItems(ctxCol);

        const item = new DescreteFlterItem("All", []);
        arItems.unshift(item);

        return arItems;
    }

    getSelectedItemIndex()
    {
        return this.m_nSelItemIdx;
    }

    setSelectedItemIndex(nItem)
    {
        if(nItem === this.m_nSelItemId)
         return;

        this.m_nSelItemIdx = nItem;

        let b = false;
        const nItemCount = this.getItemCount();
        for (var nI=0; nI<nItemCount; ++nI)
        {
            b = nItem === 0 ? false : nI !== nItem;
            this.setFilteredItem(nI, b, false);
        }

        this.m_bAPIIsCalling = true;
        this.m_eSelect.selectedIndex = nItem;
        this.m_bAPIIsCalling = false;
    }


    resetFilter(source)
    {
     this.setSelectedItemIndex(0);
     super.resetFilter(source);
    }

    render()
    {
        //let name = `radio_${groupId++}`;
        super.render();

        let eDiv = null;
        let eCB = null;
        let eLbl = null;

        const ctxCol = this.getColumnCtx();
        //const dframe = ctxCol.getDataFrame();
        //let idCB = "";

        const nItemCount = this.getItemCount();
        let nH = 0;
        let strTitle = "";

        let eOption = null;
        const eSelect = document.createElement("select");
        const filterThis = this;
        for (var nItem=0; nItem<nItemCount; ++nItem) {

            strTitle = this.getItemTitle(nItem);
            eOption = document.createElement("option");
            eOption.value = strTitle;
            eOption.text = strTitle;
            eSelect.appendChild(eOption);
        }


        eSelect.addEventListener('change', (event) => {

            if(filterThis.m_bAPIIsCalling)
              return;

            const eSelectTmp =  event.currentTarget;
            const nItemTmp = eSelectTmp.selectedIndex;
            filterThis.setSelectedItemIndex(nItemTmp);
            filterThis.requestFilter();
        });

        this.m_eSelect = eSelect;
        this.addUIItem(eSelect);
        //let nHH = this.root.clientHeight;
        //this.root.style.height=(20 + 30).toString() + "px";
    }
}