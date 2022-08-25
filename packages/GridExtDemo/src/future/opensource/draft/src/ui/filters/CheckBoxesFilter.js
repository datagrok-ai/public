import {DescreteItemsFilter} from "./DescreteItemsFilter";

export class CheckBoxesFilter extends DescreteItemsFilter
{
    constructor(ctxCol)
    {
        super(ctxCol);

        this.m_bAPIIsCalling = false;
    }


    isSelectedItem(nItem)
    {
        const bFiltered = this.isFilteredItem(nItem);
        return !bFiltered;
    }

    setSelectedItem(nItem, bSelected)
    {
     const bFilteredOld = this.isFilteredItem(nItem);
     if((bFilteredOld && !bSelected) || !bFilteredOld && bSelected)
         return;

     this.setFilteredItem(nItem, !bSelected, false);

     const eCB = this.getUIItem(nItem);
        this.m_bAPIIsCalling = true;
     eCB.checked = bSelected;
        this.m_bAPIIsCalling = false;

    }


    resetFilter(source)
    {
        const nItemCount = this.getItemCount();
        for(var nItem=0; nItem<nItemCount; ++nItem)
        {
         if(!this.isSelectedItem(nItem))
             this.setSelectedItem(nItem, true);
        }

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
        const dframe = ctxCol.getDataFrame();
        let idCB = "";

        const nItemCount = this.getItemCount();
        let nH = 0;
        let strTitle = "";

        const filterThis = this;

        for (var nItem=0; nItem<nItemCount; ++nItem) {

            strTitle = this.getItemTitle(nItem);

            eDiv = document.createElement("div");
            eDiv.classList.add('y-axis-layout');

            eCB = document.createElement("INPUT");
            idCB = dframe.name + "_" + ctxCol.getName() + "_" + CheckBoxesFilter.CheckInstanceCount.toString();
            ++CheckBoxesFilter.CheckInstanceCount;
            eCB.id = idCB;
            eCB.style.display = "inline-block";
            eCB.setAttribute("type", "checkbox");
            eCB.setAttribute("item_id", nItem.toString());

            eCB.checked = !this.isFilteredItem(nItem);

            eCB.addEventListener('change', (event) => {

                if(filterThis.m_bAPIIsCalling)
                    return;

                const strItem = event.currentTarget.getAttribute("item_id");
                const nItemTmp = parseInt(strItem);
                const bSelected = event.currentTarget.checked;
                filterThis.setSelectedItem(nItemTmp, bSelected);
                //filterThis.setFilteredItem(nItemTmp, !event.currentTarget.checked, true);
                filterThis.requestFilter();
            });

            this.setUIItem(nItem, eCB);

            eDiv.appendChild(eCB);

            eLbl = document.createElement("Label");
            eLbl.style.display = "inline-block";
            eLbl.setAttribute("for",idCB);
            eLbl.innerHTML = strTitle;
            eDiv.appendChild(eLbl);

            this.addUIItem(eDiv);

        }

    }
}

CheckBoxesFilter.CheckInstanceCount = 0;