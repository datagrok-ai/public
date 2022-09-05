import {SingleColumnFilter} from "./SingleColumnFilter";
import {MathUtils} from "../../utils/MathUtils";

export class DescreteItemsFilter extends SingleColumnFilter
{
    constructor(ctxCol)
    {
        super(ctxCol);
        this.m_arItems = this.createItems(ctxCol);
    }

    isCompatible(ctxCol)
    {
        const arCategories = ctxCol.getCategories();

        const b = arCategories.length < 10;
        return b;
    }


    isFiltered(nRecord)
    {
        let arIdxs = null;
        const arItems = this.m_arItems;
        const nItemCount = this.getItemCount();
        for(var nItem=0; nItem<nItemCount; ++nItem)
        {
            if(!arItems[nItem].isFiltered())
                continue;

            arIdxs = arItems[nItem].getIndices();
            for(var n=0; n<arIdxs.length; ++n)
            {
                if(arIdxs[n] === nRecord)
                    return true;
            }
        }

        return false;
    }

    hasFiltered()
    {
        const arItems = this.m_arItems;
        const nItemCount = this.getItemCount();
        for(var nItem=0; nItem<nItemCount; ++nItem)
        {
            if(arItems[nItem].isFiltered())
                return true;
        }

        return false;
    }

    createItems(ctxCol)
    {
        const arCategories = ctxCol.getCategories();
        const arItems = [];

        let item = null;
        let strTitle = "";
        let arIdxs = null;
        let obVal = null;
        const nRecordCount = ctxCol.getLength();

        let bHasNull = false;
        let nC = 0;
        for (;  nC<arCategories.length; ++nC)
        {
            obVal = arCategories[nC];
            if(MathUtils.isNullValue(obVal))
            {
                bHasNull = true;
                break;
            }
        }

        if(bHasNull)
        {
            let obValTmp = null;
            arIdxs = [];
            for (nC=0;  nC<arCategories.length; ++nC)
            {
                obVal = arCategories[nC];
                if(!MathUtils.isNullValue(obVal))
                 continue;

                for(var nR = 0; nR<nRecordCount; ++nR) {

                    obValTmp = ctxCol.getValue(nR);
                    if(obVal === obValTmp || (MathUtils.isNullValue(obVal) && MathUtils.isNullValue(obValTmp)))
                      arIdxs.push(nR);
                }
            }

            strTitle = "NULL";
            item = new DescreteFlterItem(strTitle, arIdxs);
            arItems.push(item);
        }



        for (nC=0;  nC<arCategories.length; ++nC) {
                obVal = arCategories[nC];

                if(MathUtils.isNullValue(obVal))
                continue;

                arIdxs = [];
                for(var nR = 0; nR<nRecordCount; ++nR) {

                    if(obVal === ctxCol.getValue(nR))
                      arIdxs.push(nR);
                }
            strTitle = obVal === true ? "Yes" : obVal === false ? "No" : obVal.toString();
            item = new DescreteFlterItem(strTitle, arIdxs);
            arItems.push(item);

        }

        return arItems;
    }


    getItemCount()
    {
        return this.m_arItems.length;
    }

    getItemTitle(nItem)
    {
        const item = this.m_arItems[nItem];
        return item.getTitle();
    }

    isFilteredItem(nItem)
    {
        const item = this.m_arItems[nItem];
        return item.isFiltered();
    }

    setFilteredItem(nItem, bFiltered, bNotify)
    {
        const item = this.m_arItems[nItem];
        item.setFiltered(bFiltered);

        if(bNotify !== undefined && bNotify)
            this.requestFilter();
    }


    getUIItem(nItem)
    {
        const item = this.m_arItems[nItem];
        return item.getUIItem();
    }

    setUIItem(nItem, eItem)
    {
        const item = this.m_arItems[nItem];
        return item.setUIItem(eItem);
    }


}

export class DescreteFlterItem
{
    constructor(strTitle, arRecordIdxs)
    {
        this.m_strTitle = strTitle;
        this.arRecordIdxs = arRecordIdxs;
        this.m_bFiltered = false;
        this.m_etem = null;
    }

    getTitle() {return this.m_strTitle;}
    getIndices() {return this.arRecordIdxs;}

    isFiltered()
    {
        return this.m_bFiltered;
    }

    setFiltered(bFiltered)
    {
        this.m_bFiltered = bFiltered;
    }

    getUIItem() {return this.m_etem;}

    setUIItem(eItem)
    {
        this.m_etem = eItem;
    }
}