import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {SingleColumnFilter} from "./SingleColumnFilter";
import {TextUtils} from "../../utils/TextUtils";


export class MultilineTextFilter extends SingleColumnFilter
{
    constructor(ctxCol)
    {
        super(ctxCol);
        this.m_strText = "";
        this.m_arTextFrags = [];

        this.m_eTextArea = null;

        this.m_bHasFiltered = false;
    }

    isCompatible(ctxCol)
    {
        const b = ctxCol.getType() === DG.COLUMN_TYPE.STRING;
        return b;
    }


    getText()
    {
        return this.m_strText;
    }

    setText(strText)
    {
        if(strText === this.m_strText)
            return;

        this.m_strText = strText;
        const ar = TextUtils.splitTextIntoWords(strText);

        const nWordCount = ar.length;
        for(var nWord=0; nWord<nWordCount; ++nWord) {
            ar[nWord] = ar[nWord].toLowerCase();
        }

        this.m_arTextFrags = ar;
        this.m_eTextArea.value = strText;
    }

    isFiltered(nRecord)
    {
        const ctxCol = this.getColumnCtx();
        let ob = ctxCol.getValue(nRecord);
        if(ob === undefined || ob === null)
            ob = "";
        else ob = ob.toString().toLowerCase();

        const nWordCount = this.m_arTextFrags.length;
        let str = "";

        let b = nWordCount > 0;
        for(var nWord=0; nWord<nWordCount; ++nWord) {
            str =  this.m_arTextFrags[nWord];
            if (ob.includes(str)) {
                b = false;
                break;
            }
        }

        if(b)
            this.m_bHasFiltered = true;

        return b;
    }

    onBeforeApplyFilter()
    {
        this.m_bHasFiltered = false;
    }

    hasFiltered()
    {
      return this.m_bHasFiltered;
    }


    resetFilter(source)
    {
        const strTextOld = this.getText();
        this.setText("");
        super.resetFilter(source);

        //just keep the old text only in the edit box
        this.m_eTextArea.value = strTextOld;
    }

    render()
    {
        super.render();

        this.m_eTextArea = document.createElement("textarea");
        this.m_eTextArea.setAttribute('maxlength', 5000);
        this.m_eTextArea.setAttribute('cols',30);
        this.m_eTextArea.setAttribute('rows', 7);
        this.m_eTextArea.style.resize = "none";
        this.m_eTextArea.style.width = "100%";
        this.m_eTextArea.style.padding = "0px";

        const filterThis = this;
        const eButton = ui.button("Search", function(e){
            const strText = filterThis.m_eTextArea.value;
            filterThis.setText(strText);
            filterThis.requestFilter();
        });

        eButton.style.cssFloat = "right";

        this.addUIItem(this.m_eTextArea);
        this.addUIItem(eButton);

        //this.root.style.height = "170px";
    }
}
