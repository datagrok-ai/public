import {ProgressIndicator} from "./ProgressIndicator";

export class ProgressDialog {
    constructor(strTitle, bSkipBusy, fnStop)
    {
        this.m_bSkipBusy = bSkipBusy;
        let eCanvas = document.createElement("canvas");
        eCanvas.width  = 400;
        eCanvas.height = 400;
        this.m_eCanvas = eCanvas;
        this.m_progress = null;
        this.m_dial = ui.dialog(strTitle);
        this.m_dial.add(eCanvas);

        if(fnStop !== undefined) {
            const eButton = ui.bigButton("Cancel", fnStop);
            this.m_dial.add(eButton);
        }

        this.m_bCancelled = false;
    }

    isCancelled() {return this.m_bCancelled;}
    getCanvas() {return this.m_eCanvas;}

    async show()
    {
        this.m_progress = await ProgressIndicator.create(this.m_eCanvas, this.m_bSkipBusy);

        //this.m_dial.show({modal: true, x: (document.body.clientWidth - 300) / 2, y: (document.body.clientHeight - 300) / 2});
        this.m_dial.show({modal: true, center:true});
        let d = this;
        this.m_dial.onCancel(async function () {
            d.m_bCancelled = true;
            await d.m_progress.dispose();
            d.m_progress = null;
        });
    }


    async close()
    {
        if(this.m_progress === null)
            throw new Error("Dialog is not showing");

        await this.m_progress.dispose();
        this.m_progress = null;
        this.m_dial.close();
    }

    async setProgress(nPct)
    {
        if(this.m_progress === null)
            throw new Error("Dialog is not showing");

        await this.m_progress.setProgress(nPct);
    }

    async setUpperText(strText)
    {
        if(this.m_progress === null)
            throw new Error("Dialog is not showing");

        await this.m_progress.setText(strText, true);
    }

   async setBottomText(strText)
   {
      if(this.m_progress === null)
       throw new Error("Dialog is not showing");

       await this.m_progress.setText(strText, false);
   }
}
