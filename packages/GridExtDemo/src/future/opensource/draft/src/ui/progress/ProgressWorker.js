import {ExecutorWorker} from "../../concurrent/ExecutorWorker";
import {InitSerializer} from "./messaging/InitSerializer";
import {InitResponseMessage} from "./messaging/InitResponseMessage";
import {DisposeResponseMessage} from "../../concurrent/messaging/DisposeResponseMessage";
import {DisposeSerializer} from "../../concurrent/messaging/DisposeSerializer";
import {ProgressPctSerializer} from "./messaging/ProgressPctSerializer";
import {ProgressPctResponseMessage} from "./messaging/ProgressPctResponseMessage";
import {TextSerializer} from "./messaging/TextSerializer";
import {TextResponseMessage} from "./messaging/TextResponseMessage";
import {TextUtils} from "../../utils/TextUtils";
import {DateUtils} from "../../utils/DateUtils";

export class ProgressWorker extends ExecutorWorker
{
    constructor(fnDispose) {
        super(fnDispose);

        this.m_nPercent = 0;
        this.m_strTextUp = "";
        this.m_strTextBottom = "";
        this.m_strFontText = "900 12px Arial";
        this.m_strFontPrg = "bold 16px Dialog";
        this.m_bStop = false;
    }


    dispose()
    {
        this.m_bStop = true;
        super.dispose();
    }

    async init(eCanvas)
    {
        const nW = eCanvas.width;
        const nH = eCanvas.height;
        const nHSect = nH/3;
        const nRadius = nHSect/2;
        const nThickness = 10;

        let fEndAngle = 0.0;
        let nRotate = 0;
        let worker = this;
        let timestampStart = undefined;
        function render(timestamp) {
                if(timestampStart === undefined)
                    timestampStart = timestamp;

                const elapsed = timestamp - timestampStart;

                const nMinCount = DateUtils.getMinCount(elapsed);
                const nSecCount = DateUtils.getSecCount(elapsed - nMinCount*60000);

                const g = eCanvas.getContext('2d');
                //g.save();
                //Clear the whole area
                g.fillStyle = "white";
                g.fillRect(0, 0, nW, nH);

                //Draw Upper Text
                g.fillStyle = "black";
                g.strokeStyle = "black";
                //g.strokeRect(0, 0, nW, nHSect);

                let str = worker.m_strTextUp;
                g.font = worker.m_strFontText;
                str = TextUtils.trimText(str, g, nW);

                let tm = g.measureText(str);
                let nWText  = tm.width;
                let nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent;
                let nYOffsetText = 2*nHFont;

                g.translate(eCanvas.width / 2, nHSect);
                g.fillText(str, -nWText/2, -tm.actualBoundingBoxDescent - nYOffsetText);
                g.translate(-eCanvas.width / 2, -nHSect);

                //Draw Bottom Text
               // g.strokeRect(0, 2*nHSect, nW, nHSect);

                str = worker.m_strTextBottom;
                g.font = worker.m_strFontText;
                str = TextUtils.trimText(str, g, nW);
                tm = g.measureText(str);
                nWText  = tm.width;
                nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent;
                nYOffsetText = 2*nHFont;

                g.translate(eCanvas.width / 2, 2*nHSect);
                g.fillText(str, -nWText/2, Math.abs(tm.actualBoundingBoxAscent) + nYOffsetText);
                g.translate(-eCanvas.width / 2, -2*nHSect);


                //Draw Time Elapsed
                str = nMinCount + ":" + (nSecCount < 10 ? "0" : "") + nSecCount;
                tm = g.measureText(str);
                nWText  = tm.width;
                nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent;
                g.translate(eCanvas.width / 2, 2*nHSect + 5*nHFont);
                g.fillText(str, -nWText/2, Math.abs(tm.actualBoundingBoxAscent) + nYOffsetText);
                g.translate(-eCanvas.width / 2, -2*nHSect -5*nHFont);

                //Progress
                g.font = worker.m_strFontPrg;
                str = worker.m_nPercent.toString() + "%";
                tm = g.measureText(str);
                nWText  = tm.width;
                nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent;

                g.translate(eCanvas.width / 2, eCanvas.height / 2);

                nRotate = (elapsed*Math.PI)/1000;
                g.rotate(nRotate);
                g.fillText(str, -nWText/2, nHFont/2);

                fEndAngle =  2*Math.PI * (worker.m_nPercent === 0 ? 1 : worker.m_nPercent/100);

                g.strokeStyle = worker.m_nPercent === 0 ? "LightPink" : "DodgerBlue";
                g.lineWidth = nThickness;
                g.beginPath();
                g.arc(0, 0, nRadius - (nThickness/2), 0, fEndAngle);//1.5 * Math.PI);
                g.stroke();
                g.beginPath();
                g.strokeStyle = worker.m_nPercent === 0 ? "LightPink" : "LightGray";
                g.arc(0, 0, nRadius - (nThickness/2), fEndAngle, 2.0*Math.PI);
                g.stroke();

                g.lineWidth = 1;

                g.rotate(-nRotate);
                g.translate(-eCanvas.width / 2, -eCanvas.height / 2);
            //g.restore();

            if(!worker.m_bStop)
             requestAnimationFrame(render);
            else
            {
                let ccc = 0;
            }
        }
        requestAnimationFrame(render);

    }

    setProgress(nPercent)
    {
        this.m_nPercent = nPercent;
    }

    setUpperText(strText)
    {
        this.m_strTextUp = strText;
    }

    setBottomText(strText)
    {
        this.m_strTextBottom = strText;
    }


    async onMessage(scope, worker, event)
    {
        let strMsg = event.data.msg;

     if(strMsg === InitSerializer.Instance.getName())
     {
        let msgReqInit = InitSerializer.Instance.toRequestMessage(event.data);
        let eCanvas = msgReqInit.getCanvas();

          worker.init(eCanvas);
        //await worker.init(msgReqInit.getWebRoot(), msgReqInit.getSmilesOrMolFiles(), msgReqInit.isMolFile());

         let msgRes = new InitResponseMessage();
         let argsRes = InitSerializer.Instance.toResponseArgs(msgRes);
         scope.postMessage(argsRes);
     }
     else if(strMsg === DisposeSerializer.Instance.getName())
     {
         worker.dispose();

         let messageResponse = new DisposeResponseMessage();
         let argsResponse = DisposeSerializer.Instance.toResponseArgs(messageResponse);

         scope.postMessage(argsResponse);
     }
     else if(strMsg === ProgressPctSerializer.Instance.getName())
     {
         let msgReq = ProgressPctSerializer.Instance.toRequestMessage(event.data);
         let nPct = msgReq.getPct();
         worker.setProgress(nPct);

         let msgResponse = new ProgressPctResponseMessage();
         let argsResponse = ProgressPctSerializer.Instance.toResponseArgs(msgResponse);

         scope.postMessage(argsResponse);
     }

     else if(strMsg === TextSerializer.Instance.getName())
     {
         let msgReq = TextSerializer.Instance.toRequestMessage(event.data);
         let strText = msgReq.getText();
         let bUp = msgReq.isUp();
         if(bUp)
             worker.setUpperText(strText);
         else
             worker.setBottomText(strText);

         let msgResponse = new TextResponseMessage();
         let argsResponse = TextSerializer.Instance.toResponseArgs(msgResponse);

         scope.postMessage(argsResponse);
     }

    }
}