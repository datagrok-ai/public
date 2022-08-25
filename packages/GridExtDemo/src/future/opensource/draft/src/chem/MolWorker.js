import {ExecutorWorker} from "../concurrent/ExecutorWorker";
import {InitSerializer} from "./messaging/InitSerializer";
import {InitResponseMessage} from "./messaging/InitResponseMessage";
import {ImageSerializer} from "./messaging/ImageSerializer";
import {ImageResponseMessage} from "./messaging/ImageResponseMessage";
import {TanimotoSerializer} from "./messaging/TanimotoSerializer";
import {TanimotoResponseMessage} from "./messaging/TanimotoResponseMessage";
import {StructSearchSerializer} from "./messaging/StructSearchSerializer";
import {StructSearchResponseMessage} from "./messaging/StructSearchResponseMessage";
import {ProgressSerializer} from "../concurrent/messaging/ProgressSerializer";
import {ProgressRequestMessage} from "../concurrent/messaging/ProgressRequestMessage";
import {DisposeSerializer} from "../concurrent/messaging/DisposeSerializer";
import {DisposeResponseMessage} from "../concurrent/messaging/DisposeResponseMessage";
import {NotImplementedError} from "../lang/NotImplementedError";

export class MolWorker extends ExecutorWorker
{
    constructor(fnDispose)
    {
        super(fnDispose);
        this.m_nMolCount = -1;
        this.m_bMolFile = false;
        this.m_arTmpSmiles = null;
        this.m_nIdxMol = 0;
    }

   async init(strWebRoot, arSmiles, bMolFile) { throw new NotImplementedError();}
   continueInit(arSmiles, nIdxSmiles, bMolFile) { throw new NotImplementedError();}
   tanimotoScores(strSmilesFrag) {throw new NotImplementedError();}
   structureSearch(strSmilesFrag) {throw new NotImplementedError();}
   paintStructure(ofscreen, nMol, nX, nY, nW, nH) {throw new NotImplementedError();}

   async onMessage(scope, worker, event)
   {
        const strMsg = event.data.msg;

        if(strMsg === ProgressSerializer.Instance.getName())
        {
            const msgRespPrg = ProgressSerializer.Instance.toResponseMessage(event.data);
            const bCancel = msgRespPrg.isCancel();
            if(bCancel) {
                worker.m_arTmpSmiles = null;
                const msgRes = new InitResponseMessage(worker.m_nIdxMol, true);
                const argsRes = InitSerializer.Instance.toResponseArgs(msgRes);
                scope.postMessage(argsRes);
                return;
            }

            const nMolCount = worker.m_nMolCount;

            const nProcMolCount = worker.continueInit(worker.m_arTmpSmiles, worker.m_nIdxMol, worker.m_bMolFile);
            worker.m_nIdxMol = nProcMolCount;

            if(nProcMolCount !== nMolCount)
            {
                const msgReq = new ProgressRequestMessage(nProcMolCount);
                const argsReq = ProgressSerializer.Instance.toRequestArgs(msgReq);
                scope.postMessage(argsReq);//{progress: nProcMolCount});
            }
            else //all done
            {
              worker.m_arTmpSmiles = null;
              const msgRes = new InitResponseMessage(nMolCount);
              const argsRes = InitSerializer.Instance.toResponseArgs(msgRes);
              scope.postMessage(argsRes);
            }
        }

        else if(strMsg === InitSerializer.Instance.getName())
        {
         const msgReqInit = InitSerializer.Instance.toRequestMessage(event.data);
         const nMolCount = msgReqInit.getSmilesOrMolFiles().length;
         worker.m_nMolCount = nMolCount;
         worker.m_arTmpSmiles = msgReqInit.getSmilesOrMolFiles();

         await worker.init(msgReqInit.getWebRoot(), msgReqInit.getSmilesOrMolFiles(), msgReqInit.isMolFile());

         const msgReqPrg = new ProgressRequestMessage(0);
         const argsReqPrg = ProgressSerializer.Instance.toRequestArgs(msgReqPrg);
         scope.postMessage(argsReqPrg);
        }

        else if(strMsg === DisposeSerializer.Instance.getName())
        {
          worker.dispose();

          const messageResponse = new DisposeResponseMessage();
          const argsResponse = DisposeSerializer.Instance.toResponseArgs(messageResponse);

          scope.postMessage(argsResponse);
        }

        else if(strMsg === TanimotoSerializer.Instance.getName())
        {
            const message = TanimotoSerializer.Instance.toRequestMessage(event.data);
            const strSmilesFrag = message.getFragmentSmiles();
            const arScores = worker.tanimotoScores(strSmilesFrag);

            const messageResponse = new TanimotoResponseMessage(arScores);
            const argsResponse = TanimotoSerializer.Instance.toResponseArgs(messageResponse);
            scope.postMessage(argsResponse);
        }

        else if(strMsg === StructSearchSerializer.Instance.getName())
        {
            const message = StructSearchSerializer.Instance.toRequestMessage(event.data);
            const strSmilesFrag = message.getFragmentSmiles();

            const arFlags = worker.structureSearch(strSmilesFrag);
            const messageResponse = new StructSearchResponseMessage(arFlags);
            const argsResponse = StructSearchSerializer.Instance.toResponseArgs(messageResponse);
            scope.postMessage(argsResponse);
        }


        else if(strMsg === ImageSerializer.Instance.getName())
        {
            const message = ImageSerializer.Instance.toRequestMessage(event.data);

            const nW = message.getW();
            const nH = message.getH();
            const nMol= message.getMolIndex();

            const offscreen = new OffscreenCanvas(nW, nH);
            worker.paintStructure(offscreen, nMol, 0,0, nW, nH);
            const bitmap = offscreen.transferToImageBitmap();

            /*
            for(var nM=0; nM<worker.m_nMolCount; ++nM)
            {
                worker.paintStructure(offscreen, nM, 0,0, nW, nH)
                bitmap = offscreen.transferToImageBitmap();
            } */


            const messageResponse = new ImageResponseMessage(bitmap);
            const argsResponse = ImageSerializer.Instance.toResponseArgs(messageResponse);
            scope.postMessage(argsResponse);
        }
        else throw new Error("Unknown message name: " + strMsg)
   }
}
