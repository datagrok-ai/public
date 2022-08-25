import Worker from "worker-loader!./ProgressWorkerFile.js";
import {ExecutorService} from "../../concurrent/ExecutorService";
import {InitRequestMessage} from "./messaging/InitRequestMessage";
import {InitSerializer} from "./messaging/InitSerializer";
import {ProgressPctRequestMessage} from "./messaging/ProgressPctRequestMessage";
import {ProgressPctSerializer} from "./messaging/ProgressPctSerializer";
import {TextRequestMessage} from "./messaging/TextRequestMessage";
import {TextSerializer} from "./messaging/TextSerializer";

export class ProgressIndicator
{
  constructor(bSkipBusy)
  {
   this.m_service = new ExecutorService(Worker, bSkipBusy);
  }

  async dispose() {
    await this.m_service.dispose();
    this.m_service = null;
  }

  async init(eCanvas)
  {
    await this.m_service.init(1);

    function CreatePostArgs(nCPU, nCPUCount, nJobFr, nJobTo)
    {
      const msg = new InitRequestMessage(eCanvas);
      const arArgs = InitSerializer.Instance.toRequestArgs(msg);

      return arArgs;
    }

    function OnResult(service, nCPU, nFr, nTo, obData)
    {
      const args = obData;
      const message = InitSerializer.Instance.toResponseMessage(args);
    }

    function OnError(service, nCPU, nFr, nTo, error)
    {

    }

    await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError);
  }

  async setProgress(nPct)
  {
    function CreatePostArgs(nCPU, nCPUCount)
    {
      const msg = new ProgressPctRequestMessage(nPct);
      const arArgs = ProgressPctSerializer.Instance.toRequestArgs(msg);

      return arArgs;
    }

    function OnResult(service, nCPU, nFr, nTo, obData)
    {
      let args = obData;
      let msg = ProgressPctSerializer.Instance.toResponseMessage(args);
    }

    function OnError(service, nCPU, nFr, nTo, error)
    {
      let CCC = 0;
    }

    await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError);
  }


  async setText(strText, bUpper)
  {
   if(bUpper === undefined)
      bUpper = false;

    function CreatePostArgs(nCPU, nCPUCount)
    {
      const msg = new TextRequestMessage(strText, bUpper);
      const arArgs = TextSerializer.Instance.toRequestArgs(msg);

      return arArgs;
    }

    function OnResult(service, nCPU, nFr, nTo, obData)
    {
      const args = obData;
      const msg = TextSerializer.Instance.toResponseMessage(args);


    }

    function OnError(service, nCPU, nFr, nTo, error)
    {
      let CCC = 0;
    }

    await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError);
  }
}

ProgressIndicator.create = async function(eCanvas, bSkipBusy)
{
  let indicator = new ProgressIndicator(bSkipBusy);

  await indicator.init(eCanvas);
  return indicator;
}