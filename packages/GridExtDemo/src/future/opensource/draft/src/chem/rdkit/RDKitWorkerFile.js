import {RDKitWorker} from "./RDKitWorker";

(function () {

let m_worker = new RDKitWorker(() => {m_worker = null;});
self.onmessage = async function (event)
{
  await m_worker.onMessage(this, m_worker, event);
}
}()); // self-executing anonymous function