import {ProgressWorker} from "./ProgressWorker";

(function () {

    let m_worker = new ProgressWorker(() => {m_worker = null;});
    self.onmessage = async function (event)
    {
        await m_worker.onMessage(this, m_worker, event);
    }
}()); // self-executing anonymous function