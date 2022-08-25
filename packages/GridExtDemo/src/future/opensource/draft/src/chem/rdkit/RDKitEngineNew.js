import {ChemEngine} from "../ChemEngine.js";
import Worker from "worker-loader!./RDKitWorkerFileNew.js";

class _RDKitEngine extends ChemEngine {
    getName() {
        return "RDKit";
    }

    createWorker() {
        return new Worker();
    }
    getWorkerClass() {
        return Worker;
    }

}

export const RDKitEngineNew = new _RDKitEngine();
