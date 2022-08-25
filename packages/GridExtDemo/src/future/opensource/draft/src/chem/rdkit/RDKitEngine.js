import {ChemEngine} from "../ChemEngine.js";
import Worker from "worker-loader!./RDKitWorkerFile.js";

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

export const RDKitEngine = new _RDKitEngine();