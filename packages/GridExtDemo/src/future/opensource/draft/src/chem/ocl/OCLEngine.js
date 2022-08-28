import {ChemEngine} from "../ChemEngine.js";
import Worker from "worker-loader!./OCLWorkerFile.js";

class _OCLEngine extends ChemEngine {
    getName() {
        return "OpenChemLib";
    }

    createWorker() {
        return new Worker();
    }

    getWorkerClass() {
        return Worker;
    }
}

export const OCLEngine = new _OCLEngine();