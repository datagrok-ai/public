import { DockNode } from "./DockNode";

export class Utils {

    private static _counter: number = 0;

    static getPixels(pixels: string): number {
        if (pixels === null) {
            return 0;
        }

        return parseInt(pixels.replace('px', ''));
    }

    static disableGlobalTextSelection(element: HTMLElement) {
        element.classList.add('disable-selection');
    }

    static enableGlobalTextSelection(element: HTMLElement) {
        element.classList.remove('disable-selection');
    }

    static isPointInsideNode(px: number, py: number, node: DockNode): boolean {
        let element = node.container.containerElement;
        let rect = element.getBoundingClientRect();

        return (
            px >= rect.left &&
            px <= rect.left + rect.width &&
            py >= rect.top &&
            py <= rect.top + rect.height
        );
    }

    static getNextId(prefix: string): string {
        return prefix + Utils._counter++;
    }

    static removeNode(node: Node): boolean {
        if (node.parentNode === null) {
            return false;
        }

        node.parentNode.removeChild(node);

        return true;
    }

    static orderByIndexes<T>(array: T[], indexes: number[]) {
        let sortedArray = [];
        for (let i = 0; i < indexes.length; i++) {
            sortedArray.push(array[indexes[i]]);
        }
        return sortedArray;
    }

    static arrayRemove<T>(array: T[], value: any): T[] | false {
        let idx = array.indexOf(value);
        if (idx !== -1) {
            return array.splice(idx, 1);
        }
        return false;
    }

    static arrayContains<T>(array: T[], value: T): boolean {
        let i = array.length;
        while (i--) {
            if (array[i] === value) {
                return true;
            }
        }
        return false;
    }

    static arrayEqual<T>(a: T[], b: T[]): boolean {
        if (a === b) return true;
        if (a == null || b == null) return false;
        if (a.length != b.length) return false;

        for (let i = 0; i < a.length; ++i) {
            if (a[i] !== b[i]) return false;
        }
        return true;
    }
}
