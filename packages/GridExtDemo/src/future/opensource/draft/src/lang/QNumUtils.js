
export class QNumUtils {
    constructor() {
        throw new Error("Never create instances of this class");
    }
}


let _qnumBuf = new DataView(new ArrayBuffer(8));
QNumUtils.getValue = function(x)
{
    _qnumBuf.setFloat64(0, x);
    let last =  _qnumBuf.getInt8(7) & 0xFC;
    _qnumBuf.setInt8(7, last);
    return  _qnumBuf.getFloat64(0);
};
