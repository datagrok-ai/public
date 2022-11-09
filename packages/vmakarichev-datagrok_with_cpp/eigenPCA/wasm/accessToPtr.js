var factory = require('./newLib');

factory().then((instance) => {

    let sarr1 = [11, -22, 23, 14];
    let arr1 = new Int32Array(sarr1);

    let sarr2 = [1,2,3,4];
    let arr2 = new Int32Array(sarr2);    

    let buf1 = /*arr1.buffer;*/instance._malloc(16);

    let buf2 = /*arr2.buffer;*/ instance._malloc(16);

    console.log(buf1);
    console.log(buf2);

    instance.HEAPF32.set(arr1, buf1 >> 2);
    instance.HEAPF32.set(arr2, buf2 >> 2);

    console.log(instance._error(buf1, buf2, 4));

    instance._free(buf1);
    instance._free(buf2);
});