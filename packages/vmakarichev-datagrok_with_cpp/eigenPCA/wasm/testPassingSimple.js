console.log("Hi!");

let num = 129;

console.log(num >> 2);

var factory = require('./newLib');

factory().then((instance) => {

    let a = 15;
    let b = 19;

    let res = instance.ccall('sumOfInts', 
    'number', ['number', 'number'], [a, b]);

    console.log(a + ' + ' + b + ' = ' + res);

    let arr1 = [1, 2, 3, 4];
    let tArr1 = new Int32Array(arr1);
    console.log(tArr1);

    let arr2 = [5, 6, 7, -8];
    let tArr2 = new Int32Array(arr2);
    console.log(tArr2);

    let arr3 = [15, 26, 17, -10];
    let tArr3 = new Int32Array(arr3);
    console.log(tArr3);

    let buf = instance._malloc(3 * tArr1.length * tArr1.BYTES_PER_ELEMENT);

    console.log('buf is ' + buf);
    console.log('buf >> 2 is ' + (buf >> 2));

    instance.HEAP32.set(tArr1, buf >> 2);
    instance.HEAP32.set(tArr2, (buf + 4 * tArr1.BYTES_PER_ELEMENT) >> 2);
    instance.HEAP32.set(tArr3, (buf + 8 * tArr1.BYTES_PER_ELEMENT) >> 2);

    res = instance.ccall('minOfArray', 
    'number', ['number', 'number'], [buf, 3 * arr1.length]);

    for(let i = 0; i < 3 * arr1.length; i++){
        console.log('  ' + instance.HEAP32[buf / tArr1.BYTES_PER_ELEMENT + i]);
    }

    console.log('min is ' + res);

    instance._free(buf);
});
