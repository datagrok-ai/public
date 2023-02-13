var factory = require('./sum.js');

factory().then((instance) => {
    let len = 5;
    let arr = new Int32Array(len);
    arr.set([1,2,3,4,5]);
    console.log(arr);
    
    console.log(instance._sum(arr.buffer, arr.length));
});