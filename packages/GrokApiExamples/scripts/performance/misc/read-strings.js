// read 2 millions 50-byte strings from 100MB bytes array

let t1 = new Date().getTime();
let bytesView = new Uint8Array(100000000);
for (let i = 0; i < bytesView.length; i++)
    bytesView[i] = 65 + i % 20;

let t2 = new Date().getTime();
let strings = [];
let decoder = new TextDecoder();
for (let i = 0; i < bytesView.length / 50; i++) {
    let view = new Uint8Array(bytesView.buffer, i * 50, 50);
    strings[i] = decoder.decode(view);
}

console.log(`init: ${t2 - t1}`);
console.log(`read: ${new Date().getTime() - t2}`);