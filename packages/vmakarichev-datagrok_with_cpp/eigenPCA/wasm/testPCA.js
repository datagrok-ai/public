var factory = require('./newLib');

function pca(module, sourceArrays, numOfPrincipalComponents, principalComponents, approxData ){
    
    /*console.log('\n\nSource data is:');
    for(let arr of sourceArrays){
        console.log(arr);
    }*/

    let height = sourceArrays.length;
    let width = sourceArrays[0].length;
    let sizeOfData = height * width;
    let sizeOfPrincComp = numOfPrincipalComponents * width; 
    let numOfBytes = sourceArrays[0].BYTES_PER_ELEMENT;

    let shift = 2;
    let heap = module.HEAPF32;

    let bufForData = module._malloc(sizeOfData * numOfBytes);

    for(let i = 0; i < height; i++)
        heap.set(sourceArrays[i], (bufForData + i * width * numOfBytes) >> shift);

    console.log('\nData:');
    for(let i = 0; i < height; i++)
        console.log(new Float32Array(heap.buffer, bufForData + i * width * numOfBytes, width));

    let bufForPrincComp = module._malloc(sizeOfPrincComp * numOfBytes);
    let bufForApprox = module._malloc(sizeOfData * numOfBytes);

    console.log(bufForData);
    console.log(bufForPrincComp);

    let res = module.ccall('principalComponentAnalysis', 
    'number', ['number', 'number', 'number', 'number', 'number', 'number'], 
    [bufForData, height, width, 
        numOfPrincipalComponents, bufForPrincComp, bufForApprox]);

    for(let i = 0; i < numOfPrincipalComponents; i++)
        principalComponents.push(new Float32Array(heap.buffer, 
            bufForPrincComp + i * width * numOfBytes, width));

    for(let i = 0; i < height; i++)
        approxData.push(new Float32Array(heap.buffer, 
            bufForApprox + i * width * numOfBytes, width));
    
    /*console.log('\nData again:');
    for(let i = 0; i < height; i++)
        console.log(new Float32Array(heap.buffer, bufForData + i * width * numOfBytes, width));*/

    module._free(bufForData);
    module._free(bufForPrincComp);
    module._free(bufForApprox);
}

factory().then((instance) => {
    
    //let arr1 = [1, 2, 3, 4, 5];
    let arr1 = [11, -22, 23, 14, -50, 34, 21];
    let tArr1 = new Float32Array(arr1);
    console.log(tArr1);

    //let arr2 = [6, 7, 8, 9, 10];
    let arr2 = [26, 17, 18, 29, -10, -3, 32];
    let tArr2 = new Float32Array(arr2);
    console.log(tArr2);

    //let arr3 = [11, 12, 13, 14, 15];
    let arr3 = [11, -12, -13, -14, 15, 11, 32];
    let tArr3 = new Float32Array(arr3);
    console.log(tArr3);

    //let arr4 = [16, 17, 18, 19, 20];
    let arr4 = [-16, -17, 18, -19, -20, 14, 71];
    let tArr4 = new Float32Array(arr4);
    console.log(tArr4);

    /*let height = 4;
    let width = 5;
    let size = height * width;
    let numOfComponents = height;

    let numOfBytes = tArr1.BYTES_PER_ELEMENT;

    let buf1 = instance._malloc(size * numOfBytes);
 
    instance.HEAPF32.set(tArr1, buf1 >> 2);
    instance.HEAPF32.set(tArr2, (buf1 + width * numOfBytes) >> 2);
    instance.HEAPF32.set(tArr3, (buf1 + 2 * width * numOfBytes) >> 2);
    instance.HEAPF32.set(tArr4, (buf1 + 3 * width * numOfBytes) >> 2);       
    
    let data = new Float32Array(size);
    
    for(let i = 0; i < size; i++){
        data[i] = instance.HEAPF32[buf1 / numOfBytes + i];
    }

    console.log('Data is ' + data);

    let buf2 = instance._malloc(size * numOfBytes);
    let buf3 = instance._malloc(width * numOfComponents * numOfBytes);

    let approximation = new Float32Array(size);
    let principalComponents = new Float32Array(width * numOfComponents);

    console.log('\nApproximation (before): ' + approximation);
    console.log('\nPrincipal components (before): ' + principalComponents);

    let res = instance.ccall('principalComponentAnalysis', 
    'number', ['number', 'number', 'number', 'number', 'number', 'number'], 
    [buf1, height, width, numOfComponents, buf3, buf2]);

    for(let i = 0; i < size; i++){
        approximation[i] = instance.HEAPF32[buf2 / numOfBytes + i];
    }

    for(let i = 0; i < width * numOfComponents; i++){
        principalComponents[i] = instance.HEAPF32[buf3 / numOfBytes + i];
    }

    console.log('\nApproximation (after): ' + approximation);
    console.log('\nPrincipal components (after): ' + principalComponents);

    instance._free(buf1);
    instance._free(buf2);
    instance._free(buf3);*/

    let princComp = [];
    let approx = [];
    
    let data = [tArr1, tArr2, tArr3, tArr4];
    
    pca(instance, data, 3, princComp, approx);

    console.log('\nPrincipal Components:');
    console.log(princComp);

    console.log('\nData vs Approximation:');
    for(let i = 0; i < data.length; i++){
        console.log('\nNo. ' + i);
        console.log('  source: ' + data[i]);
        console.log('  approximation: ' + approx[i]);
    }        
});
