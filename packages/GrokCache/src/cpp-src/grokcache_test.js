let grokCache = {};


async function grokCacheModuleInit() {
  const htmlMStatus = document.querySelector('#gk_module_status');

  if (window.createGrokCache) {
    grokCache = await createGrokCache();
  } else {
    let err = "GrokCache WASM Module is not available!";
    htmlMStatus.textContent = err;
    console.error(err);
  }
}

// Test entry point
//
function grokCacheTestStart() {
  grokCacheModuleInit().then(() => {
    {
      const htmlMStatus = document.querySelector('#gk_module_status');
      let greeting = grokCache.wasm_version();
      htmlMStatus.textContent = greeting;
    }
    {
      const htmlEStatus = document.querySelector('#gk_enums_status');
      if (grokCache.CompressMethod) {
        console.log(grokCache.CompressMethod);

        if (grokCache.CompressMethod.NONE.value == 0 && grokCache.CompressMethod.ZIP.value == 1)
          htmlEStatus.textContent = "Ok";
        else
          htmlEStatus.textContent = "CompressMethod enum failure";
      }
    }
    {
      let cr = grokCache.getCacheRecord();
      console.log(cr);
    }

    // try to create byte buffer store
    //
    {
      let bufferStore = grokCache.create_ByteBufferStore();
      bufferStore.reset();
      if (bufferStore != null) {
        let buf_ptr = bufferStore.add_buffer(4);
        let storeSize = bufferStore.size();
        if (buf_ptr && storeSize == 1) {
          const str = "ATGC";
          let size = 4;
          let buf_arr = new Uint8ClampedArray(grokCache.HEAPF32.buffer, buf_ptr, size);
          buf_arr[0] = 65;
          buf_arr[1] = 84;
          buf_arr[2] = 71;
          buf_arr[3] = 67;
          var s2 = new String();

          // read it back from the byte buffer cache
          //
          let buf_ptr2 = bufferStore.get_ptr(0);
          let sz2 =  bufferStore.get_size(0);
          if (sz2 != size) {
            console.error(sz2);
            throw "Buffer size check failed (1)";
          }
          let buf_arr2 = new Uint8ClampedArray(grokCache.HEAPF32.buffer, buf_ptr2, sz2);
          for (var i = 0; i < size; i++) {
            var ch = buf_arr2[i];
            s2 += String.fromCharCode(ch);
          }
          console.log(s2);
          if (str == s2) {
            const htmlEStatus = document.querySelector('#gk_byte_buffer');
            htmlEStatus.textContent = "Byte Buffer OK";
          }
        }
      }
      delete bufferStore;
    }

    // draw canvas
    {
      const canvas = document.getElementById('canvas1');
      const ctx = canvas.getContext('2d');
      ctx.lineWidth = 10;

      ctx.strokeRect(75, 140, 150, 110);
      ctx.fillRect(130, 190, 40, 60);

      ctx.beginPath();
      ctx.moveTo(50, 140);
      ctx.lineTo(150, 60);
      ctx.lineTo(250, 140);
      ctx.closePath();
      ctx.stroke();

      // capture canvas content
      //window.location = canvas.toDataURL("image/png");

      var sizeWidth = ctx.canvas.clientWidth;
      var sizeHeight = ctx.canvas.clientHeight;

      let imageData = ctx.getImageData(0, 0, sizeWidth, sizeHeight);

      var buf1 = imageData.data;
      console.log(buf1.length);



      {
        let bufferStore = grokCache.create_ByteBufferStore();
        {
          let buf_ptr = bufferStore.add_buffer(buf1.length);
          let buf_arr = new Uint8ClampedArray(grokCache.HEAPF32.buffer, buf_ptr, imageData.data.length);
          buf_arr.set(imageData.data);
        }

        const canvas2 = document.getElementById('canvas2');
        const ctx2 = canvas2.getContext('2d');
        let imageData2 = ctx2.getImageData(0, 0, sizeWidth, sizeHeight);

        {
          let buf_ptr2 = bufferStore.get_ptr(0);
          let sz2 = bufferStore.get_size(0);
          if (sz2 != imageData.data.length) {
            console.error(sz2);
            throw "Buffer size check failed (1)";
          }
          let buf_arr2 = new Uint8ClampedArray(grokCache.HEAPF32.buffer, buf_ptr2, sz2);
          imageData2.data.set(buf_arr2); // memory copy the image
        }

        //imageData2.data.set(buf1); // memory copy the image
        ctx2.putImageData(imageData2, 0, 0);

        delete bufferStore;
      }

    }




  });
}