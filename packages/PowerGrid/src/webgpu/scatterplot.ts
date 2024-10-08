import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let gpuViewbox: DG.Rect;

export function scWebGPURender(sc: DG.ScatterPlotViewer, show: boolean) {
    // Getting WebGPU canvas or creating it if absent
    const canvasName = "webGPUCanvas";
    var canvas = sc.canvas.parentElement?.children.namedItem(canvasName) as HTMLCanvasElement;
    if (!show) {
      if (canvas)
        canvas.hidden = true;
      return;
    }
    if (!canvas) {
      canvas = ui.canvas(sc.viewBox.width ?? 100, sc.viewBox.height ?? 100);
      canvas.style.position = 'absolute';
      canvas.style.pointerEvents = 'none';
      canvas.id = canvasName;
      sc.canvas.parentElement?.insertBefore(canvas, sc.overlay);
    }
    canvas.hidden = false;
  
    if (!gpuViewbox || !areEqual(sc.viewBox, gpuViewbox)) {
      gpuViewbox = sc.viewBox;
      canvas.width = sc.viewBox.width;
      canvas.height = sc.viewBox.height;
      canvas.style.width = `${sc.viewBox.width}px`;
      canvas.style.height = `${sc.viewBox.height}px`;
      canvas.style.left = `${sc.viewBox.left}px`;
      canvas.style.top = `${sc.viewBox.top}px`;
    }
  
    webGPUInit(canvas, sc);
}

function areEqual(rc1: DG.Rect, rc2: DG.Rect): boolean {
    return Math.floor(rc1.left) == Math.floor(rc2.left) && Math.floor(rc1.top) == Math.floor(rc2.top)
        && Math.floor(rc1.width) == Math.floor(rc2.width) && Math.floor(rc1.height) == Math.floor(rc2.height);
}

async function webGPUInit(webGPUCanvas: HTMLCanvasElement, sc: DG.ScatterPlotViewer) {
    const filteredArray = sc.filter.getSelectedIndexes();
    var xColumnData = sc.table.col(sc.props.xColumnName)?.getRawData();
    var yColumnData = sc.table.col(sc.props.yColumnName)?.getRawData();
    if (!filteredArray.length || !xColumnData?.length || !yColumnData?.length)
        return;

    const adapter = await navigator.gpu?.requestAdapter();
    const gpuDevice = await adapter?.requestDevice() ?? null;  
    if (!gpuDevice) {
        fail('need a browser that supports WebGPU');
        return;
    }

    // Get a WebGPU context from the canvas and configure it
    const gpuContext = webGPUCanvas.getContext('webgpu');
    if (!gpuContext) {
        //this.fail('failed to get WebGPU context');
        return;
    }

    const presentationFormat = navigator.gpu.getPreferredCanvasFormat();
    gpuContext.configure({
        device: gpuDevice,
        format: presentationFormat,
        alphaMode: 'premultiplied',
    });

    const module = gpuDevice.createShaderModule({
        code: `
        struct Vertex {
            @location(0) index: i32,
        };

        struct Uniforms {
            resolution: vec2f,
        };

        struct VSOutput {
            @builtin(position) position: vec4f,
            @location(0) texcoord: vec2f,
        };

        struct Rect {
            left: f32,
            top: f32,
            width: f32,
            height: f32,
        };

        struct Props {
            inverseX: u32,
            inverseY: u32,
            linearX: u32,
            linearY: u32,
        };

        struct SC {
            alignment: vec4<u32>,
            viewBox: Rect,
            viewport: Rect,
            props: Props,
            xColumnData: array<f32, ${xColumnData.length}>,
            yColumnData: array<f32, ${yColumnData.length}>,

        };

        @group(0) @binding(0) var<uniform> uni: Uniforms;
        @group(0) @binding(1) var s: sampler;
        @group(0) @binding(2) var t: texture_2d<f32>;
        @group(0) @binding(3) var<storage, read> sc: SC;

        @vertex fn vs(
            vert: Vertex,
            @builtin(vertex_index) vNdx: u32,
        ) -> VSOutput {
            let points = array(
            vec2f(-1, -1),
            vec2f( 1, -1),
            vec2f(-1,  1),
            vec2f(-1,  1),
            vec2f( 1, -1),
            vec2f( 1,  1),
            );
            var vsOut: VSOutput;
            let pos = points[vNdx];

            let screenPoint = pointToScreen(vert.index);
            let normalizedPos = convertPointToNormalizedCoords(screenPoint);
            vsOut.position = vec4f(normalizedPos + pos * 10 / uni.resolution, 0, 1);
            vsOut.texcoord = pos * 0.5 + 0.5;
            return vsOut;
        }

        @fragment fn fs(vsOut: VSOutput) -> @location(0) vec4f {
            return textureSample(t, s, vsOut.texcoord);
        }

        fn pointToScreen(index: i32) -> vec2<f32> {
            var sx = sc.viewBox.left + worldToScreen(sc.xColumnData[index], sc.viewBox.width, sc.viewport.left, sc.viewport.left + sc.viewport.width, sc.props.inverseX, sc.props.linearX);
            var sy = sc.viewBox.top + worldToScreen(sc.yColumnData[index], sc.viewBox.height, sc.viewport.top, sc.viewport.top + sc.viewport.height, sc.props.inverseY, sc.props.linearY);
            return vec2<f32>(sx, sy);
        }

        fn convertPointToNormalizedCoords(point: vec2<f32>) -> vec2<f32> {
            let centerX = sc.viewBox.left + sc.viewBox.width / 2.0;
            let centerY = sc.viewBox.top + sc.viewBox.height / 2.0;

            // Map the x and y coordinates to the normalized coordinate system (-1 to 1)
            let nx = (point.x - centerX) / (sc.viewBox.width / 2.0);
            let ny = -((point.y - centerY) / (sc.viewBox.height / 2.0));

            return vec2<f32>(nx, ny);
        }

        fn worldToScreen(world: f32, length: f32, min: f32, max: f32, inverse: u32, linear: u32) -> f32 {
            return select(worldToScreenLogarithmic(world, length, min, max, inverse), worldToScreenLinear(world, length, min, max, inverse), linear != 0);
        }

        fn worldToScreenLinear(world: f32, length: f32, min: f32, max: f32, inverse: u32) -> f32 {
            return select(((world - min) / (max - min)) * length, (1.0 - ((world - min) / (max - min))) * length, inverse != 0);
        }

        fn worldToScreenLogarithmic(world: f32, length: f32, min: f32, max: f32, inverse: u32) -> f32 {
            // Define a very small float value to avoid precision issues
            let minLogFloat: f32 = 1e-30;

            // Handle edge cases when world is <= 0 or near the boundary
            if (world <= 0.0) {
                return 0.0;
            }
            if (abs(world - min) < minLogFloat) {
                return select(0.0, length, inverse != 0);
            }
            if (abs(world - max) < minLogFloat) {
                return select(length, 0.0, inverse != 0);
            }

            // Ensure minimum and world values are clamped to minLogFloat
            var adjustedMin = min;
            if (min < minLogFloat) { 
                adjustedMin = minLogFloat;
            }
            var adjustedWorld = world;
            if (world < minLogFloat) { 
                adjustedWorld = minLogFloat;
            }

            // Perform the logarithmic transformation
            let res = (length * log(adjustedWorld / adjustedMin) / log(max / adjustedMin));

            return select(res, length - res, inverse != 0);
        }
        `,
    });

    const pipeline = gpuDevice.createRenderPipeline({
        label: 'sizeable points with texture',
        layout: 'auto',
        vertex: {
        module,
        buffers: [
            {
            arrayStride: 4, // 1 int, 4 bytes
            stepMode: 'instance',
            attributes: [
                {shaderLocation: 0, offset: 0, format: 'sint32'},  // position
            ],
            },
        ],
        },
        fragment: {
        module,
        targets: [
            {
            format: presentationFormat,
            blend: {
                color: {
                srcFactor: 'one',
                dstFactor: 'one-minus-src-alpha',
                operation: 'add',
                },
                alpha: {
                srcFactor: 'one',
                dstFactor: 'one-minus-src-alpha',
                operation: 'add',
                },
            },
            },
        ],
        },
    });

    // Defining textures that will be drawn
    const textureSize = 32;
    const circleCanvas = createCircleCanvas(textureSize, 'blue', '#000033');
    const texture = gpuDevice.createTexture({
        size: [textureSize, textureSize],
        format: 'rgba8unorm',
        usage: GPUTextureUsage.TEXTURE_BINDING |
            GPUTextureUsage.COPY_DST |
            GPUTextureUsage.RENDER_ATTACHMENT,
    });
    gpuDevice.queue.copyExternalImageToTexture(
        { source: circleCanvas, flipY: true },
        { texture: texture, premultipliedAlpha: true },
        [textureSize, textureSize],
    );

    const sampler = gpuDevice.createSampler({
        minFilter: 'linear',
        magFilter: 'linear',
    });

    const vertexBuffer = gpuDevice.createBuffer({
        label: 'vertex buffer vertices',
        size: filteredArray.byteLength,
        usage: GPUBufferUsage.VERTEX,
        mappedAtCreation: true,
    });
    new Float32Array(vertexBuffer.getMappedRange())
        .set(new Float32Array(filteredArray.buffer));
    vertexBuffer.unmap();

    const uniformValues = new Float32Array(2);
    const uniformBuffer = gpuDevice.createBuffer({
        size: uniformValues.byteLength,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    const kResolutionOffset = 0;
    const resolutionValue = uniformValues.subarray(
        kResolutionOffset, kResolutionOffset + 2);

    var viewport = sc.viewport;
    var inverseX = sc.props.invertXAxis;
    var inverseY = !sc.props.invertYAxis;
    var linearX = sc.props.xAxisType == "linear";
    var linearY = sc.props.yAxisType == "linear";

    const kAlignmentOffset = 16;
    const kViewBoxByteOffset = 16;
    const kViewPortByteOffset = 16;
    const kPropsByteOffset = 16;
    //const kFilteredArrayByteOffset = getPaddedSize(filteredArray.length);
    const kXColumnDataSize = getPaddedSize(xColumnData.length);
    const kYColumnDataSize = getPaddedSize(yColumnData.length);

    const scBuffer = gpuDevice.createBuffer({
        size: kAlignmentOffset + kViewBoxByteOffset + kViewPortByteOffset + kPropsByteOffset + kXColumnDataSize + kYColumnDataSize,
        usage: GPUBufferUsage.STORAGE,
        mappedAtCreation: true,
    });
    var scBufferArray = scBuffer.getMappedRange() // get full range
    // copy the data into the buffer at correct places
    let scBufferOffset = kAlignmentOffset;
    new Float32Array(scBufferArray, scBufferOffset, 4 + 4)
        .set([gpuViewbox.left, gpuViewbox.top, gpuViewbox.width, gpuViewbox.height,
                        viewport.left, viewport.top, viewport.width, viewport.height]);
    scBufferOffset += kViewBoxByteOffset + kViewPortByteOffset;
    new Uint32Array(scBufferArray, scBufferOffset, 4)
        .set([+!!inverseX, +!!inverseY, +!!linearX, +!!linearY]);
    scBufferOffset += kPropsByteOffset;
    new Float32Array(scBufferArray, scBufferOffset, xColumnData.length).set(xColumnData);
    scBufferOffset += kXColumnDataSize;
    new Float32Array(scBufferArray, scBufferOffset, yColumnData.length).set(yColumnData);
    scBufferOffset += kYColumnDataSize;
    scBuffer.unmap();

    const bindGroup = gpuDevice.createBindGroup({
        layout: pipeline.getBindGroupLayout(0),
        entries: [
        { binding: 0, resource: { buffer: uniformBuffer }},
        { binding: 1, resource: sampler },
        { binding: 2, resource: texture.createView() },
        { binding: 3, resource: { buffer: scBuffer }},
        ],
    });

    const renderPassDescriptor = {
        label: 'our basic canvas renderPass',
        colorAttachments: [
        {
            // view: <- to be filled out when we render
            //clearValue: [0.9, 0.9, 0.9, 0],
            loadOp: 'clear',
            storeOp: 'store',
            view: gpuContext.getCurrentTexture().createView(),
        },
        ],
    };

    // Get the current texture from the canvas context and
    // set it as the texture to render to.
    const canvasTexture = gpuContext.getCurrentTexture();
    (renderPassDescriptor.colorAttachments as any)[0].view =
        canvasTexture.createView();

    // Update the resolution in the uniform buffer
    resolutionValue.set([canvasTexture.width, canvasTexture.height]);
    gpuDevice.queue.writeBuffer(uniformBuffer, 0, uniformValues);

    const encoder = gpuDevice.createCommandEncoder();
    const pass = encoder.beginRenderPass(renderPassDescriptor as GPURenderPassDescriptor);
    pass.setPipeline(pipeline);
    pass.setVertexBuffer(0, vertexBuffer);
    pass.setBindGroup(0, bindGroup);
    pass.draw(6, filteredArray.length);
    pass.end();

    const commandBuffer = encoder.finish();
    gpuDevice.queue.submit([commandBuffer]);

    // const observer = new ResizeObserver(entries => {
    //   for (const entry of entries) {
    //     const canvas = entry.target as C;
    //     const width = entry.contentBoxSize[0].inlineSize;
    //     const height = entry.contentBoxSize[0].blockSize;
    //     if (canvas) {
    //       canvas.width = Math.max(1, Math.min(width, device.limits.maxTextureDimension2D));
    //       canvas.height = Math.max(1, Math.min(height, device.limits.maxTextureDimension2D));
    //     }
    //     // re-render
    //     render();
    //   }
    // });
    // observer.observe(canvas);
}

function getPaddedSize(length: number): number {
  // The GPU struct must be padded to 16 bytes, so we need to calculate the size of the struct in 32bit values
  const bufferSize = length * Uint32Array.BYTES_PER_ELEMENT;
  let paddedComputeInfoBufferSize = bufferSize;
  const remainder = bufferSize & 15; // check if the size is a multiple of 16
  if (remainder !== 0)
    paddedComputeInfoBufferSize += 16 - remainder; // pad the size accordingly
  return paddedComputeInfoBufferSize;
}

function createCircleCanvas(size: number, fillColor: string, strokeColor: string): OffscreenCanvas {
    const lineWidth = 4;
    const canvas = new OffscreenCanvas(size, size);
    const ctx = canvas.getContext('2d');
    const centerX = size / 2;
    const centerY = size / 2;
    const radius = size / 2 - lineWidth;
    if (ctx) {
        ctx.beginPath();
        ctx.arc(centerX, centerY, radius, 0, 2 * Math.PI, false);
        ctx.fillStyle = fillColor;
        ctx.fill();
        ctx.lineWidth = lineWidth;
        ctx.strokeStyle = strokeColor;
        ctx.stroke();
    }

    return canvas;
}

function fail(msg: string) {
    alert(msg);
}