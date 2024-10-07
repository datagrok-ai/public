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
  
    var xColumnData = sc.table.col(sc.props.xColumnName)?.getRawData();
    var yColumnData = sc.table.col(sc.props.yColumnName)?.getRawData();
    if (!xColumnData || !yColumnData)
      return false;
    var viewport = sc.viewport;
    var inverseX = sc.props.invertXAxis;
    var inverseY = !sc.props.invertYAxis;
    var linearX = sc.props.xAxisType == "linear";
    var linearY = sc.props.yAxisType == "linear";
  
    const filteredArray = sc.filter.getSelectedIndexes();
    const pointsCount = filteredArray.length;
    const points = new Float32Array(pointsCount * 3);
    var pointsIndex = 0;
    for (let i = 0; i < pointsCount; i++) {
      var index = filteredArray[i];
      var sx = gpuViewbox.left + worldToScreen(xColumnData[index], gpuViewbox.width, viewport.left, viewport.right, inverseX, linearX);
      var sy = gpuViewbox.top + worldToScreen(yColumnData[index], gpuViewbox.height, viewport.top, viewport.bottom, inverseY, linearY);
      points[pointsIndex++] = sx;
      points[pointsIndex++] = sy;
      points[pointsIndex++] = sc.getMarkerSize(index); // size
    }
  
    webGPUInit(canvas, points, sc);
}
  
function worldToScreen(world: number, length: number, min: number, max: number, inverse: boolean, linear: boolean) {
if (linear)
    return worldToScreenLinear(world, length, min, max, inverse);
else
    return worldToScreenLogarithmic(world, length, min, max, inverse);
}

function worldToScreenLinear(world: number, length: number, min: number, max: number, inverse: boolean) {
return inverse
    ? ((1 - ((world - min)) / (max - min)) * length)
    : ((world - min) / (max - min) * length);
}

function worldToScreenLogarithmic(world: number, length: number, min: number, max: number, inverse: boolean) {
if (world <= 0) return NaN;

//eliminates problems with rounding
const minLogFloat = 1e-30;
if (Math.abs(world - min) < minLogFloat) return inverse ? length : 0.0;
if (Math.abs(world - max) < minLogFloat) return inverse ? 0.0 : length;

min = (min < minLogFloat ? minLogFloat : min);
world = (world < minLogFloat ? minLogFloat : world);

const res = (length * Math.log(world / min) / Math.log(max / min));
return inverse ? length - res : res;
}

function areEqual(rc1: DG.Rect, rc2: DG.Rect): boolean {
return Math.floor(rc1.left) == Math.floor(rc2.left) && Math.floor(rc1.top) == Math.floor(rc2.top)
    && Math.floor(rc1.width) == Math.floor(rc2.width) && Math.floor(rc1.height) == Math.floor(rc2.height);
}

async function webGPUInit(webGPUCanvas: HTMLCanvasElement, points: Float32Array, sc: DG.ScatterPlotViewer) {
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
        @location(0) position: vec2f,
        @location(1) size: f32,
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
        indexes: ${wgslArrayString(filteredArray)},
        xColumnData: ${wgslArrayString(xColumnData)},
        yColumnData: ${wgslArrayString(yColumnData)},
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
        let normalizedPos = convertPointToNormalizedCoords(vert.position.x, vert.position.y);
        vsOut.position = vec4f(normalizedPos + pos * vert.size / uni.resolution, 0, 1);
        vsOut.texcoord = pos * 0.5 + 0.5;
        return vsOut;
    }

    @fragment fn fs(vsOut: VSOutput) -> @location(0) vec4f {
        return textureSample(t, s, vsOut.texcoord);
    }

    fn convertPointToNormalizedCoords(x: f32, y: f32) -> vec2<f32> {
        let centerX = sc.viewBox.left + sc.viewBox.width / 2.0;
        let centerY = sc.viewBox.top + sc.viewBox.height / 2.0;

        // Map the x and y coordinates to the normalized coordinate system (-1 to 1)
        let nx = (x - centerX) / (sc.viewBox.width / 2.0);
        let ny = -((y - centerY) / (sc.viewBox.height / 2.0));

        return vec2<f32>(nx, ny);
    }

    fn worldToScreen(
    world: f32, 
    length: f32, 
    min: f32, 
    max: f32, 
    inverse: bool, 
    linear: bool
    ) -> f32 {
        if (linear) {
            return worldToScreenLinear(world, length, min, max, inverse);
        } else {
            return worldToScreenLogarithmic(world, length, min, max, inverse);
        }
    }

    fn worldToScreenLinear(
    world: f32, 
    length: f32, 
    min: f32, 
    max: f32, 
    inverse: bool
    ) -> f32 {
        if (inverse) {
            return (1.0 - ((world - min) / (max - min))) * length;
        } else {
            return ((world - min) / (max - min)) * length;
        }
    }

    fn worldToScreenLogarithmic(
    world: f32, 
    length: f32, 
    min: f32, 
    max: f32, 
    inverse: bool
    ) -> f32 {
        // Define a very small float value to avoid precision issues
        let minLogFloat: f32 = 1e-30;

        // Handle edge cases when world is <= 0 or near the boundary
        if (world <= 0.0) {
            return 0.0;
        }
        if (abs(world - min) < minLogFloat) {
            if (inverse) { 
                return length;
            } else { 
                return 0.0; 
            }
        }
        if (abs(world - max) < minLogFloat) {
            if (inverse) { 
                return 0.0; 
            } else {
                return length;
            }
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

        if (inverse) {
            return length - res;
        }
        return res;
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
        arrayStride: (2 + 1) * 4, // 3 floats, 4 bytes each
        stepMode: 'instance',
        attributes: [
            {shaderLocation: 0, offset: 0, format: 'float32x2'},  // position
            {shaderLocation: 1, offset: 8, format: 'float32'},  // size
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
    size: points.byteLength,
    usage: GPUBufferUsage.VERTEX,
    mappedAtCreation: true,
});
new Float32Array(vertexBuffer.getMappedRange())
    .set(new Float32Array(points.buffer));
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
const kFilteredArrayByteOffset = getPaddedSize(filteredArray.length);
const kXColumnDataSize = getPaddedSize(xColumnData.length);
const kYColumnDataSize = getPaddedSize(yColumnData.length);

const scBuffer = gpuDevice.createBuffer({
    size: kAlignmentOffset + kViewBoxByteOffset + kViewPortByteOffset + kPropsByteOffset + kFilteredArrayByteOffset + kXColumnDataSize + kYColumnDataSize,
    usage: GPUBufferUsage.STORAGE,
    mappedAtCreation: true,
});
var scBufferArray = scBuffer.getMappedRange() // get full range
  // copy the data into the buffer at correct places
let scBufferOffset = kAlignmentOffset;
const viewBoxPortArray = new Float32Array(scBufferArray, scBufferOffset, 4 + 4);
viewBoxPortArray.set([gpuViewbox.left, gpuViewbox.top, gpuViewbox.width, gpuViewbox.height,
                      viewport.left, viewport.top, viewport.width, viewport.height]);
scBufferOffset += kViewBoxByteOffset + kViewPortByteOffset;
const propsArray = new Uint32Array(scBufferArray, scBufferOffset, 4);
propsArray.set([+!!inverseX, +!!inverseY, +!!linearX, +!!linearY]);
scBufferOffset += kPropsByteOffset;
const filteredBufferArray = new Int32Array(scBufferArray, scBufferOffset, filteredArray.length);
filteredBufferArray.set(filteredArray);
scBufferOffset += kFilteredArrayByteOffset;
if (xColumnData instanceof Int32Array) {
    new Int32Array(scBufferArray, scBufferOffset, xColumnData.length).set(xColumnData);
} else if (xColumnData instanceof Uint32Array) {
    new Uint32Array(scBufferArray, scBufferOffset, xColumnData.length).set(xColumnData);
} else {
    new Float32Array(scBufferArray, scBufferOffset, xColumnData.length).set(xColumnData);
}
scBufferOffset += kXColumnDataSize;
if (yColumnData instanceof Int32Array) {
    new Int32Array(scBufferArray, scBufferOffset, yColumnData.length).set(yColumnData);
} else if (yColumnData instanceof Uint32Array) {
    new Uint32Array(scBufferArray, scBufferOffset, yColumnData.length).set(yColumnData);
} else {
    new Float32Array(scBufferArray, scBufferOffset, yColumnData.length).set(yColumnData);
}
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
pass.draw(6, points.length / 3);
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

function wgslArrayString(array: Int32Array | Float32Array | Float64Array | Uint32Array): string {
    return `array<${array instanceof Int32Array ? "i32" : array instanceof Uint32Array ? "u32" : "f32"}, ${array.length}>`;
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
