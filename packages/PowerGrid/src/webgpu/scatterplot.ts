import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let gpuViewboxRect: DG.Rect;

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
  
    if (!gpuViewboxRect || !areEqual(sc.viewBox, gpuViewboxRect)) {
      gpuViewboxRect = sc.viewBox;
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
      var sx = gpuViewboxRect.left + worldToScreen(xColumnData[index], gpuViewboxRect.width, viewport.left, viewport.right, inverseX, linearX);
      var sy = gpuViewboxRect.top + worldToScreen(yColumnData[index], gpuViewboxRect.height, viewport.top, viewport.bottom, inverseY, linearY);
      var normalized = convertPointToNormalizedCoords(sx, sy);
      points[pointsIndex++] = normalized.nx;
      points[pointsIndex++] = normalized.ny;
      points[pointsIndex++] = sc.getMarkerSize(index); // size
    }
  
    webGPUInit(canvas, points);
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

function convertPointToNormalizedCoords(
x: number,
y: number
): { nx: number; ny: number } {
// Calculate the center of the viewBox
const centerX = gpuViewboxRect.left + gpuViewboxRect.width / 2;
const centerY = gpuViewboxRect.top + gpuViewboxRect.height / 2;

// Map the x and y coordinates to the normalized coordinate system (-1 to 1)
const nx = (x - centerX) / (gpuViewboxRect.width / 2);
const ny = -((y - centerY) / (gpuViewboxRect.height / 2));

return { nx, ny };
}

async function webGPUInit(webGPUCanvas: HTMLCanvasElement, points: Float32Array) {
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

    @group(0) @binding(0) var<uniform> uni: Uniforms;

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
        vsOut.position = vec4f(vert.position + pos * vert.size / uni.resolution, 0, 1);
        vsOut.texcoord = pos * 0.5 + 0.5;
        return vsOut;
    }

    @group(0) @binding(1) var s: sampler;
    @group(0) @binding(2) var t: texture_2d<f32>;

    @fragment fn fs(vsOut: VSOutput) -> @location(0) vec4f {
        return textureSample(t, s, vsOut.texcoord);
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
    usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
});
gpuDevice.queue.writeBuffer(vertexBuffer, 0, points);

const uniformValues = new Float32Array(2);
const uniformBuffer = gpuDevice.createBuffer({
    size: uniformValues.byteLength,
    usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
});
const kResolutionOffset = 0;
const resolutionValue = uniformValues.subarray(
    kResolutionOffset, kResolutionOffset + 2);

const bindGroup = gpuDevice.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
    { binding: 0, resource: { buffer: uniformBuffer }},
    { binding: 1, resource: sampler },
    { binding: 2, resource: texture.createView() },
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
