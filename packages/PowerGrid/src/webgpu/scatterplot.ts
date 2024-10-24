/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';


class WebGPUCache {
  viewBox: DG.Rect = new DG.Rect(0, 0, 0, 0);
  viewBuffer: GPUBuffer | null = null;
  indexBuffer: GPUBuffer | null = null;
  vertexBuffer: GPUBuffer | null = null;
  indexBufferVersion = -1;
  indexBufferLength = -1;
  columnBuffer: GPUBuffer | null = null;
  xCol: DG.Column | null = null;
  yCol: DG.Column | null = null;
  xColVersion = -1;
  yColVersion = -1;
  xColLength = -1;
  yColLength = -1;
  markerSizesBuffer: GPUBuffer | null = null;
  markerSizesLength = -1;
  markerDefaultSize = -1;
  sizeColumnName = '';
  colorBuffer: GPUBuffer | null = null;
  colorColumnName = '';
  colorLength = -1;
  colorAxisType: keyof typeof DG.AxisType | null = null;
  colorMin = -1;
  colorMax = -1;
  selectedRowsColor = -1;
  filteredRowsColor = -1;
  filteredOutRowsColor = -1;
  invertColorScheme = true;
  linearColorScheme: Array<number> = [];
  categoricalColorScheme: Array<number> = [];
  selectionVersion = -1;
  texture: OffscreenCanvas | null = null;
  minTextureSize = 2;
  maxTextureSize = 100;
  textureGridSize = Math.ceil(Math.sqrt((this.maxTextureSize - this.minTextureSize) / 2 + 1));
  texturePadding = 2;
  gpuTexture: GPUTexture | null = null;

  updateAndValidate(sc: DG.ScatterPlotViewer, device: GPUDevice, pt: DG.Point = new DG.Point(0, 0)) {
    const xCol = sc.table.col(sc.props.xColumnName);
    const yCol = sc.table.col(sc.props.yColumnName);
    if (!xCol || !yCol)
      return false;

    // Storing check results first not to rewrite anything beforehand
    const indexBufferChanged = this.isIndexBufferChanged(sc.filter);
    const columnChanged = this.isColumnChanged(xCol, yCol);
    const markerSizesParamsChanged = this.isMarkerSizesParamsChanged(sc);
    const colorChanged = this.isColorChanged(sc);
  
    if (indexBufferChanged)
        this.setIndexBuffer(sc.filter, device);
  
    if (columnChanged)
        this.setColumns(xCol, yCol, device);
  
    // We'll set the viewBox and viewPort each time as this is cheap
    this.setViewBuffer(sc, device, pt);
  
    if (markerSizesParamsChanged)
        this.setMarkerSizes(sc, device);

    if (colorChanged)
      this.setColor(sc, device);
  
    return this.isValid();
  }

  isViewBoxChanged(viewBox: DG.Rect) {
    return !areEqual(viewBox, this.viewBox);
  }

  isIndexBufferChanged(indexBuffer: DG.BitSet) {
    return this.indexBufferVersion != indexBuffer.version;
  }

  isColumnChanged(xCol: DG.Column, yCol: DG.Column) {
    return !this.columnBuffer || xCol != this.xCol || yCol != this.yCol ||
            !this.xCol || this.xColVersion != xCol.version ||
            !this.yCol || this.yColVersion != yCol.version;
  }

  isMarkerSizesParamsChanged(sc: DG.ScatterPlotViewer) {
    return sc.props.markerDefaultSize != this.markerDefaultSize || sc.props.sizeColumnName != this.sizeColumnName;
  }

  isColorChanged(sc: DG.ScatterPlotViewer) {
    return sc.props.colorColumnName != this.colorColumnName ||
     sc.props.selectedRowsColor != this.selectedRowsColor ||
     sc.props.filteredRowsColor != this.filteredRowsColor ||
     sc.props.filteredOutRowsColor != this.filteredOutRowsColor ||
     sc.props.colorAxisType != this.colorAxisType ||
     sc.props.invertColorScheme != this.invertColorScheme ||
     sc.props.colorMin != this.colorMin ||
     sc.props.colorMax != this.colorMax ||
     sc.dataFrame.selection.version != this.selectionVersion ||
     sc.props.linearColorScheme.length != this.linearColorScheme.length ||
     !sc.props.linearColorScheme.every((c, i) => c == this.linearColorScheme[i]) ||
     sc.props.categoricalColorScheme.length != this.categoricalColorScheme.length ||
     !sc.props.categoricalColorScheme.every((c, i) => c == this.categoricalColorScheme[i]);
  }

  isValid() {
    return this.indexBufferLength > 0 && this.xColLength > 0 && this.yColLength > 0;
  }

  setViewBox(viewBox: DG.Rect) {
    this.viewBox = viewBox;
  }

  setViewBuffer(sc: DG.ScatterPlotViewer, device: GPUDevice, pt: DG.Point = new DG.Point(0, 0)) {
    const kViewBoxByteOffset = 16;
    const kViewPortByteOffset = 16;
    const kPropsByteOffset = 16;
    const kPointByteOffset = 16;

    const inverseX = sc.props.invertXAxis;
    const inverseY = !sc.props.invertYAxis;
    const linearX = sc.props.xAxisType == DG.AxisType.linear;
    const linearY = sc.props.yAxisType == DG.AxisType.linear;

    this.viewBuffer = device.createBuffer({
      size: kViewBoxByteOffset + kViewPortByteOffset + kPropsByteOffset + kPointByteOffset,
      usage: GPUBufferUsage.STORAGE,
      mappedAtCreation: true,
    });
    const viewBufferArray = this.viewBuffer.getMappedRange();
    let viewBufferOffset = 0;
    new Float32Array(viewBufferArray, viewBufferOffset, 4 + 4)
      .set([sc.viewBox.left, sc.viewBox.top, sc.viewBox.width, sc.viewBox.height,
        sc.viewport.left, sc.viewport.top, sc.viewport.width, sc.viewport.height]);
    viewBufferOffset += kViewBoxByteOffset + kViewPortByteOffset;
    new Uint32Array(viewBufferArray, viewBufferOffset, 4)
      .set([+!!inverseX, +!!inverseY, +!!linearX, +!!linearY]);
    viewBufferOffset += kPropsByteOffset;
    new Float32Array(viewBufferArray, viewBufferOffset, 2).set([pt.x, pt.y]);
    this.viewBuffer.unmap();
  }

  setIndexBuffer(indexBuffer: DG.BitSet, device: GPUDevice) {
    this.indexBufferVersion = indexBuffer.version;
    this.indexBufferLength = indexBuffer.length;
    const kFilteredArrayByteOffset = getPaddedSize(this.indexBufferLength);
    const filteredIndexes = indexBuffer.getSelectedIndexes();

    this.indexBuffer = device.createBuffer({
      size: kFilteredArrayByteOffset,
      usage: GPUBufferUsage.STORAGE,
      mappedAtCreation: true,
    });
    new Int32Array(this.indexBuffer.getMappedRange()).set(filteredIndexes);
    this.indexBuffer.unmap();

    this.vertexBuffer = device.createBuffer({
      size: kFilteredArrayByteOffset,
      usage: GPUBufferUsage.VERTEX,
      mappedAtCreation: true,
    });
    new Int32Array(this.vertexBuffer.getMappedRange()).set(filteredIndexes);
    this.vertexBuffer.unmap();
  }

  setColumns(xCol: DG.Column, yCol: DG.Column, device: GPUDevice) {
    this.xCol = xCol;
    this.yCol = yCol;
    const xColumnData = xCol.getRawData();
    const yColumnData = yCol.getRawData();
    this.xColLength = xColumnData.length;
    this.yColLength = yColumnData.length;
    this.xColVersion = xCol.version;
    this.yColVersion = yCol.version;
    this.columnBuffer = device.createBuffer({
      size: getPaddedSize(this.xColLength + this.yColLength),
      usage: GPUBufferUsage.STORAGE,
      mappedAtCreation: true,
    });
    const scBufferArray = this.columnBuffer.getMappedRange();
    let scBufferOffset = 0;
    new Float32Array(scBufferArray, scBufferOffset, this.xColLength).set(xColumnData);
    scBufferOffset += this.xColLength * Float32Array.BYTES_PER_ELEMENT;
    new Float32Array(scBufferArray, scBufferOffset, this.yColLength).set(yColumnData);
    this.columnBuffer.unmap();
  }

  setMarkerSizes(sc: DG.ScatterPlotViewer, device: GPUDevice) {
    this.markerDefaultSize = sc.props.markerDefaultSize;
    this.sizeColumnName = sc.props.sizeColumnName;
    if (!sc.props.sizeColumnName) {
      const size = Math.max(sc.getMarkerSize(0), 2);
      this.markerSizesLength = 1;
      this.markerSizesBuffer = device.createBuffer({
          size: getPaddedSize(this.markerSizesLength),
          usage: GPUBufferUsage.STORAGE,
          mappedAtCreation: true,
        });
      const scBufferArray = this.markerSizesBuffer.getMappedRange();
      new Float32Array(scBufferArray, 0, this.markerSizesLength).set([size]);
      this.markerSizesBuffer.unmap();
    }
    else {
      const sizes = sc.getMarkerSizes();
      this.markerSizesLength = sizes.length;
      this.markerSizesBuffer = device.createBuffer({
          size: getPaddedSize(this.markerSizesLength),
          usage: GPUBufferUsage.STORAGE,
          mappedAtCreation: true,
        });
      const scBufferArray = this.markerSizesBuffer.getMappedRange();
      new Float32Array(scBufferArray, 0, this.markerSizesLength).set(sizes);
      this.markerSizesBuffer.unmap();
    }
  }

  setColor(sc: DG.ScatterPlotViewer, device: GPUDevice) {
    this.colorColumnName = sc.props.colorColumnName;
    const colors = sc.getMarkerColors();
    this.colorLength = colors.length;
    this.selectedRowsColor = sc.props.selectedRowsColor;
    this.filteredRowsColor = sc.props.filteredRowsColor;
    this.filteredOutRowsColor = sc.props.filteredOutRowsColor;
    this.invertColorScheme = sc.props.invertColorScheme;;
    this.colorAxisType = sc.props.colorAxisType;
    this.colorMin = sc.props.colorMin;
    this.colorMax = sc.props.colorMax;
    this.linearColorScheme = sc.props.linearColorScheme;
    this.categoricalColorScheme = sc.props.categoricalColorScheme;
    this.selectionVersion = sc.dataFrame.selection.version;
    this.colorBuffer = device.createBuffer({
        size: getPaddedSize(this.colorLength),
        usage: GPUBufferUsage.STORAGE,
        mappedAtCreation: true,
      });
    const colorsArray = this.colorBuffer.getMappedRange();
    new Uint32Array(colorsArray, 0, this.colorLength).set(colors);
    this.colorBuffer.unmap();
  }

  updateTexuteAtlas(sc: DG.ScatterPlotViewer, device: GPUDevice) {
    if (this.texture == null || this.gpuTexture == null 
      || this.isMarkerSizesParamsChanged(sc)
      || this.minTextureSize != roundUpToEven(sc.props.markerMinSize) 
      || this.maxTextureSize != roundUpToEven(sc.props.markerMaxSize)) {
        if (!sc.props.sizeColumnName) {
          this.minTextureSize = roundUpToEven(sc.props.markerMinSize);
          this.maxTextureSize = roundUpToEven(sc.props.markerMaxSize);
          this.textureGridSize = Math.ceil(Math.sqrt((this.maxTextureSize - this.minTextureSize) / 2 + 1));
          const size = sc.getMarkerSize(0);
          this.texture = createCircleCanvas(size + sc.props.markerBorderWidth * 2, sc);

          this.gpuTexture = device.createTexture({
              size: [this.texture.width, this.texture.height],
              format: 'rgba8unorm',
              usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT,
          });  
          device.queue.copyExternalImageToTexture(
              { source: this.texture },
              { texture: this.gpuTexture, premultipliedAlpha: true},
              [this.texture.width, this.texture.height]
          );
        }
        else {
          this.minTextureSize = roundUpToEven(sc.props.markerMinSize);
          this.maxTextureSize = roundUpToEven(sc.props.markerMaxSize);
          this.textureGridSize = Math.ceil(Math.sqrt((this.maxTextureSize - this.minTextureSize) / 2 + 1));
          this.texture = createTextureAtlas(this, sc);

          this.gpuTexture = device.createTexture({
              size: [this.texture.width, this.texture.height],
              format: 'rgba8unorm',
              usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT,
          });  
          device.queue.copyExternalImageToTexture(
              { source: this.texture },
              { texture: this.gpuTexture, premultipliedAlpha: true},
              [this.texture.width, this.texture.height]
          );
        }
    }
  }
};

const scMap: Map<DG.ScatterPlotViewer, WebGPUCache> = new Map();

function getWebGPUCache(sc: DG.ScatterPlotViewer) {
  if (!scMap.get(sc))
    scMap.set(sc, new WebGPUCache());
  return scMap.get(sc);
}

function updateCanvasPosition(cache: WebGPUCache, viewBox: DG.Rect, canvas: HTMLCanvasElement) {
  if (cache.isViewBoxChanged(viewBox)) {
    cache.setViewBox(viewBox);
    canvas.width = viewBox.width;
    canvas.height = viewBox.height;
    canvas.style.width = `${viewBox.width}px`;
    canvas.style.height = `${viewBox.height}px`;
    canvas.style.left = `${viewBox.left}px`;
    canvas.style.top = `${viewBox.top}px`;
  }
}

export async function scWebGPURender(sc: DG.ScatterPlotViewer, show: boolean) {
  // Getting WebGPU canvas or creating it if absent
  const canvasName = 'webGPUCanvas';
  let canvas = sc.canvas.parentElement?.children.namedItem(canvasName) as HTMLCanvasElement;
  //We need to hide the canvas if we are not drawing on WebGPU
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

  const cache = getWebGPUCache(sc);
  if (cache) {
    updateCanvasPosition(cache, sc.viewBox, canvas);

    try {
        if (sc.props.markerType == DG.MARKER_TYPE.DOT)
          await webGPURenderDots(canvas, sc);
        else
          await webGPURenderTexture(canvas, sc);    
    } catch (error) {
        canvas.hidden = true;
        throw error;
    }
  }
}

export async function scWebGPUPointHitTest(sc: DG.ScatterPlotViewer, pt: DG.Point) : Promise<number> {
  const cache = getWebGPUCache(sc);
  if (!cache)
    throw  'Failed to get WebGPU cache for scatter plot viewer';

  const device = await getGPUDevice();
  if (!device)
    throw  'Failed to get WebGPU device';

  if (!cache.updateAndValidate(sc, device, pt) || !cache.indexBuffer || !cache.columnBuffer || !cache.viewBuffer || !cache.markerSizesBuffer)
    throw 'Failed to update and validate cache or to initalize buffers';

  const kWorkgroupSize = 100;
  const workGroupDispatchSize = Math.ceil(Math.sqrt(Math.ceil(cache.indexBufferLength / kWorkgroupSize)));

  const hitResultBuffer = device.createBuffer({
    size: 4, // A single i32 result
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
    mappedAtCreation: true
  });
  new Int32Array(hitResultBuffer.getMappedRange()).set([-1]); // Initialize to -1 (no hit)
  hitResultBuffer.unmap();

  const module = device.createShaderModule({code: `
            ${addStructures(cache)}    

            @group(0) @binding(0) var<storage, read> indexes: array<i32, ${cache.indexBufferLength}>;
            @group(0) @binding(1) var<storage, read> sc: SC;
            @group(0) @binding(2) var<storage, read> data: Data;
            @group(0) @binding(3) var<storage, read_write> hitResult: atomic<i32>;

            @group(1) @binding(0) var<storage, read> markerSizes: array<f32, ${cache.markerSizesLength}>;

            fn getMarkerType(index: u32) -> i32 {
                return 0; // Dummy function for marker type
            }

            fn hitTest(markerSize: f32, screenPoint: vec2<f32>, markerType: i32) -> bool {
                let dx = sc.mousePoint.x - screenPoint.x;
                let dy = sc.mousePoint.y - screenPoint.y;
                return dx * dx + dy * dy <= markerSize * markerSize;
            }

            ${addPointConversionMethods()}

            @compute @workgroup_size(${kWorkgroupSize})
            fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
                let idx = global_id.y * ${workGroupDispatchSize} * ${kWorkgroupSize} + global_id.x;
                if (idx >= ${cache.indexBufferLength}) {
                    return;
                }

                let filteredIndex = indexes[idx];
                let markerSize = markerSizes[ ${!sc.props.sizeColumnName ? 0 : `filteredIndex`}];
                let screenPoint = pointToScreen(filteredIndex);
                let markerType = getMarkerType(idx);
                if (hitTest(ceil(markerSize) / 2.0, screenPoint, markerType)) {
                    atomicMax(&hitResult, i32(filteredIndex));
                }
            }
        `});
  const pipeline = device.createComputePipeline({
    layout: 'auto',
    compute: {
      module: module,
      entryPoint: 'main',
    },
  });

  const bindGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: cache.indexBuffer}},
      {binding: 1, resource: {buffer: cache.viewBuffer}},
      {binding: 2, resource: {buffer: cache.columnBuffer}},
      {binding: 3, resource: {buffer: hitResultBuffer}}
    ],
  });

  const dataGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(1),
    entries: [
      {binding: 0, resource: {buffer: cache.markerSizesBuffer}},
    ],
  });

  const commandEncoder = device.createCommandEncoder();
  const passEncoder = commandEncoder.beginComputePass();
  passEncoder.setPipeline(pipeline);
  passEncoder.setBindGroup(0, bindGroup);
  passEncoder.setBindGroup(1, dataGroup);
  passEncoder.dispatchWorkgroups(workGroupDispatchSize, workGroupDispatchSize);
  passEncoder.end();

  // Copy hit result to a readable buffer
  const readbackBuffer = device.createBuffer({
    size: 4, // size for i32 result
    usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
  });
  commandEncoder.copyBufferToBuffer(hitResultBuffer, 0, readbackBuffer, 0, 4);

  const commandBuffer = commandEncoder.finish();
  device.queue.submit([commandBuffer]);

  // Wait for the result and read hit point index
  await readbackBuffer.mapAsync(GPUMapMode.READ);
  const resultArray = new Int32Array(readbackBuffer.getMappedRange());
  const hitIndex = resultArray[0];
  readbackBuffer.unmap();

  const info = await module.getCompilationInfo();
  for (const message of info.messages) {
    if (message.type === 'error') {
      throw `WebGPU hit test shader module error: ${message.message}`;
    }
  }

  return hitIndex;
}

function areEqual(rc1: DG.Rect, rc2: DG.Rect): boolean {
  return Math.floor(rc1.left) == Math.floor(rc2.left) && Math.floor(rc1.top) == Math.floor(rc2.top) &&
        Math.floor(rc1.width) == Math.floor(rc2.width) && Math.floor(rc1.height) == Math.floor(rc2.height);
}

async function webGPURenderDots(webGPUCanvas: HTMLCanvasElement, sc: DG.ScatterPlotViewer) {
  const cache = getWebGPUCache(sc);
  if (!cache)
    throw  'Failed to get WebGPU cache for scatter plot viewer';

  const device = await getGPUDevice();
  if (!device)
    throw  'Failed to get WebGPU device';

  if (!cache.updateAndValidate(sc, device) || !cache.indexBuffer || !cache.columnBuffer || !cache.viewBuffer || !cache.colorBuffer)
    throw 'Failed to update and validate cache or to initalize buffers';

  const gpuContext = webGPUCanvas.getContext('webgpu');
  if (!gpuContext)
    throw 'Failed to get gpu context from canvas';

  const presentationFormat = navigator.gpu.getPreferredCanvasFormat();
  gpuContext.configure({
    device: device,
    format: presentationFormat,
    alphaMode: 'premultiplied',
  });

  // Defining textures that will be drawn
  const module = device.createShaderModule({
    code: `
        struct Vertex {
            @location(0) index: i32,
        };

        struct VSOutput {
            @builtin(position) position: vec4f,
            @location(0) @interpolate(flat) markerIndex: u32, // This stores the marker index for the fragment shader
        };

        ${addStructures(cache)}

        @group(0) @binding(0) var<storage, read> sc: SC;
        @group(0) @binding(1) var<storage, read> data: Data;
        @group(0) @binding(2) var<storage, read> markerColors: array<u32, ${cache.colorLength}>;

        ${addDotsRendering()}

        ${addPointConversionMethods()}

        ${addColorMethods()}
        `,
  });

  const pipeline = device.createRenderPipeline({
    label: '1 pixel points',
    layout: 'auto',
    vertex: {
      module,
      buffers: [
        {
          arrayStride: 4, // 1 int, 4 bytes
          stepMode: 'instance',
          attributes: [
            {shaderLocation: 0, offset: 0, format: 'sint32'}, // position
          ],
        },
      ],
    },
    fragment: {
      module,
      targets: [{ format: presentationFormat }],
    },
    primitive: {
      topology: 'point-list',
    },
  });

  const dataGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: cache.viewBuffer}},
      {binding: 1, resource: {buffer: cache.columnBuffer}},
      {binding: 2, resource: {buffer: cache.colorBuffer}},
    ],
  });

  const renderPassDescriptor = {
    label: 'Dot canvas renderPass',
    colorAttachments: [
      {
        loadOp: 'clear',
        storeOp: 'store',
        view: gpuContext.getCurrentTexture().createView(),
      },
    ],
  };

  const encoder = device.createCommandEncoder();
  const pass = encoder.beginRenderPass(renderPassDescriptor as GPURenderPassDescriptor);
  pass.setPipeline(pipeline);
  pass.setVertexBuffer(0, cache.vertexBuffer);
  pass.setBindGroup(0, dataGroup);
  pass.draw(1, cache.indexBufferLength);
  pass.end();

  const encoderBuffer = encoder.finish();
  device.queue.submit([encoderBuffer]);
  await device.queue.onSubmittedWorkDone();

  const info = await module.getCompilationInfo();
  for (const message of info.messages) {
    if (message.type === 'error') {
      throw `WebGPURender shader module error: ${message.message}`;
    }
  }
}

async function webGPURenderTexture(webGPUCanvas: HTMLCanvasElement, sc: DG.ScatterPlotViewer) {
  const cache = getWebGPUCache(sc);
  if (!cache)
    throw  'Failed to get WebGPU cache for scatter plot viewer';

  const device = await getGPUDevice();
  if (!device)
    throw  'Failed to get WebGPU device';

  cache.updateTexuteAtlas(sc, device);
  if (!cache.gpuTexture)
    throw 'Failed to update texture atlas';

  if (!cache.updateAndValidate(sc, device) || !cache.indexBuffer || !cache.columnBuffer || !cache.viewBuffer || !cache.markerSizesBuffer || !cache.colorBuffer)
    throw 'Failed to update and validate cache or to initalize buffers';

  const gpuContext = webGPUCanvas.getContext('webgpu');
  if (!gpuContext)
    throw 'Failed to get gpu context from canvas';

  const presentationFormat = navigator.gpu.getPreferredCanvasFormat();
  gpuContext.configure({
    device: device,
    format: presentationFormat,
    alphaMode: 'premultiplied',
  });

  // Defining textures that will be drawn
  const module = device.createShaderModule({
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
            @location(1) @interpolate(flat) markerIndex: u32, // This stores the marker index for the fragment shader
            @location(2) @interpolate(flat) sizeIndex: u32, // This stores the marker size index for the fragment shader
        };

        ${addStructures(cache)}

        @group(0) @binding(0) var<uniform> uni: Uniforms;
        @group(0) @binding(1) var s: sampler;
        @group(0) @binding(2) var t: texture_2d<f32>;
        
        @group(1) @binding(0) var<storage, read> sc: SC;
        @group(1) @binding(1) var<storage, read> data: Data;
        @group(1) @binding(2) var<storage, read> markerSizes: array<f32, ${cache.markerSizesLength}>;
        @group(1) @binding(3) var<storage, read> markerColors: array<u32, ${cache.colorLength}>;

        ${!sc.props.sizeColumnName ? addSingleMarkerSizeRendering(cache, sc) : addDifferentMarkerSizesRendering(cache, sc)}

        ${addPointConversionMethods()}

        ${addColorMethods()}
        `,
  });

  const pipeline = device.createRenderPipeline({
    label: 'sizeable points with texture',
    layout: 'auto',
    vertex: {
      module,
      buffers: [
        {
          arrayStride: 4, // 1 int, 4 bytes
          stepMode: 'instance',
          attributes: [
            {shaderLocation: 0, offset: 0, format: 'sint32'}, // position
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

  const sampler = device.createSampler({
    minFilter: 'nearest',
    magFilter: 'nearest',
  });

  const uniformValues = new Float32Array(2);
  const uniformBuffer = device.createBuffer({
    size: uniformValues.byteLength,
    usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
  });
  const kResolutionOffset = 0;
  const resolutionValue = uniformValues.subarray(
    kResolutionOffset, kResolutionOffset + 2);

  const dataGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(1),
    entries: [
      {binding: 0, resource: {buffer: cache.viewBuffer}},
      {binding: 1, resource: {buffer: cache.columnBuffer}},
      {binding: 2, resource: {buffer: cache.markerSizesBuffer}},
      {binding: 3, resource: {buffer: cache.colorBuffer}},
    ],
  });

  const renderPassDescriptor = {
    label: 'our basic canvas renderPass',
    colorAttachments: [
      {
        loadOp: 'clear',
        storeOp: 'store',
        view: gpuContext.getCurrentTexture().createView(),
      },
    ],
  };

  // Get the current texture from the canvas context and
  // set it as the texture to render to.
  const canvasTexture = gpuContext.getCurrentTexture();
  // Update the resolution in the uniform buffer
  resolutionValue.set([canvasTexture.width, canvasTexture.height]);
  device.queue.writeBuffer(uniformBuffer, 0, uniformValues);

  const renderGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: uniformBuffer}},
      {binding: 1, resource: sampler},
      {binding: 2, resource: cache.gpuTexture.createView()},
    ],
  });

  const encoder = device.createCommandEncoder();
  const pass = encoder.beginRenderPass(renderPassDescriptor as GPURenderPassDescriptor);
  pass.setPipeline(pipeline);
  pass.setVertexBuffer(0, cache.vertexBuffer);
  pass.setBindGroup(0, renderGroup);
  pass.setBindGroup(1, dataGroup);
  pass.draw(6, cache.indexBufferLength);
  pass.end();

  const encoderBuffer = encoder.finish();
  device.queue.submit([encoderBuffer]);
  await device.queue.onSubmittedWorkDone();

  const info = await module.getCompilationInfo();
  for (const message of info.messages) {
    if (message.type === 'error') {
      throw `WebGPURender shader module error: ${message.message}`;
    }
  }
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

function createTextureAtlas(cache: WebGPUCache, sc: DG.ScatterPlotViewer): OffscreenCanvas {
    const minSize = cache.minTextureSize + sc.props.markerBorderWidth * 2;
    const maxSize = cache.maxTextureSize + sc.props.markerBorderWidth * 2;
    // We'll take the minimum size as 2 while adding the border width to maintain the whole size
    const sizes = Array.from({ length: (maxSize - minSize) / 2 + 1 }, (_, i) => 2 * i + minSize);
  
    // Calculate the size of the atlas canvas
    const atlasSize = (maxSize + cache.texturePadding) * cache.textureGridSize; // Each cell will have a size of maxSize + padding
    const atlasCanvas = new OffscreenCanvas(atlasSize, atlasSize);
    const ctx = atlasCanvas.getContext('2d');
  
    if (!ctx)
        return atlasCanvas;
  
    // Draw each circle in the grid
    sizes.forEach((size, index) => {
        const x = index % cache.textureGridSize;
        const y = Math.floor(index / cache.textureGridSize);
    
        // Get the canvas for the current circle size
        const circleCanvas = createCircleCanvas(size, sc);
    
        // Calculate position in the atlas
        const cellSize = maxSize + cache.texturePadding;
        // Center the texture in the cell
        const posX = x * cellSize + (cellSize - size) / 2; 
        const posY = y * cellSize + (cellSize - size) / 2;
    
        // Draw the texture onto the atlas at the calculated position
        ctx.drawImage(circleCanvas, posX, posY);
    });
  
    return atlasCanvas;
  }

function createCircleCanvas(size: number, sc: DG.ScatterPlotViewer): OffscreenCanvas {
  const lineWidth = sc.props.markerBorderWidth;
  const canvas = new OffscreenCanvas(size, size);
  const ctx = canvas.getContext('2d');
  // Align center to the pixel grid
  const centerX = Math.floor(size / 2) + 0.5;
  const centerY = Math.floor(size / 2) + 0.5;
  // Ensure radius fits within pixel grid
  const radius = Math.floor((size - lineWidth) / 2);
  if (ctx) {
    ctx.beginPath();
    ctx.arc(centerX, centerY, radius, 0, 2 * Math.PI, false);
    ctx.save();
    ctx.clip();
    ctx.fillStyle = "#0000ff";
    ctx.fill();
    ctx.lineWidth = lineWidth * 2;
    ctx.strokeStyle = "#00003f";
    ctx.stroke();
    ctx.restore();
  }

  return canvas;
}

function roundUpToEven(num: number): number {
  return Math.ceil(num / 2) * 2;
}

function addPointConversionMethods() {
  return `
        fn pointToScreen(index: i32) -> vec2<f32> {
            var sx = sc.viewBox.left + worldToScreen(data.xColumnData[index], sc.viewBox.width, sc.viewport.left, sc.viewport.left + sc.viewport.width, sc.props.inverseX, sc.props.linearX);
            var sy = sc.viewBox.top + worldToScreen(data.yColumnData[index], sc.viewBox.height, sc.viewport.top, sc.viewport.top + sc.viewport.height, sc.props.inverseY, sc.props.linearY);
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
    `;
}

function addStructures(cache: WebGPUCache) {
  return `
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
            viewBox: Rect,
            viewport: Rect,
            props: Props,
            mousePoint: vec2<f32>,
        };

        struct Data {
            xColumnData: array<f32, ${cache.xColLength}>,
            yColumnData: array<f32, ${cache.yColLength}>,
        };
    `;
}

function addDifferentMarkerSizesRendering(cache: WebGPUCache, sc: DG.ScatterPlotViewer) {
  return `
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

          // Get the marker size for the current vertex
          let markerSize = markerSizes[vert.index];
          
          let minSize = f32(${cache.minTextureSize});
          let maxSize = f32(${cache.maxTextureSize});

          // The textures are made of even sizes to avoid blur and artefacts
          // Rounding up to the nearest even index and dividing by two to get the needed index in the texture atlas
          let sizeIndex = u32(ceil(ceil(markerSize) / 2.0) * 2.0 - minSize) / 2;

          var vsOut: VSOutput;
          let pos = points[vNdx];

          let screenPoint = pointToScreen(vert.index);
          let normalizedPos = convertPointToNormalizedCoords(screenPoint);
          // Making a pixel perfect position, to avoid artefacts and blurring
          vsOut.position = vec4f(floor((normalizedPos + pos * (maxSize + ${sc.props.markerBorderWidth * 2 + cache.texturePadding}) / uni.resolution) * uni.resolution) / uni.resolution, 0, 1);
          vsOut.texcoord = pos * 0.5 + 0.5;
          vsOut.markerIndex = u32(vert.index);
          vsOut.sizeIndex = sizeIndex;
          return vsOut;
      }

      @fragment fn fs(vsOut: VSOutput) -> @location(0) vec4f {
          let gridSize: u32 = ${cache.textureGridSize};

          // Get the size index based on the marker index (you might want a mapping function here)
          let sizeIndex: u32 = vsOut.sizeIndex;

          // Calculate (x, y) in the texture atlas grid
          let x: u32 = sizeIndex % gridSize;
          let y: u32 = sizeIndex / gridSize;

          // Calculate UV offset for the selected size
          let uvOffset = vec2f(f32(x) * (1.0 / f32(gridSize)), f32(y) * (1.0 / f32(gridSize)));
          let uvScale = 1.0 / f32(gridSize);

          // Adjust the texcoords to the right portion of the atlas
          let texCoords = vsOut.texcoord * uvScale + uvOffset;

          let color = textureSample(t, s, texCoords);
          var c = markerColors[vsOut.markerIndex];
          
          // color.b is from 0 to 1 - basically it measures intencity, as stroke is darker
          if (color.b > 0.0) {
            return vec4f(r(c) * color.b, g(c) * color.b, b(c) * color.b, color.a);  
          }
          return color;
      }
  `;
}

function addSingleMarkerSizeRendering(cache: WebGPUCache, sc: DG.ScatterPlotViewer) {
  return `
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

          // Get the marker size for the current vertex
          let markerSize = markerSizes[0];
          
          var vsOut: VSOutput;
          let pos = points[vNdx];

          let screenPoint = pointToScreen(vert.index);
          let normalizedPos = convertPointToNormalizedCoords(screenPoint);
          // Making a pixel perfect position, to avoid artefacts and blurring
          vsOut.position = vec4f(normalizedPos + pos * (markerSize + ${sc.props.markerBorderWidth * 2}) / uni.resolution, 0, 1);
          vsOut.texcoord = pos * 0.5 + 0.5;
          vsOut.markerIndex = u32(vert.index);   // Pass marker index to fragment shader
          vsOut.sizeIndex = 0;
          return vsOut;
      }

      @fragment fn fs(vsOut: VSOutput) -> @location(0) vec4f {
          let color = textureSample(t, s, vsOut.texcoord);
          var c = markerColors[vsOut.markerIndex];

          // color.b is from 0 to 1 - basically it measures intencity, as stroke is darker
          if (color.b > 0.0) {
            return vec4f(r(c) * color.b, g(c) * color.b, b(c) * color.b, color.a);  
          }
          return color;
      }
  `;
}

function addDotsRendering() {
  return `
      @vertex fn vs(vert: Vertex) -> VSOutput {
          var vsOut: VSOutput;

          let screenPoint = pointToScreen(vert.index) + 0.5;
          let normalizedPos = convertPointToNormalizedCoords(screenPoint);
          // Making a pixel perfect position, to avoid artefacts and blurring
          vsOut.position = vec4f(normalizedPos, 0, 1);
          vsOut.markerIndex = u32(vert.index);   // Pass marker index to fragment shader
          return vsOut;
      }

      @fragment fn fs(vsOut: VSOutput) -> @location(0) vec4f {
          var c = markerColors[vsOut.markerIndex];
          return vec4f(r(c), g(c), b(c), a(c));
      }
  `;
}

function addColorMethods() {
  return `
      fn a(c: u32) -> f32 {    
        return f32((c >> 24u) & 0xFFu) / 255.0;
      }

      fn r(c: u32) -> f32 {
        return f32((c >> 16u) & 0xFFu) / 255.0;
      }

      fn g(c: u32) -> f32 {
        return f32((c >> 8u) & 0xFFu) / 255.0;
      }

      fn b(c: u32) -> f32 {
        return f32(c & 0xFFu) / 255.0;
      }

      fn clamp(value: i32, minVal: i32, maxVal: i32) -> u32 {
          return u32(max(min(value, maxVal), minVal));
      }

      fn darken(c: u32, diff: i32) -> u32 {
          let a: u32 = clamp(i32((c >> 24) & 0xFF), 0, 255);
          let r: u32 = clamp(i32((c >> 16) & 0xFF) - diff, 0, 255);
          let g: u32 = clamp(i32((c >> 8) & 0xFF) - diff, 0, 255);
          let b: u32 = clamp(i32(c & 0xFF) - diff, 0, 255);

          return u32((a << 24) | (r << 16) | (g << 8) | b);
      }
  `;
}