import time
import psutil
import os

# Setup
n_chunks = 100
chunk_size = 3_024_000  # ~3 MB per chunk
chunks = [bytes([i % 256]) * chunk_size for i in range(n_chunks)]
total_size = n_chunks * chunk_size

process = psutil.Process(os.getpid())

# Methods
def method_bytes():
    ch = []
    for chunk in chunks:
        ch.append(chunk)
    return b''.join(ch)

def method_extend():
    result = bytearray()
    for chunk in chunks:
        result.extend(chunk)
    return result

def method_prealloc_slice():
    buffer = bytearray(total_size)
    offset = 0
    for chunk in chunks:
        buffer[offset:offset + len(chunk)] = chunk
        offset += len(chunk)
    return buffer

def method_prealloc_memoryview():
    buffer = bytearray(total_size)
    mv = memoryview(buffer)
    offset = 0
    for chunk in chunks:
        mv[offset:offset + len(chunk)] = chunk
        offset += len(chunk)
    return buffer

# Benchmarking with psutil memory usage
def benchmark(func, name):
    mem_before = process.memory_info().rss / (1024 * 1024)  # MiB
    start_time = time.time()
    result = func()
    elapsed = time.time() - start_time
    mem_after = process.memory_info().rss / (1024 * 1024)   # MiB
    mem_diff = mem_after - mem_before
    return name, elapsed, mem_diff

results = [
    benchmark(method_bytes, "Bytes join"),
    benchmark(method_extend, "Extend"),
    benchmark(method_prealloc_slice, "Prealloc + Slice"),
    benchmark(method_prealloc_memoryview, "Prealloc + MemoryView")
]

results.sort(key=lambda x: x[1])  # Sort by elapsed time

for name, elapsed, mem_diff in results:
    print(f"{name:25s} | Time: {elapsed:.4f}s | Memory change: {mem_diff:+.2f} MiB")
