var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import * as grok from 'datagrok-api/grok';
export function ensureContainerRunning(containerName, timeout = 300000) {
    return __awaiter(this, void 0, void 0, function* () {
        const container = yield grok.dapi.docker.dockerContainers.filter(containerName).first();
        if (!(container.status.startsWith('started') || container.status.startsWith('checking'))) {
            console.log(`starting container ${container.name}`);
            yield grok.dapi.docker.dockerContainers.run(container.id, false);
        }
        const deadline = Date.now() + timeout;
        let restartCount = 0;
        while (Date.now() < deadline) {
            const cont = yield grok.dapi.docker.dockerContainers.find(container.id);
            if (cont.status.startsWith('started') || cont.status.startsWith('checking'))
                return;
            if (cont.status === 'error' || cont.status.includes('stopped')) {
                if (restartCount < 3) {
                    restartCount++;
                    yield grok.dapi.docker.dockerContainers.run(container.id, false);
                }
                else
                    throw new Error(`Container ${containerName} failed: ${cont.status} after ${restartCount} restart attempts`);
            }
            yield new Promise((r) => setTimeout(r, 5000));
        }
        throw new Error(`${containerName} hasn't been started after ${timeout / 60000} minutes`);
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC1jb250YWluZXItdXRpbHMuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyJ0ZXN0LWNvbnRhaW5lci11dGlscy50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7Ozs7Ozs7QUFBQSxPQUFPLEtBQUssSUFBSSxNQUFNLG1CQUFtQixDQUFDO0FBRTFDLE1BQU0sVUFBZ0Isc0JBQXNCLENBQUMsYUFBcUIsRUFBRSxVQUFrQixNQUFNOztRQUMxRixNQUFNLFNBQVMsR0FBRyxNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsTUFBTSxDQUFDLGdCQUFnQixDQUFDLE1BQU0sQ0FBQyxhQUFhLENBQUMsQ0FBQyxLQUFLLEVBQUUsQ0FBQztRQUN4RixJQUFJLENBQUMsQ0FBQyxTQUFTLENBQUMsTUFBTSxDQUFDLFVBQVUsQ0FBQyxTQUFTLENBQUMsSUFBSSxTQUFTLENBQUMsTUFBTSxDQUFDLFVBQVUsQ0FBQyxVQUFVLENBQUMsQ0FBQyxFQUFFO1lBQ3hGLE9BQU8sQ0FBQyxHQUFHLENBQUMsc0JBQXNCLFNBQVMsQ0FBQyxJQUFJLEVBQUUsQ0FBQyxDQUFDO1lBQ3BELE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxNQUFNLENBQUMsZ0JBQWdCLENBQUMsR0FBRyxDQUFDLFNBQVMsQ0FBQyxFQUFFLEVBQUUsS0FBSyxDQUFDLENBQUM7U0FDbEU7UUFFRCxNQUFNLFFBQVEsR0FBRyxJQUFJLENBQUMsR0FBRyxFQUFFLEdBQUcsT0FBTyxDQUFDO1FBQ3RDLElBQUksWUFBWSxHQUFHLENBQUMsQ0FBQztRQUNyQixPQUFPLElBQUksQ0FBQyxHQUFHLEVBQUUsR0FBRyxRQUFRLEVBQUU7WUFDNUIsTUFBTSxJQUFJLEdBQUcsTUFBTSxJQUFJLENBQUMsSUFBSSxDQUFDLE1BQU0sQ0FBQyxnQkFBZ0IsQ0FBQyxJQUFJLENBQUMsU0FBUyxDQUFDLEVBQUUsQ0FBUSxDQUFDO1lBQy9FLElBQUksSUFBSSxDQUFDLE1BQU0sQ0FBQyxVQUFVLENBQUMsU0FBUyxDQUFDLElBQUksSUFBSSxDQUFDLE1BQU0sQ0FBQyxVQUFVLENBQUMsVUFBVSxDQUFDO2dCQUN6RSxPQUFPO1lBQ1QsSUFBSSxJQUFJLENBQUMsTUFBTSxLQUFLLE9BQU8sSUFBSSxJQUFJLENBQUMsTUFBTSxDQUFDLFFBQVEsQ0FBQyxTQUFTLENBQUMsRUFBRTtnQkFDOUQsSUFBSSxZQUFZLEdBQUcsQ0FBQyxFQUFFO29CQUNwQixZQUFZLEVBQUUsQ0FBQztvQkFDZixNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsTUFBTSxDQUFDLGdCQUFnQixDQUFDLEdBQUcsQ0FBQyxTQUFTLENBQUMsRUFBRSxFQUFFLEtBQUssQ0FBQyxDQUFDO2lCQUNsRTs7b0JBRUMsTUFBTSxJQUFJLEtBQUssQ0FBQyxhQUFhLGFBQWEsWUFBWSxJQUFJLENBQUMsTUFBTSxVQUFVLFlBQVksbUJBQW1CLENBQUMsQ0FBQzthQUMvRztZQUNELE1BQU0sSUFBSSxPQUFPLENBQUMsQ0FBQyxDQUFDLEVBQUUsRUFBRSxDQUFDLFVBQVUsQ0FBQyxDQUFDLEVBQUUsSUFBSSxDQUFDLENBQUMsQ0FBQztTQUMvQztRQUNELE1BQU0sSUFBSSxLQUFLLENBQUMsR0FBRyxhQUFhLDhCQUE4QixPQUFPLEdBQUcsS0FBSyxVQUFVLENBQUMsQ0FBQztJQUMzRixDQUFDO0NBQUEiLCJzb3VyY2VzQ29udGVudCI6WyJpbXBvcnQgKiBhcyBncm9rIGZyb20gJ2RhdGFncm9rLWFwaS9ncm9rJztcclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBlbnN1cmVDb250YWluZXJSdW5uaW5nKGNvbnRhaW5lck5hbWU6IHN0cmluZywgdGltZW91dDogbnVtYmVyID0gMzAwMDAwKSB7XHJcbiAgY29uc3QgY29udGFpbmVyID0gYXdhaXQgZ3Jvay5kYXBpLmRvY2tlci5kb2NrZXJDb250YWluZXJzLmZpbHRlcihjb250YWluZXJOYW1lKS5maXJzdCgpO1xyXG4gIGlmICghKGNvbnRhaW5lci5zdGF0dXMuc3RhcnRzV2l0aCgnc3RhcnRlZCcpIHx8IGNvbnRhaW5lci5zdGF0dXMuc3RhcnRzV2l0aCgnY2hlY2tpbmcnKSkpIHtcclxuICAgIGNvbnNvbGUubG9nKGBzdGFydGluZyBjb250YWluZXIgJHtjb250YWluZXIubmFtZX1gKTtcclxuICAgIGF3YWl0IGdyb2suZGFwaS5kb2NrZXIuZG9ja2VyQ29udGFpbmVycy5ydW4oY29udGFpbmVyLmlkLCBmYWxzZSk7XHJcbiAgfVxyXG5cclxuICBjb25zdCBkZWFkbGluZSA9IERhdGUubm93KCkgKyB0aW1lb3V0O1xyXG4gIGxldCByZXN0YXJ0Q291bnQgPSAwO1xyXG4gIHdoaWxlIChEYXRlLm5vdygpIDwgZGVhZGxpbmUpIHtcclxuICAgIGNvbnN0IGNvbnQgPSBhd2FpdCBncm9rLmRhcGkuZG9ja2VyLmRvY2tlckNvbnRhaW5lcnMuZmluZChjb250YWluZXIuaWQpIGFzIGFueTtcclxuICAgIGlmIChjb250LnN0YXR1cy5zdGFydHNXaXRoKCdzdGFydGVkJykgfHwgY29udC5zdGF0dXMuc3RhcnRzV2l0aCgnY2hlY2tpbmcnKSlcclxuICAgICAgcmV0dXJuO1xyXG4gICAgaWYgKGNvbnQuc3RhdHVzID09PSAnZXJyb3InIHx8IGNvbnQuc3RhdHVzLmluY2x1ZGVzKCdzdG9wcGVkJykpIHtcclxuICAgICAgaWYgKHJlc3RhcnRDb3VudCA8IDMpIHtcclxuICAgICAgICByZXN0YXJ0Q291bnQrKztcclxuICAgICAgICBhd2FpdCBncm9rLmRhcGkuZG9ja2VyLmRvY2tlckNvbnRhaW5lcnMucnVuKGNvbnRhaW5lci5pZCwgZmFsc2UpO1xyXG4gICAgICB9XHJcbiAgICAgIGVsc2VcclxuICAgICAgICB0aHJvdyBuZXcgRXJyb3IoYENvbnRhaW5lciAke2NvbnRhaW5lck5hbWV9IGZhaWxlZDogJHtjb250LnN0YXR1c30gYWZ0ZXIgJHtyZXN0YXJ0Q291bnR9IHJlc3RhcnQgYXR0ZW1wdHNgKTtcclxuICAgIH1cclxuICAgIGF3YWl0IG5ldyBQcm9taXNlKChyKSA9PiBzZXRUaW1lb3V0KHIsIDUwMDApKTtcclxuICB9XHJcbiAgdGhyb3cgbmV3IEVycm9yKGAke2NvbnRhaW5lck5hbWV9IGhhc24ndCBiZWVuIHN0YXJ0ZWQgYWZ0ZXIgJHt0aW1lb3V0IC8gNjAwMDB9IG1pbnV0ZXNgKTtcclxufVxyXG4iXX0=