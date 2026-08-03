var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import { delay, expect, testEvent, testEventAsync } from "./test";
export function createViewer(tv, v, packageName) {
    return __awaiter(this, void 0, void 0, function* () {
        let res;
        if (packageName) {
            res = (yield tv.dataFrame.plot.fromType(v));
            tv.dockManager.dock(res);
        }
        else
            res = tv.addViewer(v);
        return res;
    });
}
export function testViewerInternal(tv, viewerName, packageName, event, actions, awaitViewer, layout, actionArgs, actionsWithDelay = true) {
    return __awaiter(this, void 0, void 0, function* () {
        let actionsRes = null;
        yield testEventAsync(event, (e) => __awaiter(this, void 0, void 0, function* () {
            let viewer = null;
            for (const v1 of tv.viewers) {
                if (v1.type === viewerName)
                    viewer = v1;
            }
            if (!viewer)
                throw Error('Viewer hasn\'t been added');
            yield Promise.resolve(); //re schedules subsequent commands into microtask
            if (awaitViewer)
                yield awaitViewer(viewer);
            if (actions) {
                const args = actionArgs !== null && actionArgs !== void 0 ? actionArgs : {};
                args.tv = tv;
                args.viewer = viewer;
                actionsRes = yield actions(args, actionsWithDelay);
            }
            //check that there are no active subscriptions in the viewer after close
            yield testEvent(grok.events.onViewerClosed, () => {
                expect(viewer.subs.some((s) => !s.closed), false);
            }, () => viewer.close(), 3000);
        }), () => __awaiter(this, void 0, void 0, function* () {
            layout ? tv.loadLayout(layout) : yield createViewer(tv, viewerName, packageName);
        }), 60000, 'TEST_EVENT_ASYNC');
        if (actionsRes)
            return actionsRes;
    });
}
export function selectFilterChangeCurrent(args, withDelay = true) {
    return __awaiter(this, void 0, void 0, function* () {
        const currentDf = args.tv.dataFrame;
        const dfSaved = currentDf.clone();
        //remove values in the first row
        Array.from(currentDf.row(0).cells).forEach((c) => c.value = null);
        //selection
        const num = currentDf.rowCount < 20 ? Math.floor(currentDf.rowCount / 2) : 10;
        currentDf.rows.select((row) => row.idx >= 0 && row.idx < num);
        if (withDelay)
            yield delay(50);
        //filter
        for (let i = num; i < num * 2; i++)
            currentDf.filter.set(i, false);
        if (withDelay)
            yield delay(50);
        //change current row
        currentDf.currentRowIdx = 1;
        //remove columns
        currentDf.columns.names().slice(0, Math.ceil(currentDf.columns.length / 2))
            .forEach((c) => currentDf.columns.remove(c));
        if (withDelay)
            yield delay(100);
        //set back initial df with whole set of columns and preserved data
        args.tv.dataFrame = dfSaved;
        yield delay(50);
    });
}
export function filterAsync(args, withDelay = true) {
    return __awaiter(this, void 0, void 0, function* () {
        const currentDf = args.tv.dataFrame;
        setTimeout(() => currentDf.filter.set(0, !currentDf.filter.get(0)), 0);
    });
}
export function changeOptionsSaveLayout(args, withDelay = true) {
    return __awaiter(this, void 0, void 0, function* () {
        //get current options and properties
        let optns;
        try {
            optns = args.viewer.getOptions(true).look;
        }
        catch (err) {
            //@ts-ignore
            throw new Error(`Viewer's .getOptions() error.`, { cause: err });
        }
        let props;
        try {
            props = args.viewer.getProperties();
        }
        catch (err) {
            //@ts-ignore
            throw new Error(`Viewer's .getProperties() error.`, { cause: err });
        }
        //change options and properties
        const newProps = {};
        Object.keys(optns).filter((k) => typeof optns[k] === 'boolean').forEach((k) => newProps[k] = !optns[k]);
        props.filter((p) => p.choices !== null)
            .forEach((p) => newProps[p.name] = p.choices.find((c) => c !== optns[p.name]));
        //set new options
        args.viewer.setOptions(newProps);
        yield delay(300);
        const layout = args.tv.saveLayout();
        const savedLook = args.viewer.getOptions().look;
        return { layout: layout, savedProps: savedLook };
    });
}
export function loadLayout(args, withDelay = true) {
    return __awaiter(this, void 0, void 0, function* () {
        expect(JSON.stringify(args.viewer.getOptions().look), JSON.stringify(args.savedProps));
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC12aWV3ZXItdXRpbHMuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyJ0ZXN0LXZpZXdlci11dGlscy50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7Ozs7Ozs7QUFHQSxPQUFPLEVBQUUsS0FBSyxFQUFFLE1BQU0sRUFBRSxTQUFTLEVBQUUsY0FBYyxFQUFFLE1BQU0sUUFBUSxDQUFDO0FBR2xFLE1BQU0sVUFBZ0IsWUFBWSxDQUFDLEVBQWlCLEVBQUUsQ0FBUyxFQUFFLFdBQW9COztRQUNqRixJQUFJLEdBQWUsQ0FBQztRQUNwQixJQUFJLFdBQVcsRUFBRTtZQUNiLEdBQUcsSUFBRyxNQUFNLEVBQUUsQ0FBQyxTQUFTLENBQUMsSUFBSSxDQUFDLFFBQVEsQ0FBQyxDQUFDLENBQWUsQ0FBQSxDQUFDO1lBQ3hELEVBQUUsQ0FBQyxXQUFXLENBQUMsSUFBSSxDQUFDLEdBQUcsQ0FBQyxDQUFDO1NBQzVCOztZQUNHLEdBQUcsR0FBRyxFQUFFLENBQUMsU0FBUyxDQUFDLENBQUMsQ0FBQyxDQUFDO1FBQzFCLE9BQU8sR0FBRyxDQUFDO0lBQ2YsQ0FBQztDQUFBO0FBRUQsTUFBTSxVQUFnQixrQkFBa0IsQ0FBQyxFQUFpQixFQUFFLFVBQWtCLEVBQUUsV0FBbUIsRUFDL0YsS0FBc0IsRUFBRSxPQUF5RCxFQUNqRixXQUFtRCxFQUFFLE1BQXVCLEVBQzVFLFVBQWdCLEVBQUUsZ0JBQWdCLEdBQUcsSUFBSTs7UUFDekMsSUFBSSxVQUFVLEdBQVEsSUFBSSxDQUFDO1FBQzNCLE1BQU0sY0FBYyxDQUFDLEtBQUssRUFBRSxDQUFPLENBQU0sRUFBRSxFQUFFO1lBQ3pDLElBQUksTUFBTSxHQUFzQixJQUFJLENBQUM7WUFDckMsS0FBSyxNQUFNLEVBQUUsSUFBSSxFQUFHLENBQUMsT0FBTyxFQUFFO2dCQUMxQixJQUFJLEVBQUUsQ0FBQyxJQUFJLEtBQUssVUFBVTtvQkFDdEIsTUFBTSxHQUFHLEVBQUUsQ0FBQzthQUNuQjtZQUNELElBQUksQ0FBQyxN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TUFBTyxDQUFDLFVBQVUsRUFBRSxDQUFDLElBQUksQ0FBQztRQUNqRCxPQUFPLEVBQUUsTUFBTSxFQUFFLE1BQU0sRUFBRSxVQUFVLEVBQUUsU0FBUyxFQUFFLENBQUM7SUFDckQsQ0FBQztDQUFBO0FBRUQsTUFBTSxVQUFnQixVQUFVLENBQUMsSUFBUyxFQUFFLFNBQVMsR0FBRyxJQUFJOztRQUN4RCxNQUFNLENBQUMsSUFBSSxDQUFDLFNBQVMsQ0FBQyxJQUFJLENBQUMsTUFBTyxDQUFDLFVBQVUsRUFBRSxDQUFDLElBQUksQ0FBQyxFQUFFLElBQUksQ0FBQyxTQUFTLENBQUMsSUFBSSxDQUFDLFVBQVUsQ0FBQyxDQUFDLENBQUM7SUFDNUYsQ0FBQztDQUFBIiwic291cmNlc0NvbnRlbnQiOlsiaW1wb3J0IHR5cGUgKiBhcyBfZ3JvayBmcm9tICdkYXRhZ3Jvay1hcGkvZ3Jvayc7XHJcbmltcG9ydCB0eXBlICogYXMgX0RHIGZyb20gJ2RhdGFncm9rLWFwaS9kZyc7XHJcbmRlY2xhcmUgbGV0IGdyb2s6IHR5cGVvZiBfZ3JvaywgREc6IHR5cGVvZiBfREc7XHJcbmltcG9ydCB7IGRlbGF5LCBleHBlY3QsIHRlc3RFdmVudCwgdGVzdEV2ZW50QXN5bmMgfSBmcm9tIFwiLi90ZXN0XCI7XHJcbmltcG9ydCB7IE9ic2VydmFibGUgfSBmcm9tICdyeGpzJztcclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBjcmVhdGVWaWV3ZXIodHY6IF9ERy5UYWJsZVZpZXcsIHY6IHN0cmluZywgcGFja2FnZU5hbWU/OiBzdHJpbmcpOiBQcm9taXNlPF9ERy5WaWV3ZXI+IHtcclxuICAgIGxldCByZXM6IF9ERy5WaWV3ZXI7XHJcbiAgICBpZiAocGFja2FnZU5hbWUpIHtcclxuICAgICAgICByZXMgPSBhd2FpdCB0di5kYXRhRnJhbWUucGxvdC5mcm9tVHlwZSh2KSBhcyBfREcuVmlld2VyO1xyXG4gICAgICAgIHR2LmRvY2tNYW5hZ2VyLmRvY2socmVzKTtcclxuICAgIH0gZWxzZVxyXG4gICAgICAgIHJlcyA9IHR2LmFkZFZpZXdlcih2KTtcclxuICAgIHJldHVybiByZXM7XHJcbn1cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiB0ZXN0Vmlld2VySW50ZXJuYWwodHY6IF9ERy5UYWJsZVZpZXcsIHZpZXdlck5hbWU6IHN0cmluZywgcGFja2FnZU5hbWU6IHN0cmluZyxcclxuICAgIGV2ZW50OiBPYnNlcnZhYmxlPGFueT4sIGFjdGlvbnM/OiAoYXJnczogYW55LCB3aXRoRGVsYXk6IGJvb2xlYW4pID0+IFByb21pc2U8YW55PixcclxuICAgIGF3YWl0Vmlld2VyPzogKHZpZXdlcjogX0RHLlZpZXdlcikgPT4gUHJvbWlzZTx2b2lkPiwgbGF5b3V0PzogX0RHLlZpZXdMYXlvdXQsXHJcbiAgICBhY3Rpb25BcmdzPzogYW55LCBhY3Rpb25zV2l0aERlbGF5ID0gdHJ1ZSk6IFByb21pc2U8YW55PiB7XHJcbiAgICBsZXQgYWN0aW9uc1JlczogYW55ID0gbnVsbDtcclxuICAgIGF3YWl0IHRlc3RFdmVudEFzeW5jKGV2ZW50LCBhc3luYyAoZTogYW55KSA9PiB7XHJcbiAgICAgICAgbGV0IHZpZXdlcjogX0RHLlZpZXdlciB8IG51bGwgPSBudWxsO1xyXG4gICAgICAgIGZvciAoY29uc3QgdjEgb2YgdHYhLnZpZXdlcnMpIHtcclxuICAgICAgICAgICAgaWYgKHYxLnR5cGUgPT09IHZpZXdlck5hbWUpXHJcbiAgICAgICAgICAgICAgICB2aWV3ZXIgPSB2MTtcclxuICAgICAgICB9XHJcbiAgICAgICAgaWYgKCF2aWV3ZXIpXHJcbiAgICAgICAgICAgIHRocm93IEVycm9yKCdWaWV3ZXIgaGFzblxcJ3QgYmVlbiBhZGRlZCcpO1xyXG4gICAgICAgIGF3YWl0IFByb21pc2UucmVzb2x2ZSgpOyAvL3JlIHNjaGVkdWxlcyBzdWJzZXF1ZW50IGNvbW1hbmRzIGludG8gbWljcm90YXNrXHJcbiAgICAgICAgaWYgKGF3YWl0Vmlld2VyKSBhd2FpdCBhd2FpdFZpZXdlcih2aWV3ZXIpO1xyXG4gICAgICAgIGlmIChhY3Rpb25zKSB7XHJcbiAgICAgICAgICAgIGNvbnN0IGFyZ3MgPSBhY3Rpb25BcmdzID8/IHt9O1xyXG4gICAgICAgICAgICBhcmdzLnR2ID0gdHY7XHJcbiAgICAgICAgICAgIGFyZ3Mudmlld2VyID0gdmlld2VyO1xyXG4gICAgICAgICAgICBhY3Rpb25zUmVzID0gYXdhaXQgYWN0aW9ucyhhcmdzLCBhY3Rpb25zV2l0aERlbGF5KTtcclxuICAgICAgICB9XHJcbiAgICAgICAgIC8vY2hlY2sgdGhhdCB0aGVyZSBhcmUgbm8gYWN0aXZlIHN1YnNjcmlwdGlvbnMgaW4gdGhlIHZpZXdlciBhZnRlciBjbG9zZVxyXG4gICAgICAgICBhd2FpdCB0ZXN0RXZlbnQoZ3Jvay5ldmVudHMub25WaWV3ZXJDbG9zZWQsICgpID0+IHtcclxuICAgICAgICAgICAgIGV4cGVjdCh2aWV3ZXIhLnN1YnMuc29tZSgocykgPT4gIXMuY2xvc2VkKSwgZmFsc2UpO1xyXG4gICAgICAgICB9LCAoKSA9PiB2aWV3ZXIhLmNsb3NlKCksIDMwMDApO1xyXG4gICAgfSwgYXN5bmMgKCkgPT4ge1xyXG4gICAgICAgIGxheW91dCA/IHR2LmxvYWRMYXlvdXQobGF5b3V0KSA6IGF3YWl0IGNyZWF0ZVZpZXdlcih0diEsIHZpZXdlck5hbWUsIHBhY2thZ2VOYW1lKTtcclxuICAgIH0sIDYwMDAwLCAnVEVTVF9FVkVOVF9BU1lOQycpO1xyXG4gICAgaWYgKGFjdGlvbnNSZXMpXHJcbiAgICAgICAgcmV0dXJuIGFjdGlvbnNSZXM7XHJcbn1cclxuXHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gc2VsZWN0RmlsdGVyQ2hhbmdlQ3VycmVudChhcmdzOiBhbnksIHdpdGhEZWxheSA9IHRydWUpOiBQcm9taXNlPHZvaWQ+IHtcclxuICAgIGNvbnN0IGN1cnJlbnREZiA9IGFyZ3MudHYuZGF0YUZyYW1lO1xyXG4gICAgY29uc3QgZGZTYXZlZCA9IGN1cnJlbnREZi5jbG9uZSgpO1xyXG4gICAgLy9yZW1vdmUgdmFsdWVzIGluIHRoZSBmaXJzdCByb3dcclxuICAgIEFycmF5LmZyb20oY3VycmVudERmLnJvdygwKS5jZWxscykuZm9yRWFjaCgoYzogYW55KSA9PiBjLnZhbHVlID0gbnVsbCk7XHJcbiAgICAvL3NlbGVjdGlvblxyXG4gICAgY29uc3QgbnVtID0gY3VycmVudERmLnJvd0NvdW50IDwgMjAgPyBNYXRoLmZsb29yKGN1cnJlbnREZi5yb3dDb3VudCAvIDIpIDogMTA7XHJcbiAgICBjdXJyZW50RGYucm93cy5zZWxlY3QoKHJvdzogX0RHLlJvdykgPT4gcm93LmlkeCA+PSAwICYmIHJvdy5pZHggPCBudW0pO1xyXG4gICAgaWYgKHdpdGhEZWxheSlcclxuICAgICAgICBhd2FpdCBkZWxheSg1MCk7XHJcbiAgICAvL2ZpbHRlclxyXG4gICAgZm9yIChsZXQgaSA9IG51bTsgaSA8IG51bSAqIDI7IGkrKykgY3VycmVudERmLmZpbHRlci5zZXQoaSwgZmFsc2UpO1xyXG4gICAgaWYgKHdpdGhEZWxheSlcclxuICAgICAgICBhd2FpdCBkZWxheSg1MCk7XHJcbiAgICAvL2NoYW5nZSBjdXJyZW50IHJvd1xyXG4gICAgY3VycmVudERmLmN1cnJlbnRSb3dJZHggPSAxO1xyXG4gICAgLy9yZW1vdmUgY29sdW1uc1xyXG4gICAgY3VycmVudERmLmNvbHVtbnMubmFtZXMoKS5zbGljZSgwLCBNYXRoLmNlaWwoY3VycmVudERmLmNvbHVtbnMubGVuZ3RoIC8gMikpXHJcbiAgICAgICAgLmZvckVhY2goKGM6IGFueSkgPT4gY3VycmVudERmLmNvbHVtbnMucmVtb3ZlKGMpKTtcclxuICAgIGlmICh3aXRoRGVsYXkpXHJcbiAgICAgICAgYXdhaXQgZGVsYXkoMTAwKTtcclxuICAgIC8vc2V0IGJhY2sgaW5pdGlhbCBkZiB3aXRoIHdob2xlIHNldCBvZiBjb2x1bW5zIGFuZCBwcmVzZXJ2ZWQgZGF0YVxyXG4gICAgYXJncy50diEuZGF0YUZyYW1lID0gZGZTYXZlZDtcclxuICAgIGF3YWl0IGRlbGF5KDUwKTtcclxufVxyXG5cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBmaWx0ZXJBc3luYyhhcmdzOiBhbnksIHdpdGhEZWxheSA9IHRydWUpOiBQcm9taXNlPHZvaWQ+IHtcclxuICAgIGNvbnN0IGN1cnJlbnREZjogX0RHLkRhdGFGcmFtZSA9IGFyZ3MudHYuZGF0YUZyYW1lO1xyXG4gICAgc2V0VGltZW91dCgoKSA9PiBjdXJyZW50RGYuZmlsdGVyLnNldCgwLCAhY3VycmVudERmLmZpbHRlci5nZXQoMCkpLCAwKTtcclxufVxyXG5cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dChhcmdzOiBhbnksIHdpdGhEZWxheSA9IHRydWUpOiBQcm9taXNlPHsgbGF5b3V0OiBfREcuVmlld0xheW91dCwgc2F2ZWRQcm9wczogYW55IH0+IHtcclxuICAgIC8vZ2V0IGN1cnJlbnQgb3B0aW9ucyBhbmQgcHJvcGVydGllc1xyXG4gICAgbGV0IG9wdG5zOiB7IFtwOiBzdHJpbmddOiBhbnkgfTtcclxuICAgIHRyeSB7XHJcbiAgICAgICAgb3B0bnMgPSBhcmdzLnZpZXdlciEuZ2V0T3B0aW9ucyh0cnVlKS5sb29rO1xyXG4gICAgfSBjYXRjaCAoZXJyOiBhbnkpIHtcclxuICAgICAgICAvL0B0cy1pZ25vcmVcclxuICAgICAgICB0aHJvdyBuZXcgRXJyb3IoYFZpZXdlcidzIC5nZXRPcHRpb25zKCkgZXJyb3IuYCwgeyBjYXVzZTogZXJyIH0pO1xyXG4gICAgfVxyXG4gICAgbGV0IHByb3BzOiBfREcuUHJvcGVydHlbXTtcclxuICAgIHRyeSB7XHJcbiAgICAgICAgcHJvcHMgPSBhcmdzLnZpZXdlciEuZ2V0UHJvcGVydGllcygpO1xyXG4gICAgfSBjYXRjaCAoZXJyOiBhbnkpIHtcclxuICAgICAgICAvL0B0cy1pZ25vcmVcclxuICAgICAgICB0aHJvdyBuZXcgRXJyb3IoYFZpZXdlcidzIC5nZXRQcm9wZXJ0aWVzKCkgZXJyb3IuYCwgeyBjYXVzZTogZXJyIH0pO1xyXG4gICAgfVxyXG4gICAgLy9jaGFuZ2Ugb3B0aW9ucyBhbmQgcHJvcGVydGllc1xyXG4gICAgY29uc3QgbmV3UHJvcHM6IFJlY29yZDxzdHJpbmcsIHN0cmluZyB8IGJvb2xlYW4+ID0ge307XHJcbiAgICBPYmplY3Qua2V5cyhvcHRucykuZmlsdGVyKChrKSA9PiB0eXBlb2Ygb3B0bnNba10gPT09ICdib29sZWFuJykuZm9yRWFjaCgoaykgPT4gbmV3UHJvcHNba10gPSAhb3B0bnNba10pO1xyXG4gICAgcHJvcHMuZmlsdGVyKChwOiBfREcuUHJvcGVydHkpID0+IHAuY2hvaWNlcyAhPT0gbnVsbClcclxuICAgICAgICAuZm9yRWFjaCgocDogX0RHLlByb3BlcnR5KSA9PiBuZXdQcm9wc1twLm5hbWVdID0gcC5jaG9pY2VzLmZpbmQoKGM6IGFueSkgPT4gYyAhPT0gb3B0bnNbcC5uYW1lXSkhKTtcclxuICAgIC8vc2V0IG5ldyBvcHRpb25zXHJcbiAgICBhcmdzLnZpZXdlciEuc2V0T3B0aW9ucyhuZXdQcm9wcyk7XHJcbiAgICBhd2FpdCBkZWxheSgzMDApO1xyXG4gICAgY29uc3QgbGF5b3V0ID0gYXJncy50diEuc2F2ZUxheW91dCgpO1xyXG4gICAgY29uc3Qgc2F2ZWRMb29rID0gYXJncy52aWV3ZXIhLmdldE9wdGlvbnMoKS5sb29rO1xyXG4gICAgcmV0dXJuIHsgbGF5b3V0OiBsYXlvdXQsIHNhdmVkUHJvcHM6IHNhdmVkTG9vayB9O1xyXG59XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gbG9hZExheW91dChhcmdzOiBhbnksIHdpdGhEZWxheSA9IHRydWUpOiBQcm9taXNlPHZvaWQ+IHtcclxuICAgIGV4cGVjdChKU09OLnN0cmluZ2lmeShhcmdzLnZpZXdlciEuZ2V0T3B0aW9ucygpLmxvb2spLCBKU09OLnN0cmluZ2lmeShhcmdzLnNhdmVkUHJvcHMpKTtcclxufVxyXG4iXX0=