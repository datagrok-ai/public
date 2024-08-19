var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import * as DG from 'datagrok-api/dg';
import './test/demo-project-test';
import './test/projects-test';
import { runTests, tests } from '@datagrok-libraries/utils/src/test';
export const _package = new DG.Package();
export { tests };
//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export function test(category, test, testContext) {
    return __awaiter(this, void 0, void 0, function* () {
        const data = yield runTests({ category, test, testContext });
        return DG.DataFrame.fromObjects(data);
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoicGFja2FnZS10ZXN0LmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXMiOlsicGFja2FnZS10ZXN0LnRzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiI7Ozs7Ozs7OztBQUFBLE9BQU8sS0FBSyxFQUFFLE1BQU0saUJBQWlCLENBQUM7QUFDdEMsT0FBTywwQkFBMEIsQ0FBQztBQUNsQyxPQUFPLHNCQUFzQixDQUFDO0FBRTlCLE9BQU8sRUFBQyxRQUFRLEVBQUUsS0FBSyxFQUFjLE1BQU0sb0NBQW9DLENBQUM7QUFFaEYsTUFBTSxDQUFDLE1BQU0sUUFBUSxHQUFHLElBQUksRUFBRSxDQUFDLE9BQU8sRUFBRSxDQUFDO0FBQ3pDLE9BQU8sRUFBQyxLQUFLLEVBQUMsQ0FBQztBQUVmLFlBQVk7QUFDWix5Q0FBeUM7QUFDekMscUNBQXFDO0FBQ3JDLDRDQUE0QztBQUM1QywwQkFBMEI7QUFDMUIsTUFBTSxVQUFnQixJQUFJLENBQUMsUUFBZ0IsRUFBRSxJQUFZLEVBQUUsV0FBd0I7O1FBQ2pGLE1BQU0sSUFBSSxHQUFHLE1BQU0sUUFBUSxDQUFDLEVBQUMsUUFBUSxFQUFFLElBQUksRUFBRSxXQUFXLEVBQUMsQ0FBQyxDQUFDO1FBQzNELE9BQU8sRUFBRSxDQUFDLFNBQVMsQ0FBQyxXQUFXLENBQUMsSUFBSSxDQUFFLENBQUM7SUFDekMsQ0FBQztDQUFBIiwic291cmNlc0NvbnRlbnQiOlsiaW1wb3J0ICogYXMgREcgZnJvbSAnZGF0YWdyb2stYXBpL2RnJztcbmltcG9ydCAnLi90ZXN0L2RlbW8tcHJvamVjdC10ZXN0JztcbmltcG9ydCAnLi90ZXN0L3Byb2plY3RzLXRlc3QnO1xuXG5pbXBvcnQge3J1blRlc3RzLCB0ZXN0cywgVGVzdENvbnRleHR9IGZyb20gJ0BkYXRhZ3Jvay1saWJyYXJpZXMvdXRpbHMvc3JjL3Rlc3QnO1xuXG5leHBvcnQgY29uc3QgX3BhY2thZ2UgPSBuZXcgREcuUGFja2FnZSgpO1xuZXhwb3J0IHt0ZXN0c307XG5cbi8vbmFtZTogdGVzdFxuLy9pbnB1dDogc3RyaW5nIGNhdGVnb3J5IHtvcHRpb25hbDogdHJ1ZX1cbi8vaW5wdXQ6IHN0cmluZyB0ZXN0IHtvcHRpb25hbDogdHJ1ZX1cbi8vaW5wdXQ6IG9iamVjdCB0ZXN0Q29udGV4dCB7b3B0aW9uYWw6IHRydWV9XG4vL291dHB1dDogZGF0YWZyYW1lIHJlc3VsdFxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRlc3QoY2F0ZWdvcnk6IHN0cmluZywgdGVzdDogc3RyaW5nLCB0ZXN0Q29udGV4dDogVGVzdENvbnRleHQpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICBjb25zdCBkYXRhID0gYXdhaXQgcnVuVGVzdHMoe2NhdGVnb3J5LCB0ZXN0LCB0ZXN0Q29udGV4dH0pO1xuICByZXR1cm4gREcuRGF0YUZyYW1lLmZyb21PYmplY3RzKGRhdGEpITtcbn1cbiJdfQ==