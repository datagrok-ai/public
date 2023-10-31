var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import { category, test } from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import { _package } from '../package-test';
category('Connections', () => {
    test('queriesTest', () => __awaiter(void 0, void 0, void 0, function* () {
        const queries = yield grok.dapi.queries.filter(`options.testExpectedRows != null and package.shortName = "${_package.name}"`).include('params').list();
        for (const query of queries) {
            const call = query.prepare();
            console.log("NAME: " + query.name);
            for (const property of query.inputs)
                property.set(call, property.defaultValue == null ? property.defaultValue : yield grok.functions.eval(`${property.defaultValue}`));
            yield call.call();
            const t = call.getOutputParamValue();
            if (t == null)
                throw 'Result of ' + query.name + 'is not a DataFrame';
            if (t.rowCount.toString() != query.options['testExpectedRows'])
                // eslint-disable-next-line no-throw-literal
                throw 'Rows number in' + query.name + 'table is not as expected';
        }
    }));
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiZGVtby1wcm9qZWN0LXRlc3QuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyJkZW1vLXByb2plY3QtdGVzdC50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7Ozs7Ozs7QUFBQSxPQUFPLEVBQUMsUUFBUSxFQUFFLElBQUksRUFBd0IsTUFBTSxvQ0FBb0MsQ0FBQztBQUN6RixPQUFPLEtBQUssSUFBSSxNQUFNLG1CQUFtQixDQUFDO0FBRTFDLE9BQU8sRUFBRSxRQUFRLEVBQUUsTUFBTSxpQkFBaUIsQ0FBQztBQUUzQyxRQUFRLENBQUMsYUFBYSxFQUFFLEdBQUcsRUFBRTtJQUMzQixJQUFJLENBQUMsYUFBYSxFQUFFLEdBQVMsRUFBRTtRQUM3QixNQUFNLE9BQU8sR0FBRyxNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLE1BQU0sQ0FBQyw2REFBNkQsUUFBUSxDQUFDLElBQUksR0FBRyxDQUFDLENBQUMsT0FBTyxDQUFDLFFBQVEsQ0FBQyxDQUFDLElBQUksRUFBRSxDQUFDO1FBRXZKLEtBQUssTUFBTSxLQUFLLElBQUksT0FBTyxFQUFFO1lBQzNCLE1BQU0sSUFBSSxHQUFHLEtBQUssQ0FBQyxPQUFPLEVBQUUsQ0FBQztZQUU3QixPQUFPLENBQUMsR0FBRyxDQUFDLFFBQVEsR0FBRyxLQUFLLENBQUMsSUFBSSxDQUFDLENBQUM7WUFFbkMsS0FBSyxNQUFNLFFBQVEsSUFBSSxLQUFLLENBQUMsTUFBTTtnQkFDakMsUUFBUSxDQUFDLEdBQUcsQ0FBQyxJQUFJLEVBQUUsUUFBUSxDQUFDLFlBQVksSUFBSSxJQUFJLENBQUMsQ0FBQyxDQUFDLFFBQVEsQ0FBQyxZQUFZLENBQUMsQ0FBQyxDQUFDLE1BQU0sSUFBSSxDQUFDLFNBQVMsQ0FBQyxJQUFJLENBQUMsR0FBRyxRQUFRLENBQUMsWUFBWSxFQUFFLENBQUMsQ0FBQyxDQUFDO1lBRXBJLE1BQU0sSUFBSSxDQUFDLElBQUksRUFBRSxDQUFDO1lBRWxCLE1BQU0sQ0FBQyxHQUFHLElBQUksQ0FBQyxtQkFBbUIsRUFBa0IsQ0FBQztZQUNyRCxJQUFJLENBQUMsSUFBSSxJQUFJO2dCQUNYLE1BQU0sWUFBWSxHQUFHLEtBQUssQ0FBQyxJQUFJLEdBQUcsb0JBQW9CLENBQUM7WUFFekQsSUFBSSxDQUFDLENBQUMsUUFBUSxDQUFDLFFBQVEsRUFBRSxJQUFJLEtBQUssQ0FBQyxPQUFPLENBQUMsa0JBQWtCLENBQUM7Z0JBQzVELDRDQUE0QztnQkFDNUMsTUFBTSxnQkFBZ0IsR0FBRyxLQUFLLENBQUMsSUFBSSxHQUFHLDBCQUEwQixDQUFDO1NBQ3BFO0lBQ0gsQ0FBQyxDQUFBLENBQUMsQ0FBQztBQUNMLENBQUMsQ0FBQyxDQUFDIiwic291cmNlc0NvbnRlbnQiOlsiaW1wb3J0IHtjYXRlZ29yeSwgdGVzdCwgYmVmb3JlLCBhZnRlciwgZXhwZWN0fSBmcm9tICdAZGF0YWdyb2stbGlicmFyaWVzL3V0aWxzL3NyYy90ZXN0JztcbmltcG9ydCAqIGFzIGdyb2sgZnJvbSAnZGF0YWdyb2stYXBpL2dyb2snO1xuaW1wb3J0ICogYXMgREcgZnJvbSAnZGF0YWdyb2stYXBpL2RnJztcbmltcG9ydCB7IF9wYWNrYWdlIH0gZnJvbSAnLi4vcGFja2FnZS10ZXN0JztcblxuY2F0ZWdvcnkoJ0Nvbm5lY3Rpb25zJywgKCkgPT4ge1xuICB0ZXN0KCdxdWVyaWVzVGVzdCcsIGFzeW5jICgpID0+IHtcbiAgICBjb25zdCBxdWVyaWVzID0gYXdhaXQgZ3Jvay5kYXBpLnF1ZXJpZXMuZmlsdGVyKGBvcHRpb25zLnRlc3RFeHBlY3RlZFJvd3MgIT0gbnVsbCBhbmQgcGFja2FnZS5zaG9ydE5hbWUgPSBcIiR7X3BhY2thZ2UubmFtZX1cImApLmluY2x1ZGUoJ3BhcmFtcycpLmxpc3QoKTtcbiAgIFxuICAgIGZvciAoY29uc3QgcXVlcnkgb2YgcXVlcmllcykge1xuICAgICAgY29uc3QgY2FsbCA9IHF1ZXJ5LnByZXBhcmUoKTtcblxuICAgICAgY29uc29sZS5sb2coXCJOQU1FOiBcIiArIHF1ZXJ5Lm5hbWUpOyAgICAgIFxuXG4gICAgICBmb3IgKGNvbnN0IHByb3BlcnR5IG9mIHF1ZXJ5LmlucHV0cylcbiAgICAgICAgcHJvcGVydHkuc2V0KGNhbGwsIHByb3BlcnR5LmRlZmF1bHRWYWx1ZSA9PSBudWxsID8gcHJvcGVydHkuZGVmYXVsdFZhbHVlIDogYXdhaXQgZ3Jvay5mdW5jdGlvbnMuZXZhbChgJHtwcm9wZXJ0eS5kZWZhdWx0VmFsdWV9YCkpO1xuXG4gICAgICBhd2FpdCBjYWxsLmNhbGwoKTtcblxuICAgICAgY29uc3QgdCA9IGNhbGwuZ2V0T3V0cHV0UGFyYW1WYWx1ZSgpIGFzIERHLkRhdGFGcmFtZTtcbiAgICAgIGlmICh0ID09IG51bGwpXG4gICAgICAgIHRocm93ICdSZXN1bHQgb2YgJyArIHF1ZXJ5Lm5hbWUgKyAnaXMgbm90IGEgRGF0YUZyYW1lJztcblxuICAgICAgaWYgKHQucm93Q291bnQudG9TdHJpbmcoKSAhPSBxdWVyeS5vcHRpb25zWyd0ZXN0RXhwZWN0ZWRSb3dzJ10pXG4gICAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBuby10aHJvdy1saXRlcmFsXG4gICAgICAgIHRocm93ICdSb3dzIG51bWJlciBpbicgKyBxdWVyeS5uYW1lICsgJ3RhYmxlIGlzIG5vdCBhcyBleHBlY3RlZCc7XG4gICAgfVxuICB9KTtcbn0pO1xuIl19