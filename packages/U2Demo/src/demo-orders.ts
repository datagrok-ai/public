/* The demo orders behind the `demoOrders` package function — a leaf module, so a page may use
   them without importing package.ts (which imports the pages: a cycle the test bundle cannot
   initialize — `SRC_ROOT` was read before nav.ts finished evaluating). */
import * as DG from 'datagrok-api/dg';

const ORDERS = [
  {orderId: 1001, customer: 'Aspirin Labs', city: 'Kyiv', total: 1240, daysAgo: 2},
  {orderId: 1002, customer: 'Bayer', city: 'Lviv', total: 380, daysAgo: 5},
  {orderId: 1003, customer: 'Roche', city: 'Basel', total: 2150, daysAgo: 11},
  {orderId: 1004, customer: 'Novartis', city: 'Basel', total: 640, daysAgo: 24},
  {orderId: 1005, customer: 'Pfizer', city: 'New York', total: 1790, daysAgo: 45},
  {orderId: 1006, customer: 'Merck', city: 'Darmstadt', total: 920, daysAgo: 88},
];

/** Demo orders placed within the last `days` days. */
export function demoOrders(days: number): DG.DataFrame {
  const rows = ORDERS.filter((order) => order.daysAgo <= days);
  return DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.create(0);
}
