import * as bits from '../src/bits';
import * as big from '../src/big';

export interface equal<T> {
	equal(b: T): boolean;
}

export function expect<T extends equal<T>>(v: T) {
	return {
		toEqual(v2: T) {
			if (!v.equal(v2))
				console.log("fail");
		}
	};
}

export function test(name: string, fn: ()=>void) {
	console.log("testing: " + name);
	fn();
	console.log("finished: " + name);
}

function timed(name: string, fn: ()=>void) {
	console.log("timed: " + name);
	const start = Date.now();
	fn();
	const end = Date.now();
	console.log("finished: " + name + " in " + (end - start) + "ms");
}

const inf = big.float.Infinity;
console.log(big.float.from(2).setPrecision(60).sqrt().toString(10, Infinity));

const a = big.float.from('0.1');
console.log(a.toString(10, Infinity));
console.log(a.addPrecision(100).toString(10, Infinity));
const b = big.float.fromString('0.1', 100);
console.log(b.toString(10, Infinity));
console.log(b.addPrecision(100).toString(10, Infinity));



console.log(big.tan(inf, 1000).toString(10, Infinity));
console.log(big.atan(big.float.Infinity, 1000).toString(10, Infinity));
console.log(big.tan(big.atan(big.float.fromString('0.1e2',1000), 1000), 1000).toString(10, Infinity));

const x = big.log(42, 1000);
const y = big.exp(x, 1000);
console.log(y.toString(10, Infinity));

console.log(big.pi(150).toString(10));
console.log(big.float.from(2).addPrecision(1000).sqrt().toString(10));

const t = big.float.from('42.1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890e1000');
console.log(t.toString(10));
console.log(t.sqrt().toString(10));
console.log(t.pow(3).root(3).toString(10));

console.log(t.sub(1).toString(10));
console.log(t.mul(t).toString(10));
console.log(t.div(10).toString(10));

//const sp = new bits.DenseBits();
const sp = new bits.SparseBits2();
sp.setRange(42, 100);
sp.clearRange(64, 96);
sp.set(42);
sp.set(1000);
sp.set(10000);

for (const i of sp) {
	console.log(i);
}


/*
sp.selfNot();

for (let i = sp.next(-1, false); i !== -1; i = sp.next(i, false)) {
	console.log(i);
}
*/

sp.set(0);
sp.set(2);
sp.selfComplement();
for (const i of sp.ranges()) {
	console.log(i);
}


for (const i of sp.where(false)) {
	console.log(i);
}


for (let i = sp.next(-1); i !== -1; i = sp.next(i)) {
	console.log(i);
}
