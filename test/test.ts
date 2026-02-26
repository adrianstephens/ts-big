import big from '../dist';
import {dec} from '../dist';

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
/*
const f2 = big.from(1.01);
for (let i = 0, t = f2; i < 100; i++) {
	console.log(t.toString());
	t = t.mul(f2);
}
*/

console.log(big.pi(100).toString());
console.log(dec.pi(100).toString());
console.log(dec.sin(1, 100).toString());

const a = dec.from(1e100);
const b = dec.from(10.12).addPrecision(100);

/*
const c = a.add(b);
console.log(c.toString(10));
const d = a + b;
console.log(d.toString(10));
console.log((b * b).toString());
console.log(b.sqrt().toString());
console.log((b * 10).sqrt().toString());
console.log((b * 100).sqrt().toString());
*/

//console.log(big.tan(inf, 1000).toString(10, Infinity));
console.log(big.atan(big.Infinity, 1000).toString(10, Infinity));
console.log(big.tan(big.atan(big.from('0.1e2',1000), 1000), 1000).toString(10, Infinity));

const x = dec.log(42, 1000);
const y = dec.exp(x, 1000);
console.log(y.toString(10));

console.log(big.pi(150).toString(10));
console.log(big.from(2).addPrecision(1000).sqrt().toString(10));

const t = big.from('42.1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890');
console.log(t.toString(10));
console.log(t.sqrt().toString(10));
console.log(t.pow(3).rpow(1, 3).toString(10));

console.log(t.sub(1).toString(10));
console.log(t.mul(t).toString(10));
console.log(t.div(10).toString(10));

console.log(t.pow(-3).toString());
console.log(t.pow(big.from('-3')).toString());

