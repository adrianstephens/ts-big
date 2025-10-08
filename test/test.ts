//import * as big from 'src/big';
import * as big from '../dist/big';
import * as dec from '../dist/dec';

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

console.log(dec.pi(100).toString());
console.log(dec.sin(1, 100).toString());

const a = dec.float.from(1);
const b = dec.float.from(10.12).addPrecision(100);

const c = a.add(b);
console.log(c.toString(10));
const d = a + b;
console.log(d.toString(10));
console.log((b * b).toString());
console.log(b.sqrt().toString());
console.log((b * 10).sqrt().toString());
console.log((b * 100).sqrt().toString());

//console.log(big.tan(inf, 1000).toString(10, Infinity));
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
