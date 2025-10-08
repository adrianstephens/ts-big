import { bits } from '@isopodlabs/utilities';
import { sqrt, root, compare, randomBits, Round, RoundMode } from './big';

function round(m: bigint, n: number, mode: RoundMode): bigint {
	const	neg = m < 0n;
	let		a	= neg ? -m : m;
	const	r	= 10n ** BigInt(n);

	switch (mode) {
		case Round.trunc:
			break;

		case Round.down:
			if (neg && (a % r))
				a += r;
			break;

		case Round.up:
			if (!neg && (a % r))
				a += r;
			break;

		case Round.halfUp:
			a += r >> 1n;
			break;

		case Round.halfEven: {
			const f = a % r;
			if ((f >= r / 2n) && (a % (r << 1n)))
				a += r;
			break;
		}

	}

	a /= r;
	return neg ? -a : a;
}

export class float {
	static readonly zero	= new float(0, 0n);
	static readonly one		= new float(0, 1n);
	static readonly two		= new float(0, 2n);
	static readonly Infinity= new float(Infinity, 1n);

	static from(v: number|bigint|string|float): float {
		switch (typeof v) {
			case 'object':
				return v;

			case 'bigint':
				return new float(0, v);

			case 'number':
				if (v === Infinity)
					return this.Infinity;
				if (v === -Infinity)
					return this.Infinity.neg();
				if (v === Math.floor(v))
					return new float(0, BigInt(v));
				v = v.toString();
				//fallthrough

			case 'string': {
				const m = v.match(/^(-?)(\d+)(?:\.(\d*))?(?:e([-+]?\d+))?$/);
				if (m) {
					let x	= BigInt(m[2]);
					let	e	= (m[4] ? parseInt(m[4]) : 0);

					if (m[3]) {
						const d = m[3].length;
						x = x * 10n ** BigInt(d) + BigInt(m[3]);
						e -= d;
					}
					return new float(e, x);
				}
				return new float(0, BigInt(v));
			}
		}
	}

	constructor(public exponent: number, public mantissa: bigint) {
	}

	private _rep(exponent: number): bigint {//use a's exponent
		return this.mantissa * 10n ** BigInt(this.exponent - exponent);
	}

	addPrecision(p: number, mode: RoundMode = Round.halfEven): float {
		if (p === Infinity)
			return this.mantissa < 0 ? float.Infinity.neg() : float.Infinity;
		return new float(this.exponent - p, p < 0 ? round(this.mantissa, -p, mode) : this.mantissa * 10n ** BigInt(p));
	}
	setPrecision(p: number, mode: RoundMode = Round.halfEven): float {
		return this.addPrecision(this.exponent + p, mode);
	}
	capPrecision(p: number, mode: RoundMode = Round.halfEven): float {
		if (this.exponent + p < 0)
			return this.addPrecision(this.exponent + p, mode);
		return this;
	}
	neg(): float {
		return new float(this.exponent, -this.mantissa);
	}
	abs(): float {
		return new float(this.exponent, this.mantissa < 0n ? -this.mantissa : this.mantissa);
	}

	private static _add(a: float, b: float): float {//use a's exponent
		if (b.exponent === Infinity)
			return b;
		return new float(a.exponent, a.mantissa + (b._rep(a.exponent)));
	}

	add(other: float|number|bigint): float {
		other = float.from(other);
		return this.exponent < other.exponent
			? float._add(this, other)
			: float._add(other, this);
	}
	sub(other: float|number|bigint): float {
		other = float.from(other).neg();
		return this.exponent < other.exponent
			? float._add(this, other)
			: float._add(other, this);
	}
	mul(other: float|number|bigint): float {
		other = float.from(other);
		return new float(this.exponent + other.exponent, this.mantissa * other.mantissa);
	}
	div(other: float|number|bigint): float {
		other = float.from(other);
		return other.exponent === Infinity
			? float.zero
			: other.mantissa === 0n
			? float.Infinity
			: new float(this.exponent, this._rep(this.exponent + other.exponent) / other.mantissa);
	}
	mod(other: float|number|bigint): float {
		other = float.from(other);
		return this.lt(other)
			? this
			: this.exponent < other.exponent
			? new float(this.exponent, this.mantissa % (other.mantissa << BigInt(other.exponent - this.exponent)))
			: this.exponent === Infinity
			? float.zero
			: new float(other.exponent, (this.mantissa << BigInt(this.exponent - other.exponent)) % other.mantissa);
	}
	divmod(other: float|number|bigint) {
		other = float.from(other);
		const e = Math.min(this.exponent, other.exponent);
		const m0 = this.mantissa << BigInt(this.exponent - e);
		const m1 = other.mantissa << BigInt(other.exponent - e);
		return [m0 / m1, new float(e, m0 % m1)];
	}
	square(): float {
		return new float(this.exponent + this.exponent, this.mantissa * this.mantissa);
	}
	sqrt(): float {
		return this.exponent === Infinity
			? float.Infinity
			: new float(this.exponent >> 1, sqrt(this.exponent & 1 ? this.mantissa * 10n : this.mantissa));
	}
	pow(other: number): float {
		return other === 0
			? float.one
			: other === 1
			? this
			: new float(this.exponent * other, this.mantissa ** BigInt(other));
	}
	root(base: number): float {
		if (base < 1)
			return new float(0, 0n);
		const e1 = -this.exponent % base;
		if (e1)
			return new float(((this.exponent + e1) / base) - 1, root(this._rep(this.exponent + e1), base));
		return new float(this.exponent / base, root(this.mantissa, base));
	}

	toInt(mode: RoundMode = Round.trunc): bigint {
		return this.exponent < 0
			? round(this.mantissa, -this.exponent, mode)
			: this.mantissa << BigInt(this.exponent);
	}

	frac(): float {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			return new float(this.exponent, this.mantissa % mask);
		}
		return float.zero;
	}
	floor(): float {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			if (this.mantissa % mask)
				return new float(this.exponent, (this.mantissa < 0 ? this.mantissa - mask : this.mantissa) / mask * mask);
		}
		return this;
	}
	ceil(): float {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			if (this.mantissa % mask)
				return new float(this.exponent, (this.mantissa < 0 ? this.mantissa : this.mantissa - mask) / mask * mask);
		}
		return this;
	}
	trunc(): float {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			return new float(this.exponent, this.mantissa - (this.mantissa % mask));
		}
		return this;
	}
	round(): float {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			const m = this.mantissa + (mask >> 1n);
			return new float(this.exponent, m / mask * mask);
		}
		return this;
	}

	compare(other: float|number|bigint): number {
		other = float.from(other);
		return this.exponent === other.exponent
			? compare(this.mantissa, other.mantissa)
			: this.exponent < other.exponent ? (other.exponent === Infinity ? -1 : compare(this.mantissa, other._rep(this.exponent)))
			: this.exponent === Infinity ? 1 : compare(this._rep(other.exponent), other.mantissa);
	}

	lt0(): boolean {
		return this.mantissa < 0n;
	}
	le(other: float|number|bigint): boolean {
		other = float.from(other);
		return this.exponent === other.exponent	? this.mantissa <= other.mantissa
			: this.exponent < other.exponent	? other.exponent === Infinity || this.mantissa <= (other._rep(this.exponent))
			: this.exponent !== Infinity && this._rep(other.exponent) <= other.mantissa;
	}
	eq(other: float|number|bigint): boolean {
		other = float.from(other);
		return this.exponent === other.exponent	? this.mantissa === other.mantissa
			: this.exponent < other.exponent	? other.exponent !== Infinity && this.mantissa === (other._rep(this.exponent))
			: this.exponent !== Infinity && this._rep(other.exponent) === other.mantissa;
	}
	ne(other: float|number|bigint): boolean {
		return !this.eq(other);
	}
	ge(other: float|number|bigint): boolean {
		return float.from(other).le(this);
	}
	lt(other: float|number|bigint): boolean {
		return !this.ge(other);
	}
	gt(other: float|number|bigint): boolean {
		return !this.le(other);
	}

	toString(base = 10): string {
		if (!this.mantissa)
			return '0';
		
		const	e = this.exponent;
		const	s = this.mantissa < 0n ? '-' : '';

		if (e > 0) {
			if (e === Infinity)
				return s + 'Infinity';
			return s + this._rep(0).toString(base);
		}

		const p = 10n ** BigInt(-e);
		const m = this.mantissa < 0n ? -this.mantissa : this.mantissa;
		const x = m / p, y = m % p;
		if (!y)
			return s + x.toString(base);

		return s + x.toString(base) + '.' + (p + y).toString(base).slice(1);
	}
}

export function max(...values: float[]) {
	if (values.length === 0)
		return float.Infinity.neg();
	let maxv = values[0];
	for (let i = 1; i < values.length; i++) {
		if (values[i].gt(maxv))
			maxv = values[i];
	}
	return maxv;
}

export function min(...values: float[]) {
	if (values.length === 0)
		return float.Infinity;
	let minv = values[0];
	for (let i = 1; i < values.length; i++) {
		if (values[i].lt(minv))
			minv = values[i];
	}
	return minv;
}

export function random(digits: number) {
	const m = randomBits(Math.ceil(digits * 3.321928094887362)); // = digits * log2(10)
	return new float(-digits, m % (10n ** BigInt(digits)));
}

//Gauss-Legendre iterative algorithm for pi
export function pi(digits: number) {
	let a = float.one.addPrecision(digits);
	let b = a.div(float.two.addPrecision(digits << 1).sqrt());
	let t = float.one;
	const bits = Math.ceil((digits + 4) * 3.321928094887362); // = (digits + 4) * log2(10)

	for (let p = 0; (8 << p) < bits; ++p) {
		t = t.sub(b.sub(a).square().mul(1 << p)).capPrecision(digits + 4);

		const ab = a.mul(b);
		a = a.add(b).div(float.two);
		b = ab.sqrt();
	}

	return a.add(b).square().div(t).setPrecision(digits);
}


function sin_helper(x: float, pi: float, digits: number): float {
	// Reduce to [-pi, pi]
	const twoPi = pi.mul(float.two);
	x = x.mod(twoPi);
	if (x.lt0())
		x = x.add(twoPi);

	// Further reduce to [-pi/2, pi/2]
	if (x.gt(pi.div(float.two)))
		x = pi.sub(x);

	// Taylor series
	let result	= x;
	const x2	= x.mul(x).neg();
	const limit = new float(-digits - 1, 1n);

	for (let n = 2n, term = x.setPrecision(digits + 1); ; n += 2n) {
		term = term.mul(x2).div(n * (n + 1n)).capPrecision(digits + 1);
		if (term.abs().lt(limit))
			break;
		result	= result.add(term);
	}
	return result;
}

export function sin(x: float|number|bigint, digits: number): float {
	const π = pi(digits + 2);
	return sin_helper(float.from(x), π, digits).setPrecision(digits);
}

export function cos(x: float|number|bigint, digits: number): float {
	const π = pi(digits + 2);
	return sin_helper(float.from(x).add(π.div(float.two)), π, digits).setPrecision(digits);
}

export function tan(x: float|number|bigint, digits: number): float {
	x = float.from(x);
	const π = pi(digits + 2);
	return sin_helper(x, π, digits).div(sin_helper(x.add(π.div(float.two)), π, digits)).setPrecision(digits);
}

function asin_helper(x: float, digits: number): float {
	const limit = new float(-digits - 1, 1n);
	const x2	= x.mul(x);
	let result	= x;

	for (let n = 2n, num = x, den = 1n; ; n += 2n) {
		num = num.mul(x2).capPrecision(digits + 1).mul((n - 1n) * (n - 1n));
		den *= n * (n + 1n);
		const term = num.div(den);
		if (term.abs().lt(limit))
			break;
		result = result.add(term);
	}
	return result;
}

// High-precision arcsin using Taylor series and argument reduction
export function asin(x: float|number|bigint, digits: number): float {
	x = float.from(x);

	if (x.lt(0.707))
		return asin_helper(x, digits).setPrecision(digits);
	// arcsin(x) = pi/2 - 2*arcsin(sqrt((1-x)/2))
	const π		= pi(digits + 2);
	const asin1	= asin_helper(float.one.sub(x).div(float.two).sqrt(), digits);
	const res	= π.div(float.two).sub(asin1.mul(float.two)).setPrecision(digits);
	return x.lt0() ? res.neg() : res;
}

// Use identity arccos(x) = pi/2 - arcsin(x)
export function acos(x: float|number|bigint, digits: number): float {
	const π		= pi(digits + 2);
	x = float.from(x);

	if (x.lt(0.707))
		return π.div(float.two).sub(asin_helper(x, digits)).setPrecision(digits);

	const asin	= asin_helper(float.one.sub(x).div(float.two).sqrt(), digits);
	return x.lt0()
		? π.sub(asin.mul(float.two)).setPrecision(digits)
		: asin.mul(float.two).setPrecision(digits);
}


// High-precision arctan using Taylor series and argument reduction
export function atan(x: float|number|bigint|string, digits: number): float {
	x = float.from(x);
	const limit = new float(-digits - 1, 1n);

	function helper(x: float): float {
		// Taylor series for |x| <= 1
		let result	= x;
		const x2 = x.mul(x).neg();

		for (let n = 3n, term = x; ; n += 2n) {
			term = term.mul(x2).capPrecision(digits + 1);
			const next = term.div(n);
			if (next.abs().lt(limit))
				break;
			result = result.add(next);
		}
		return result;
	}

	if (x.abs().le(float.one))
		return helper(x).setPrecision(digits);

	// Argument reduction for |x| > 1
	const π		= pi(digits + 2);
	const atan1	= helper(float.one.addPrecision(digits + 4).div(x));
	return (x.lt0() ? π.neg() : π).div(float.two).sub(atan1).setPrecision(digits);
}

// with x in [1,2), y = (x-1)/(x+1)
function log_helper(y: float, bits: number): float {
	const limit = new float(-bits - 4, 1n);
	const y2	= y.mul(y);

	let result	= y;
	for (let n = 3n, term = y; ; n += 2n) {
		term = term.mul(y2).capPrecision(bits + 4);
		const next = term.div(n);
		if (next.abs().lt(limit))
			break;
		result = result.add(next);
	}
	return result.mul(float.two);
}

const ln2s: float[] = [];
function ln2(digits: number): float {
	const nextpow2 = bits.highestSet(digits);
	return (ln2s[nextpow2] ??= log_helper(float.one.addPrecision(1 << nextpow2).div(3n), 1 << nextpow2)).setPrecision(digits);
}

// High-precision natural logarithm (ln) and exponential (exp)
export function log(x: float|number|bigint, digits: number): float {
	x = float.from(x);
	if (x.le(0n)) {
		if (x.lt0())
			throw new Error('ln(x): x must be positive');
		return float.Infinity.neg();
	}

	// Argument reduction: x = m * 2^k, ln(x) = ln(m) + k*ln(2)
	const k = bits.highestSet(x.mantissa) + x.exponent - 1;
	const m = x.div(1 << k);

	const result = log_helper(m.sub(float.one).addPrecision(digits).div(m.add(float.one)), digits);

	// ln(x) = ln(m) + k*ln(2)
	return result.add(ln2(digits).mul(k)).setPrecision(digits);
}

export function exp(x: float|number|bigint, digits: number): float {
	const [shift, xred] = float.from(x).divmod(ln2(digits));

	// Taylor series for exp(xred)
	const limit	= new float(-digits - 4, 1n);
	let result	= float.one;
	for (let n = 1n, term = float.one; ; n++) {
		term = term.mul(xred).div(n).capPrecision(digits + 4);
		if (term.abs().lt(limit))
			break;
		result = result.add(term);
	}
	
	// Undo argument reduction
	return result.setPrecision(digits).mul(1 << Number(shift));
}
