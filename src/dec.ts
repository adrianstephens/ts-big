import { bits } from '@isopodlabs/utilities';
import { sqrt, root, compare, randomBits, Round, RoundMode } from './big';
const log2_10 = 3.321928094887362; // log2(10)

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

//-----------------------------------------------------------------------------
// float class
//-----------------------------------------------------------------------------

export class float {
	static readonly zero	= new float(0, 0n);
	static readonly one		= new float(0, 1n);
	static readonly two		= new float(0, 2n);
	static readonly Infinity= new float(Infinity, 1n);


	static max(...values: float[]) {
		let maxv = float.Infinity.neg();
		for (const i of values) {
			if (i.gt(maxv))
				maxv = i;
		}
		return maxv;
	}

	static min(...values: float[]) {
		let minv = float.Infinity;
		for (const i of values) {
			if (i.lt(minv))
				minv = i;
		}
		return minv;
	}

	static random(digits: number) {
		const m = randomBits(Math.ceil(digits * log2_10)); // = digits * log2(10)
		return new float(-digits, m % (10n ** BigInt(digits)));
	}

	static pi(digits: number): float {
		const nextpow2 = bits.highestSet(digits);
		return (pis[nextpow2] ??= pi_helper(1 << nextpow2)).setPrecision(digits);
	}

	static sin(x: float|number|bigint, digits: number): float {
		const π = this.pi(digits + 2);
		return sin_helper(float.from(x), π, digits).setPrecision(digits);
	}

	static cos(x: float|number|bigint, digits: number): float {
		const π = this.pi(digits + 2);
		return sin_helper(float.from(x).add(π.shift(-1)), π, digits).setPrecision(digits);
	}

	static tan(x: float|number|bigint, digits: number): float {
		x = float.from(x);
		const π = this.pi(digits + 2);
		return sin_helper(x, π, digits).div(sin_helper(x.add(π.shift(-1)), π, digits)).setPrecision(digits);
	}
	static asin(x: float|number|bigint, digits: number): float {
		x = float.from(x);

		if (x.lt(0.707))
			return asin_helper(x, digits).setPrecision(digits);
		
		// arcsin(x) = pi/2 - 2*arcsin(sqrt((1-x)/2))
		const π		= this.pi(digits + 2);
		const asin1	= asin_helper(float.one.sub(x).shift(-1).sqrt(), digits);
		const res	= π.shift(-1).sub(asin1.shift(1)).setPrecision(digits);
		return x.sign() < 0 ? res.neg() : res;
	}

	// Use identity arccos(x) = pi/2 - arcsin(x)
	static acos(x: float|number|bigint, digits: number): float {
		const π		= this.pi(digits + 2);
		x = float.from(x);
		if (x.lt(0.707))
			return π.shift(-1).sub(asin_helper(x, digits)).setPrecision(digits);

		const asin	= asin_helper(float.one.sub(x).shift(-1).sqrt(), digits);
		return x.sign() < 0
			? π.sub(asin.shift(1)).setPrecision(digits)
			: asin.shift(1).setPrecision(digits);
	}
	static atan(x: float|number|bigint|string, digits: number): float {
		x = float.from(x);

		if (x.abs().le(float.one))
			return atan_helper(x, digits).setPrecision(digits);

		// Argument reduction for |x| > 1
		const π_2	= this.pi(digits + 2).shift(-1);
		const a		= atan_helper(x.recip(), digits);
		return (x.sign() < 0 ? π_2.neg() : π_2).sub(a).setPrecision(digits);
	}

	static atan2(y: float|number|bigint|string, x: float|number|bigint|string, digits: number): float {
		y = float.from(y);
		x = float.from(x);
		const π_2 = this.pi(digits + 2).shift(-1);

		if (x.mantissa === 0n) {
			switch (y.sign()) {
				case 1:		return π_2;
				case -1:	return π_2.neg();
				default:	return float.zero;
			}
		}

		if (x.abs().gt(y.abs())) {
			const a = atan_helper(y.div(x), digits);
			return x.sign() > 0 ? a : a.add(y.sign() >= 0 ? π_2 : π_2.neg());
		} else {
			const a = atan_helper(x.div(y), digits);
			return y.sign() > 0 ? π_2.sub(a) : π_2.neg().sub(a);
		}
	}
	static log(x: float|number|bigint, nbits: number): float {
		return float.from(x).setPrecision(nbits).log();
	}

	static exp(x: float|number|bigint, bits: number): float {
		return float.from(x).setPrecision(Math.max(bits, 32)).exp();
	}

	static from<C extends new (exponent: number, mantissa: bigint) => any>(this: C, v: number|bigint|string|float): InstanceType<C> {
		switch (typeof v) {
			case 'object':
				return v instanceof this ? v as InstanceType<C> : new this(v.exponent, v.mantissa);

			case 'bigint':
				return new this(0, v);

			case 'number':
				if (v === Infinity)
					return new this(Infinity, 1n);
				if (v === -Infinity)
					return new this(Infinity, -1n);
				//if (Number.isInteger(v))
				//	return new float(0, BigInt(v));
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
					return new this(e, x);
				}
				return new this(0, BigInt(v));
			}
		}
	}

	log(): this {
		if (this.sign() <= 0) {
			if (this.sign() < 0)
				throw new Error('ln(x): x must be positive');
			return this.infinity().neg();
		}

		// Argument reduction: x = m * 2^k, ln(x) = ln(m) + k*ln(2)
		const k = bits.highestSet(this.mantissa) + Math.ceil(this.exponent * log2_10) - 1;
		const m = this.shift(-k);

		const digits = Math.max(-this.exponent, 10);
		const result = log_helper(m.sub(float.one).addPrecision(digits).div(m.add(float.one)), digits);

		// ln(x) = ln(m) + k*ln(2)
		return result.add(ln2(digits).mul(k)).setPrecision(digits);
	}

	exp(): this {
		if (this.exponent === Infinity && this.mantissa < 0n)
			return this.zero();

		const x		= this.create(this.exponent, this.mantissa);
		const digits = Math.max(-x.exponent, 10);
		const shift	= float.from(x).divmod(ln2(digits));

		// Taylor series for exp(xred)
		const limit	= new float(-digits - 4, 1n);
		let result	= this.one();
		for (let n = 1n, term = this.one(); ; n++) {
			term = term.mul(x).div(n).capPrecision(digits + 4);
			if (term.abs().lt(limit))
				break;
			result = result.add(term);
		}
		
		// Undo argument reduction
		return result.setPrecision(digits).shift(Number(shift));
	}

	constructor(public exponent: number, public mantissa: bigint) {
	}
	protected create(exponent: number, mantissa: bigint): this {
		return new (this.constructor as new (exponent: number, mantissa: bigint) => this)(exponent, mantissa);
	}
	from(other: number|bigint|string|float): this {
		return float.from.call(this.constructor as new (exponent: number, mantissa: bigint) => this, other) as this;
	}
	protected zero()		{ return this.create(0, 0n); }
	protected one()			{ return this.create(0, 1n); }
	protected infinity()	{ return this.create(Infinity, 1n); }

	//get mantissa with specified exponent
	private _rep(exponent: number): bigint {
		return this.mantissa * 10n ** BigInt(this.exponent - exponent);
	}

	addPrecision(p: number, mode: RoundMode = Round.halfEven): this {
		if (p === Infinity)
			return this.create(Infinity, this.mantissa < 0 ? -1n : 1n);
		return this.create(this.exponent - p, p < 0 ? round(this.mantissa, -p, mode) : this.mantissa * 10n ** BigInt(p));
	}
	setPrecision(p: number, mode: RoundMode = Round.halfEven): this {
		return this.addPrecision(this.exponent + p, mode);
	}
	capPrecision(p: number, mode: RoundMode = Round.halfEven): this {
		if (this.exponent + p < 0)
			return this.addPrecision(this.exponent + p, mode);
		return this;
	}
	shift(p: number): this {
		if (p > 0)
			return this.create(this.exponent, this.mantissa << BigInt(p));

		const lowest = bits.lowestSet(this.mantissa);
		return lowest < -p
			?	this.create(this.exponent - (-p - lowest), (this.mantissa >> BigInt(lowest)) * 5n ** BigInt(-p - lowest))
			:	this.create(this.exponent, this.mantissa >> BigInt(-p));
	}
	dup(): this {
		return this.create(this.exponent, this.mantissa);
	}
	neg(): this {
		return this.create(this.exponent, -this.mantissa);
	}
	abs(): this {
		return this.create(this.exponent, this.mantissa < 0n ? -this.mantissa : this.mantissa);
	}
	sign(): number {
		return this.mantissa === 0n ? 0 : this.mantissa < 0n ? -1 : 1;
	}
	recip(): this {
		return this.one().div(this);
	}

	add(other: float|number|bigint): this {
		other = float.from(other);
		return this.exponent < other.exponent
			? (other.exponent === Infinity ? this.infinity() : this.create(this.exponent, this.mantissa + other._rep(this.exponent)))
			: (this.exponent === Infinity ? this : this.create(other.exponent, this._rep(other.exponent) + other.mantissa));
	}
	sub(other: float|number|bigint): this {
		return this.add(float.from(other).neg());
	}
	scale(other: number) {
		return this.mul(other);
	}
	mul(other: float|number|bigint): this {
		other = float.from(other);
		return this.create(this.exponent + other.exponent, this.mantissa * other.mantissa);
	}
	div(other: float|number|bigint): this {
		other = float.from(other);
		return other.exponent === Infinity		? this.zero()
			: other.mantissa === 0n				? this.infinity()
			: this.create(this.exponent, this._rep(this.exponent + other.exponent) / other.mantissa);
	}
	mod(other: float|number|bigint): this {
		other = float.from(other);
		return this.lt(other)					? this
			: this.exponent < other.exponent	? this.create(this.exponent, this.mantissa % (other.mantissa << BigInt(other.exponent - this.exponent)))
			: this.exponent === Infinity		? this.zero()
			: this.create(other.exponent, (this.mantissa << BigInt(this.exponent - other.exponent)) % other.mantissa);
	}
	divmod(other: float|number|bigint) {
		other = float.from(other);
		if (other.mantissa === 0n)
			return Infinity;
		const e		= Math.min(this.exponent, other.exponent);
		const m0	= this._rep(e);
		const m1	= other._rep(e);
		this.exponent = e;
		this.mantissa = m0 % m1;
		return m0 / m1;
	}
	square(): this {
		return this.create(this.exponent + this.exponent, this.mantissa * this.mantissa);
	}
	sqrt(): this {
		return this.exponent === Infinity		? this.infinity()
			: this.create(this.exponent >> 1, sqrt(this.exponent & 1 ? this.mantissa * 10n : this.mantissa));
	}
	ipow(other: number): this {
		other = Math.floor(other);
		return other === 0		? this.one()
			: other === 1		? this
			: other < 0 		? this.one().div(this.ipow(-other))
			: this.create(this.exponent * other, this.mantissa ** BigInt(other));
	}
	rpow(n: number, d: number): this {
		const recip = n * d < 0;
		d = Math.floor(Math.abs(d));
		if (d === 0)
			return recip ? this.zero() : this.infinity();

		const t1	= this.ipow(Math.abs(n));
		const e1	= -t1.exponent % d;
		const t2	= e1
			?	this.create(((t1.exponent + e1) / d) - 1, root(t1._rep(t1.exponent + e1), d))
			:	this.create(t1.exponent / d, root(t1.mantissa, d));
		return recip ? this.one().div(t2) : t2;
	}
	npow(other: number): this {
		return Number.isInteger(other) ? this.ipow(other): this.log().mul(other).exp();
	}
	pow(other: float|number|bigint): this {
		return typeof other === "number" ? this.npow(other) : this.log().mul(other).exp();
	}
	/*
	introot(other: number): this {
		other = Math.floor(other);
		if (other < 1)
			return other < 0 ? this.one().div(this.introot(-other)) : this.infinity();
		const e1 = -this.exponent % other;
		if (e1)
			return this.create(((this.exponent + e1) / other) - 1, root(this._rep(this.exponent + e1), other));
		return this.create(this.exponent / other, root(this.mantissa, other));
	}
	root(other: float|number|bigint): this {
		return typeof other === "number" && Number.isInteger(other)
			? this.introot(other)
			: this.log().div(other).exp();
	}
*/
	toInt(mode: RoundMode = Round.trunc): bigint {
		return this.exponent < 0
			? round(this.mantissa, -this.exponent, mode)
			: this.mantissa << BigInt(this.exponent);
	}

	frac(): this {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			return this.create(this.exponent, this.mantissa % mask);
		}
		return this.zero();
	}
	floor(): this {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			if (this.mantissa % mask)
				return this.create(this.exponent, (this.mantissa < 0 ? this.mantissa - mask : this.mantissa) / mask * mask);
		}
		return this;
	}
	ceil(): this {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			if (this.mantissa % mask)
				return this.create(this.exponent, (this.mantissa < 0 ? this.mantissa : this.mantissa - mask) / mask * mask);
		}
		return this;
	}
	trunc(): this {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			return this.create(this.exponent, this.mantissa - (this.mantissa % mask));
		}
		return this;
	}
	round(): this {
		if (this.exponent < 0) {
			const mask = 10n ** BigInt(-this.exponent);
			const m = this.mantissa + (mask >> 1n);
			return this.create(this.exponent, m / mask * mask);
		}
		return this;
	}

	compare(other: float|number|bigint): number {
		other = float.from(other);
		return this.exponent === other.exponent	? compare(this.mantissa, other.mantissa)
			: this.exponent < other.exponent	? (other.exponent === Infinity ? -1 : compare(this.mantissa, other._rep(this.exponent)))
			: this.exponent === Infinity		? 1
			: compare(this._rep(other.exponent), other.mantissa);
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

	toString(base = 10, max_digits?: number): string {
		if (!this.mantissa)
			return '0';
		
		const	s = this.mantissa < 0n ? '-' : '';
		let		e = this.exponent;

		if (e > 0) {
			if (e === Infinity)
				return s + 'Infinity';
			return this._rep(0).toString(base);
		}

		let m = this.mantissa < 0n ? -this.mantissa : this.mantissa;

		if (max_digits !== undefined) {
			if (-e > max_digits) {
				m = this._rep(-max_digits);
				e = -max_digits;
			}
		}

		const p = 10n ** BigInt(-e);
		const x = m / p, y = m % p;
		if (!y)
			return s + x.toString(base);

		return s + x.toString(base) + '.' + (p + y).toString(base).slice(1);
	}
	[Symbol.for("debug.description")]() { return this.toString(10, 10); }
}

//-----------------------------------------------------------------------------
// pi
//-----------------------------------------------------------------------------

//Gauss-Legendre iterative algorithm for pi
function pi_helper(digits: number) {
	let a = float.one.addPrecision(digits);
	let b = a.div(float.two.addPrecision(digits << 1).sqrt());
	let t = float.one;
	const bits = Math.ceil((digits + 4) * log2_10); // = (digits + 4) * log2(10)

	for (let p = 0; (8 << p) < bits; ++p) {
		t = t.sub(b.sub(a).square().mul(1 << p)).capPrecision(digits + 4);

		const ab = a.mul(b);
		a = a.add(b).shift(-1);
		b = ab.sqrt();
	}

	return a.add(b).square().div(t).setPrecision(digits);
}

const pis: float[] = [];

function sin_helper(x: float, pi: float, digits: number): float {
	// Reduce to [-pi, pi]
	const twoPi = pi.shift(1);
	x = x.mod(twoPi);
	if (x.sign() < 0)
		x = x.add(twoPi);

	// Further reduce to [-pi/2, pi/2]
	if (x.gt(pi.shift(-1)))
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

//-----------------------------------------------------------------------------
// asin,acos,atan
//-----------------------------------------------------------------------------

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



// Taylor series for |x| <= 1
function atan_helper(x: float, digits: number): float {
	let result	= x;
	const x2	= x.mul(x).neg();
	const limit = new float(-digits - 1, 1n);

	for (let n = 3n, term = x; ; n += 2n) {
		term = term.mul(x2).capPrecision(digits + 1);
		const next = term.div(n);
		if (next.abs().lt(limit))
			break;
		result = result.add(next);
	}
	return result;
}


//-----------------------------------------------------------------------------
// log/exp
//-----------------------------------------------------------------------------

// with x in [1,2), y = (x-1)/(x+1)
function log_helper<T extends float>(y: T, digits: number): T {
	const limit = new float(-digits - 1, 1n);
	const y2	= y.mul(y);

	let result	= y;
	for (let n = 3n, term = y; ; n += 2n) {
		term = term.mul(y2).capPrecision(digits + 1);
		const next = term.div(n);
		if (next.abs().lt(limit))
			break;
		result = result.add(next);
	}
	return result.shift(1);
}

const ln2s: float[] = [];
function ln2(digits: number): float {
	const nextpow2 = bits.highestSet(digits);
	return (ln2s[nextpow2] ??= log_helper(float.one.addPrecision(1 << nextpow2).div(3n), 1 << nextpow2)).setPrecision(digits);
}

export default float;
