/* eslint-disable @typescript-eslint/no-unsafe-declaration-merging */
import { bits } from '@isopodlabs/utilities';

//-----------------------------------------------------------------------------
// bigint functions
//-----------------------------------------------------------------------------

// Progressive square root for bigints
export function sqrt(x: bigint): bigint {
	if (x <= 0n)
		return 0n;

	const	resultBits = bits.highestSet(x) >> 1;

	let k	= 16;
	let y	= BigInt(Math.floor(Math.sqrt(Number(x >> BigInt((resultBits - k) * 2))))); // Initial guess

	while (k * 2 < resultBits) {
		const xk = x >> BigInt(resultBits * 2 - k * 3);
		y = ((y << BigInt(k)) + xk / y) >> 1n;
		k <<= 1;
	}

	// Final refinement at full precision
	y = y << BigInt(resultBits - k);
	let y1 = x / y;
	y = (y + y1) >> 1n;

	while (y > y1) {
		y1 = x / y;
		y = (y + y1) >> 1n;
	};

	return y;
}

// Progressive n-th root for bigints
export function root(x: bigint, b: number): bigint {
	if (b < 1)
		return 0n;
	if ((b & 1) === 0 && x < 0n)
		return 0n;

	const resultBits = Math.floor(bits.highestSet(x) / b);

	let k = 16;
	let y = BigInt(Math.floor(Math.pow(Number(x >> BigInt((resultBits - k) * b)), 1 / b)));

	const b1 = BigInt(b - 1);
	const bb = BigInt(b);

	while (k * 2 < resultBits) {
		const xk = x >> BigInt((resultBits - k) * b - k);
		y = (((b1 * y) << BigInt(k)) + xk / (y ** b1)) / bb;
		k <<= 1;
	}

	// Final refinement at full precision
	y = y << BigInt(resultBits - k);
	y = (b1 * y + x / (y ** b1)) / bb;

	while ((y ** bb) > x)
		--y;
	return x < 0n ? -y : y;
}

export function compare(a: bigint, b: bigint) {
	return a === b ? 0 : a < b ? -1 : 1;
}

export function randomBits(bits: number) {
	let m = 0n;
	let i = bits;
	while (i > 48) {
		m = (m << 48n) | BigInt(Math.floor(Math.random() * 0x1000000000000));
		i -= 48;
	}
	return (m << BigInt(i)) | BigInt(Math.floor(Math.random() * 0x1000000000000)) & ((1n << BigInt(i)) - 1n);
}

export const Round = {
	trunc:		0,	//	Rounds towards zero.I.e. truncate, no rounding
	down:		1,	//	Rounds towards −∞	(floor)
	up:			2,	//	Rounds towards +∞	(ceil)
	halfUp:		3,	//	Rounds towards nearest neighbour.If equidistant, rounds away from zero
	halfEven:	4,	//	Rounds towards nearest neighbour.If equidistant, rounds towards even neighbour
} as const;

export type RoundMode = (typeof Round)[keyof typeof Round];

function round(m: bigint, n: number, mode: RoundMode): bigint {
	const	neg = m < 0n;
	let		a	= neg ? -m : m;
	const	r	= 1n << BigInt(n - 1);

	switch (mode) {
		case Round.trunc:
			break;

		case Round.down:
			if (neg && (a & ((r << 1n) - 1n)))
				a += r << 1n;
			break;

		case Round.up:
			if (!neg && (a & ((r << 1n) - 1n)))
				a += r << 1n;
			break;

		case Round.halfUp:
			a += r;
			break;

		case Round.halfEven:
			if ((a & r) && (a & (r << 2n) - 1n))
				a += r;
			break;

	}

	a >>= BigInt(n);
	return neg ? -a : a;
}

function shift(m: bigint, n: number) {
	return n > 0 
		? m << BigInt(n)
		: m >> BigInt(-n);
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
		let minv = this.Infinity;
		for (const i of values) {
			if (i.lt(minv))
				minv = i;
		}
		return minv;
	}

	static random(bits: number) {
		return new float(-bits, randomBits(bits));
	}
	static pi(nbits: number): float {
		const nextpow2 = bits.highestSet(nbits);
		return (pis[nextpow2] ??= pi_helper(1 << nextpow2)).setPrecision(nbits);
	}

	static sin(x: float|number|bigint, bits?: number): float {
		x = this.from(x);
		bits ??= -x.exponent;
		const π = this.pi(bits + 8);
		return sin_helper(x, π, bits).setPrecision(bits);
	}

	static cos(x: float|number|bigint, bits?: number): float {
		x = this.from(x);
		bits ??= -x.exponent;
		const π = this.pi(bits + 8);
		return sin_helper(x.add(π.shift(-1)), π, bits).setPrecision(bits);
	}

	static tan(x: float|number|bigint, bits?: number): float {
		x = this.from(x);
		bits ??= -x.exponent;
		const π = this.pi(bits + 8);
		return sin_helper(x, π, bits).div(sin_helper(x.add(π.shift(-1)), π, bits)).setPrecision(bits);
	}

	static asin(x: float|number|bigint, bits?: number): float {
		x = this.from(x);
		bits ??= -x.exponent;

		if (x.lt(0.707))
			return asin_helper(x, bits).setPrecision(bits);

		// arcsin(x) = pi/2 - 2*arcsin(sqrt((1-x)/2))
		const π		= this.pi(bits + 8);
		const asin1	= asin_helper(this.one.sub(x).shift(-1).sqrt(), bits);
		const res	= π.shift(-1).sub(asin1.shift(1)).setPrecision(bits);
		return x.sign() < 0 ? res.neg() : res;
	}

	// Use identity arccos(x) = pi/2 - arcsin(x)
	static acos(x: float|number|bigint, bits?: number): float {
		x = this.from(x);
		bits ??= -x.exponent;
		const π		= this.pi(bits + 8);
		if (x.lt(0.707))
			return π.shift(-1).sub(asin_helper(x, bits)).setPrecision(bits);

		const asin	= asin_helper(this.one.sub(x).shift(-1).sqrt(), bits);
		return x.sign() < 0
			? π.sub(asin.shift(1)).setPrecision(bits)
			: asin.shift(1).setPrecision(bits);
	}


	static atan(x: float|number|bigint|string, bits?: number): float {
		x = this.from(x);
		bits ??= -x.exponent;

		if (x.abs().le(this.one))
			return atan_helper(x, bits).setPrecision(bits);

		// Argument reduction for |x| > 1
		const π_2	= this.pi(bits + 8).shift(-1);
		const a		= atan_helper(x.recip(), bits);
		return (x.sign() < 0 ? π_2.neg() : π_2).sub(a).setPrecision(bits);
	}

	static atan2(y: float|number|bigint|string, x: float|number|bigint|string, bits?: number): float {
		y = this.from(y);
		x = this.from(x);
		bits ??= -x.exponent;
		const π_2 = this.pi(bits + 8).shift(-1);

		if (x.mantissa === 0n) {
			switch (y.sign()) {
				case 1:		return π_2;
				case -1:	return π_2.neg();
				default:	return this.zero;
			}
		}

		if (x.abs().gt(y.abs())) {
			const a = atan_helper(y.div(x), bits);
			return x.sign() > 0 ? a : a.add(y.sign() >= 0 ? π_2 : π_2.neg());
		} else {
			const a = atan_helper(x.div(y), bits);
			return y.sign() > 0 ? π_2.sub(a) : π_2.neg().sub(a);
		}
	}

	static log(x: float|number|bigint, bits?: number): float {
		x = this.from(x, bits);
		return x.log();
	}

	static exp(x: float|number|bigint, bits: number): float {
		x = this.from(x, bits);
		return x.exp();
	}

	static from<C extends new (exponent: number, mantissa: bigint) => any>(this: C, v: number|bigint|string|float, nbits?: number): InstanceType<C> {
		function withPrecision(f: float): InstanceType<C> {
			return (nbits === undefined ? f : f.setPrecision(nbits)) as InstanceType<C>;
		}
		switch (typeof v) {
			case 'object':
				return withPrecision(v instanceof this ? v : new this(v.exponent, v.mantissa));

			case 'bigint':
				return withPrecision(new this(0, v));

			case 'number': {
				if (v === Infinity)
					return new this(Infinity, 1n);
				if (v === -Infinity)
					return new this(Infinity, -1n);
				if (Number.isInteger(v))
					return withPrecision(new this(0, BigInt(v)));
				
				const a = Math.abs(v);
				let e = Math.floor(Math.log2(a)) - 52;
				let m = BigInt(Math.floor(a * Math.pow(2, -e)));
				const zeros = bits.lowestSet(m);
				e += zeros;
				m >>= BigInt(zeros);
				return withPrecision(new this(e, v < 0 ? -m : m));
			}

			case 'string': {
				const m = v.match(/^(-?)(\d+)(?:\.(\d*))?(?:e([-+]?\d+))?$/);
				if (m) {
					let 	x	= BigInt(m[2] + (m[3] ? m[3] : ''));
					const	e10	= (m[4] ? parseInt(m[4]) : 0) - (m[3] ?? '').length;
					const	e	= Math.floor(e10 * 3.321928094887362) - (nbits ?? 0); // log2(10)
					if (e10 < 0) {
						x = (x << BigInt(-e)) / (10n ** BigInt(-e10));
					} else {
						x = shift(x * 10n ** BigInt(e10), -e);
					}
					return new this(e, m[1] === '-' ? -x : x);
				}

				return withPrecision(new this(0, BigInt(v)));
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
		const k = bits.highestSet(this.mantissa) + this.exponent - 1;
		const m = this.shift(-k);

		const nbits = Math.max(-this.exponent, 32);
		const result = log_helper(m.sub(float.one).addPrecision(nbits).div(m.add(float.one)), nbits);

		// ln(x) = ln(m) + k*ln(2)
		return result.add(ln2(nbits).mul(k)).setPrecision(nbits);
	}

	exp(): this {
		if (this.exponent === Infinity && this.mantissa < 0n)
			return this.zero();

		const x		= this.create(this.exponent, this.mantissa);
		const nbits	= Math.max(-x.exponent, 32);
		const shift = x.divmod(ln2(nbits));

		// Taylor series for exp(xred)
		const limit	= this.create(-nbits - 4, 1n);
		let result	= this.one();
		for (let n = 1n, term = this.one(); ; n++) {
			term = term.mul(x).div(n).capPrecision(nbits + 4);
			if (term.abs().lt(limit))
				break;
			result = result.add(term);
		}
		
		// Undo argument reduction
		return result.setPrecision(nbits).shift(Number(shift));
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

	private _rep(exponent: number): bigint {
		return this.mantissa << BigInt(this.exponent - exponent);
	}

	addPrecision(p: number, mode: RoundMode = Round.halfEven): this {
		if (p === Infinity)
			return this.create(Infinity, this.mantissa < 0 ? -1n : 1n);
		return this.create(this.exponent - p, p < 0 ? round(this.mantissa, -p, mode) : this.mantissa << BigInt(p));
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
		return this.create(this.exponent + p, this.mantissa);
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
		const e = this.exponent - this.divPrecision;
		return other.exponent === Infinity		? this.zero()
			: other.mantissa === 0n				? this.infinity()
			: this.create(e, this._rep(e + other.exponent) / other.mantissa);
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
		if (other.mantissa === 0n || this.exponent == Infinity)
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
		return this.exponent === Infinity		? this.create(Infinity, 1n)
			: this.create(this.exponent >> 1, sqrt(this.mantissa << (this.exponent & 1 ? 1n : 0n)));
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
		let n1 = Math.floor(other), d1 = 1;
		let x = other - n1;

		for (let n2 = 1, d2 = 0; d1 < 100 && Math.abs(x) > 1e-15; ) {
			x = 1 / x;
			const f = Math.floor(x);
			[n2, n1, d2, d1] = [n1, f * n1 + n2, d1, f * d1 + d2];
			x -= f;
		}
		return d1 === 1 ? this.ipow(n1) : d1 < 100 ? this.rpow(n1, d1) : this.log().mul(other).exp();
	}
	pow(other: float|number|bigint): this {
		return typeof other === "number" ? this.npow(other) : this.log().mul(other).exp();
	}

	toInt(mode: RoundMode = Round.trunc): bigint {
		return this.exponent < 0
			? round(this.mantissa, -this.exponent, mode)
			: this.mantissa << BigInt(this.exponent);
	}

	frac(): this {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			return this.create(this.exponent, this.mantissa & (mask - 1n));
		}
		return this.zero();
	}
	floor(): this {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			if (this.mantissa & (mask - 1n))
				return this.create(this.exponent, (this.mantissa & -mask) - (this.mantissa < 0 ? mask : 0n));
		}
		return this;
	}
	ceil(): this {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			if (this.mantissa & (mask - 1n))
				return this.create(this.exponent, (this.mantissa & -mask) - (this.mantissa < 0 ? 0n : mask));
		}
		return this;
	}
	trunc(): this {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			return this.create(this.exponent, this.mantissa & -mask);
		}
		return this;
	}
	round(): this {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			const m = this.mantissa + (mask >> 1n);
			return this.create(this.exponent, m & -mask);
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
		
		let		e = this.exponent;
		let		m = this.mantissa;

		if (e > 0) {
			if (e === Infinity)
				return (m < 0n ? '-' : '') + 'Infinity';
			return this._rep(0).toString(base);
		}

		if (max_digits !== undefined) {
			const max_bits = Math.ceil(max_digits * Math.log2(base));
			if (-e > max_bits) {
				m >>= BigInt(-e - max_bits);
				e = -max_bits;
			}
		}

		let s = '';
		if (m < 0n) {
			s = '-';
			m = -m;
		}

		const bits = BigInt(-e);
		const mask = (1n << bits) - 1n;
		if (mask)
			++m;// round last digit

		s += (m >> bits).toString(base);

		m &= mask;
		if (m) {
			const digits = '0123456789abcdefghijklmnopqrstuvwxyz';
			s += '.';
		
			const bbase = BigInt(base);
			for (let e = 2n; m >= e; e *= bbase) {
				m *= bbase;
				s += digits[Number(m >> bits)];
				m &= mask;
			}
		}
		return s;
	}

	mag(): this {
		return this.abs();
	}
	valueOf(): number {
		let m = this.mantissa;
		let e = this.exponent;
		if (e < -53) {
			m = (m >> BigInt(-e - 53)) + 1n;
			e = -53;
		}
		return Number(m) * Math.pow(2, e);
	}
	[Symbol.for("debug.description")]() { return this.toString(10, 10); }

}
export interface float {
	divPrecision: number;
}
float.prototype.divPrecision = 53;

// as above but mul & intpow won't explode mantissa

export class float2 extends float {
	protected create(exponent: number, mantissa: bigint): this {
		return new (this.constructor as new (exponent: number, mantissa: bigint) => this)(exponent, mantissa);
	}

	mul(other: float|number|bigint): this {
		other = float.from(other);
		const precision = Math.max(-this.exponent, -other.exponent);
		return this.create(this.exponent + other.exponent, this.mantissa * other.mantissa).capPrecision(precision);
	}
	ipow(other: number): this {
		return super.ipow(other).capPrecision(-this.exponent);
	}
}

//-----------------------------------------------------------------------------
// pi
//-----------------------------------------------------------------------------

//Gauss-Legendre iterative algorithm for pi
export function pi_helper(bits: number) {
	let a = float.one.addPrecision(bits);
	let b = a.div(float.two.addPrecision(bits << 1).sqrt());
	let t = float.one;

	for (let p = 0; (8 << p) < bits; ++p) {
		t = t.sub(b.sub(a).square().shift(p)).capPrecision(bits + 4);

		const ab = a.mul(b);
		a = a.add(b).shift(-1);
		b = ab.sqrt();
	}

	return a.add(b).square().div(t).setPrecision(bits);
}

const pis: float[] = [];


//-----------------------------------------------------------------------------
// sin/cos/tan
//-----------------------------------------------------------------------------

function sin_helper(x: float, pi: float, bits: number): float {
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
	const limit = new float(-bits - 4, 1n);

	for (let n = 2n, term = x; ; n += 2n) {
		term = term.mul(x2).div(n * (n + 1n)).capPrecision(bits + 4);
		if (term.abs().lt(limit))
			break;
		result	= result.add(term);
	}
	return result;
}

//-----------------------------------------------------------------------------
// asin,acos,atan
//-----------------------------------------------------------------------------

function asin_helper(x: float, bits: number): float {
	const limit = new float(-bits - 4, 1n);
	const x2	= x.mul(x);
	let result	= x;

	for (let n = 2n, num = x, den = 1n; ; n += 2n) {
		num = num.mul(x2).capPrecision(bits + 4).mul((n - 1n) * (n - 1n));
		den *= n * (n + 1n);
		const term = num.div(den);
		if (term.abs().lt(limit))
			break;
		result = result.add(term);
	}
	return result;
}


// Taylor series for |x| <= 1
function atan_helper(x: float, bits: number): float {
	let result	= x;
	const x2	= x.mul(x).neg();
	const limit = new float(-bits - 4, 1n);

	for (let n = 3n, term = x; ; n += 2n) {
		term = term.mul(x2).capPrecision(bits + 4);
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
function log_helper<T extends float>(y: T, bits: number): T {
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
	return result.shift(1);
}

const ln2s: float[] = [];
function ln2(nbits: number): float {
	const nextpow2 = bits.highestSet(nbits);
	return (ln2s[nextpow2] ??= log_helper(float.one.addPrecision(1 << nextpow2).div(3n), 1 << nextpow2)).setPrecision(nbits);
}

export default float;
