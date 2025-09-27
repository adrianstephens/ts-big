import { bits } from '@isopodlabs/utilities';

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
	// Initial guess: n-th root of truncated x
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

function shift(m: bigint, n: number) {
	return n > 0 
		? m << BigInt(n)
		: m >> BigInt(-n);
}
function compare(a: bigint, b: bigint) {
	return a === b ? 0 : a < b ? -1 : 1;
}

export const Round = {
	trunc:		0,	//	Rounds towards zero.I.e. truncate, no rounding
	down:		1,	//	Rounds towards −∞	(floor)
	up:			2,	//	Rounds towards +∞	(ceil)
	halfUp:		3,	//	Rounds towards nearest neighbour.If equidistant, rounds away from zero
	halfEven:	4,	//	Rounds towards nearest neighbour.If equidistant, rounds towards even neighbour
} as const;

type RoundMode = 0 | 1 | 2 | 3 | 4;

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

export class float {
	static readonly zero	= new float(0, 0n);
	static readonly one		= new float(0, 1n);
	static readonly two		= new float(0, 2n);
	static readonly Infinity= new float(Infinity, 1n);

	static fromString(v: string, bits: number): float {
		const m = v.match(/^(-?)(\d+)(?:\.(\d*))?(?:e([-+]?\d+))?$/);
		if (m) {
			let 	x	= BigInt(m[2] + (m[3] ? m[3] : ''));
			const	e10	= (m[4] ? parseInt(m[4]) : 0) - (m[3] ?? '').length;
			const	e	= Math.floor(e10 * 3.321928094887362) - bits; // log2(10)
			if (e10 < 0) {
				x = (x << BigInt(-e)) / (10n ** BigInt(-e10));
			} else {
				x = shift(x * 10n ** BigInt(e10), -e);
			}
			return new float(e, m[1] === '-' ? -x : x);
		}

		return new float(-bits, BigInt(v) << BigInt(bits));
	}

	static from(v: number|bigint|string|float): float {
		switch (typeof v) {
			case 'object':
				return v;

			case 'bigint':
				return new float(0, v);

			case 'number': {
				if (v === Infinity)
					return this.Infinity;
				if (v === -Infinity)
					return this.Infinity.neg();
				if (v === Math.floor(v))
					return new float(0, BigInt(v));
				
				const a = Math.abs(v);
				let e = Math.floor(Math.log2(a)) - 52;
				let m = BigInt(Math.floor(a * Math.pow(2, -e)));
				const zeros = bits.lowestSet(m);
				e += zeros;
				m >>= BigInt(zeros);
				return new float(e, v < 0 ? -m : m);
			}

			case 'string':
				return this.fromString(v, 0);
		}
	}

	constructor(public exponent: number, public mantissa: bigint) {
	}

	addPrecision(p: number, mode: RoundMode = Round.halfEven): float {
		if (p === Infinity)
			return this.mantissa < 0 ? float.Infinity.neg() : float.Infinity;
		return new float(this.exponent - p, p < 0 ? round(this.mantissa, -p, mode) : this.mantissa << BigInt(p));
	}
	setPrecision(p: number, mode: RoundMode = Round.halfEven): float {
		return this.addPrecision(this.exponent + p, mode);
	}
	capPrecision(p: number, mode: RoundMode = Round.halfEven): float {
		if (this.exponent + p < 0)
			return this.addPrecision(this.exponent + p, mode);
		return this;
	}
	shift(p: number): float {
		return new float(this.exponent + p, this.mantissa);
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
		return new float(a.exponent, a.mantissa + (b.mantissa << BigInt(b.exponent - a.exponent)));
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
			: new float(this.exponent, (this.mantissa << BigInt(-other.exponent)) / other.mantissa);
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
			: new float(this.exponent >> 1, sqrt(this.mantissa << (this.exponent & 1 ? 1n : 0n)));
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
			return new float(((this.exponent + e1) / base) - 1, root(this.mantissa << BigInt(base - e1), base));
		return new float(this.exponent / base, root(this.mantissa, base));
	}

	toInt(mode: RoundMode = Round.trunc): bigint {
		return this.exponent < 0
			? round(this.mantissa, -this.exponent, mode)
			: this.mantissa << BigInt(this.exponent);
	}

	frac(): float {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			return new float(this.exponent, this.mantissa & (mask - 1n));
		}
		return float.zero;
	}
	floor(): float {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			if (this.mantissa & (mask - 1n))
				return new float(this.exponent, (this.mantissa & -mask) - (this.mantissa < 0 ? mask : 0n));
		}
		return this;
	}
	ceil(): float {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			if (this.mantissa & (mask - 1n))
				return new float(this.exponent, (this.mantissa & -mask) - (this.mantissa < 0 ? 0n : mask));
		}
		return this;
	}
	trunc(): float {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			return new float(this.exponent, this.mantissa & -mask);
		}
		return this;
	}
	round(): float {
		if (this.exponent < 0) {
			const mask = 1n << BigInt(-this.exponent);
			const m = this.mantissa + (mask >> 1n);
			return new float(this.exponent, m & -mask);
		}
		return this;
	}

	compare(other: float|number|bigint): number {
		other = float.from(other);
		return this.exponent === other.exponent
			? compare(this.mantissa, other.mantissa)
			: this.exponent < other.exponent ? (other.exponent === Infinity ? -1 : compare(this.mantissa, other.mantissa << BigInt(other.exponent - this.exponent)))
			: this.exponent === Infinity ? 1 : compare(this.mantissa << BigInt(this.exponent - other.exponent), other.mantissa);
	}

	lt0(): boolean {
		return this.mantissa < 0n;
	}
	le(other: float|number|bigint): boolean {
		other = float.from(other);
		return this.exponent === other.exponent	? this.mantissa <= other.mantissa
			: this.exponent < other.exponent	? other.exponent === Infinity || this.mantissa <= (other.mantissa << BigInt(other.exponent - this.exponent))
			: this.exponent !== Infinity && this.mantissa << BigInt(this.exponent - other.exponent) <= other.mantissa;
	}
	eq(other: float|number|bigint): boolean {
		other = float.from(other);
		return this.exponent === other.exponent	? this.mantissa === other.mantissa
			: this.exponent < other.exponent	? other.exponent !== Infinity && this.mantissa === (other.mantissa << BigInt(other.exponent - this.exponent))
			: this.exponent !== Infinity && this.mantissa << BigInt(this.exponent - other.exponent) === other.mantissa;
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

	toString(base = 10, max_digits = 10): string {
		if (!this.mantissa)
			return '0';
		
		const	e = this.exponent;
		let		m = this.mantissa;

		if (e > 0) {
			if (e === Infinity)
				return (m < 0n ? '-' : '') + 'Infinity';
			return (m << BigInt(e)).toString(base);
		}

		const max_bits = Math.ceil(max_digits * Math.log2(base));
		if (-e > max_bits)
			m >>= BigInt(-e - max_bits);

		const bits = BigInt(Math.min(-e, max_bits));
//		const bits = BigInt(-e);
		const mask = (1n << bits) - 1n;

		if (mask)
			++m;// round last digit

		const digits = '0123456789abcdefghijklmnopqrstuvwxyz';
		let s = '';
		if (m < 0n) {
			s = '-';
			m = -m;
		}
		s += (m >> bits).toString(base);

		m &= mask;
		if (m) {
			s += '.';
		
			let e = 2n;
			//m += 1n;// round last digit
			const bbase = BigInt(base);
			while (m) {//} && (max_digits-- > 0)) {
				m *= bbase;
				e *= bbase;
				s += digits[Number(m >> bits)];
				m &= mask;
				if (m < e)
					break;
			}
		}
		return s;
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

export function random(bits: number) {
	let m = 0n;
	let i = bits;
	while (i > 48) {
		m = (m << 48n) | BigInt(Math.floor(Math.random() * 0x1000000000000));
		i -= 48;
	}
	m = (m << BigInt(i)) | BigInt(Math.floor(Math.random() * 0x1000000000000)) & ((1n << BigInt(i)) - 1n);
	return new float(-bits, m);
}

//Gauss-Legendre iterative algorithm for pi
export function pi(bits: number) {
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

function sin_helper(x: float, pi: float, bits: number): float {
	// Reduce to [-pi, pi]
	const twoPi = pi.shift(1);
	x = x.mod(twoPi);
	if (x.lt0())
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

export function sin(x: float|number|bigint, bits: number): float {
	const π = pi(bits + 8);
	return sin_helper(float.from(x), π, bits).setPrecision(bits);
}

export function cos(x: float|number|bigint, bits: number): float {
	const π = pi(bits + 8);
	return sin_helper(float.from(x).add(π.shift(-1)), π, bits).setPrecision(bits);
}

export function tan(x: float|number|bigint, bits: number): float {
	x = float.from(x);
	const π = pi(bits + 8);
	return sin_helper(x, π, bits).div(sin_helper(x.add(π.shift(-1)), π, bits)).setPrecision(bits);
}

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

// High-precision arcsin using Taylor series and argument reduction
export function asin(x: float|number|bigint, bits: number): float {
	x = float.from(x);

	if (x.lt(0.707))
		return asin_helper(x, bits).setPrecision(bits);

	// arcsin(x) = pi/2 - 2*arcsin(sqrt((1-x)/2))
	const π		= pi(bits + 8);
	const asin1	= asin_helper(float.one.sub(x).shift(-1).sqrt(), bits);
	const res	= π.shift(-1).sub(asin1.shift(1)).setPrecision(bits);
	return x.lt0() ? res.neg() : res;
}

// Use identity arccos(x) = pi/2 - arcsin(x)
export function acos(x: float|number|bigint, bits: number): float {
	const π		= pi(bits + 8);
	x = float.from(x);

	if (x.lt(0.707))
		return π.shift(-1).sub(asin_helper(x, bits)).setPrecision(bits);

	const asin	= asin_helper(float.one.sub(x).shift(-1).sqrt(), bits);
	return x.lt0()
		? π.sub(asin.shift(1)).setPrecision(bits)
		: asin.shift(1).setPrecision(bits);
}


// High-precision arctan using Taylor series and argument reduction
export function atan(x: float|number|bigint|string, bits: number): float {
	x = float.from(x);
	const limit = new float(-bits - 4, 1n);

	function helper(x: float): float {
		// Taylor series for |x| <= 1
		let result	= x;
		const x2 = x.mul(x).neg();

		for (let n = 3n, term = x; ; n += 2n) {
			term = term.mul(x2).capPrecision(bits + 4);
			const next = term.div(n);
			if (next.abs().lt(limit))
				break;
			result = result.add(next);
		}
		return result;
	}

	if (x.abs().le(float.one))
		return helper(x).setPrecision(bits);

	// Argument reduction for |x| > 1
	const π		= pi(bits + 8);
	const atan1	= helper(float.one.addPrecision(bits + 4).div(x));
	return (x.lt0() ? π.neg() : π).shift(-1).sub(atan1).setPrecision(bits);
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
	return result.shift(1);
}

const ln2s: float[] = [];
function ln2(nbits: number): float {
	const nextpow2 = bits.highestSet(nbits);
	return (ln2s[nextpow2] ??= log_helper(float.one.addPrecision(1 << nextpow2).div(3n), 1 << nextpow2)).setPrecision(nbits);
}

// High-precision natural logarithm (ln) and exponential (exp)
export function log(x: float|number|bigint, nbits: number): float {
	x = float.from(x);
	if (x.le(0n)) {
		if (x.lt0())
			throw new Error('ln(x): x must be positive');
		return float.Infinity.neg();
	}

	// Argument reduction: x = m * 2^k, ln(x) = ln(m) + k*ln(2)
	const k = bits.highestSet(x.mantissa) + x.exponent - 1;
	const m = x.shift(-k);

	const result = log_helper(m.sub(float.one).addPrecision(nbits).div(m.add(float.one)), nbits);

	// ln(x) = ln(m) + k*ln(2)
	return result.add(ln2(nbits).mul(k)).setPrecision(nbits);
}

export function exp(x: float|number|bigint, bits: number): float {
	const [shift, xred] = float.from(x).divmod(ln2(bits));

	// Taylor series for exp(xred)
	const limit	= new float(-bits - 4, 1n);
	let result	= float.one;
	for (let n = 1n, term = float.one; ; n++) {
		term = term.mul(xred).div(n).capPrecision(bits + 4);
		if (term.abs().lt(limit))
			break;
		result = result.add(term);
	}
	
	// Undo argument reduction
	return result.setPrecision(bits).shift(Number(shift));
}
