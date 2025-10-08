# @isopodlabs/big
[![npm version](https://img.shields.io/npm/v/@isopodlabs/big.svg)](https://www.npmjs.com/package/@isopodlabs/big)
[![GitHub stars](https://img.shields.io/github/stars/adrianstephens/ts-big.svg?style=social)](https://github.com/adrianstephens/ts-big)
[![License](https://img.shields.io/npm/l/@isopodlabs/big.svg)](LICENSE.txt)

This module provides high-precision arithmetic and transcendental functions for big integers and custom floating-point numbers. It is designed for numerical applications requiring more precision than JavaScript's native `number` type.

## ☕ Support My Work  
If you use this package, consider [buying me a cup of tea](https://coff.ee/adrianstephens) to support future updates!  

## Features
- Arbitrary-precision integer math (`bigint`)
- Custom `float` class for high-precision floating-point arithmetic
- IEEE-754 compatible rounding modes
- High-precision implementations of:
  - Square root, n-th root
  - Trigonometric functions: sin, cos, tan, asin, acos, atan
  - Exponential and logarithmic functions: exp, ln
  - Random number generation with specified bit width
  - Pi calculation using the Gauss-Legendre algorithm
- Utility functions for max, min, shifting, and comparison

## Usage
Import the module:
```typescript
import * as big from './big';
```

### Creating Floats
```typescript
const a = big.float.from(123.456); // From number
const b = big.float.from('3.14159'); // From string
const c = big.float.from(12345678901234567890n); // From bigint
```

### Arithmetic
```typescript
const sum = a.add(b);
const product = a.mul(b);
const quotient = a.div(b);
```

### Transcendental Functions
```typescript
const s = big.sin(a, bits); // Sine
const c = big.cos(a, bits); // Cosine
const t = big.tan(a, bits); // Tangent
const l = big.log(a, bits);  // Natural logarithm
const e = big.exp(a, bits); // Exponential
```

### Random Number Generation
```typescript
const r = big.random(bits); // Random float with specified bits of precision
```

### Pi Calculation
```typescript
const pi = big.pi(bits); // High-precision pi
```

## API Reference

### bigint functions
These are used by the float class internally, but are exposed for direct roots of bigints.
 - `sqrt(x: bigint)` - Square root
 - `root(x: bigint, b: number)` — n-th root

### Rounding Modes
These can be passed to precision modifying functions and toInt.
- `Round.trunc`: Truncate toward zero
- `Round.down`: Round toward −∞
- `Round.up`: Round toward +∞
- `Round.halfUp`: Round to nearest, ties away from zero
- `Round.halfEven`: Round to nearest, ties to even (default IEEE-754)


### `float` class methods

#### Construction
- `float.from(value: number | bigint | string | float): float` — Create a float from various types.
- `float.fromString(str: string, bits: number): float` — Parse a string to float with specified precision.

#### Arithmetic
- `add(other: float | number | bigint): float` — Addition
- `sub(other: float | number | bigint): float` — Subtraction
- `mul(other: float | number | bigint): float` — Multiplication
- `div(other: float | number | bigint): float` — Division
- `mod(other: float | number | bigint): float` — Modulus
- `divmod(other: float | number | bigint): [bigint, float]` — Whole part of division and remainder
- `square(): float` — Square
- `sqrt(): float` — Square root
- `pow(exp: number): float` — Power
- `root(base: number): float` — n-th root

#### Precision and Shifting
- `addPrecision(bits: number, mode = Round.halfEven): float` — Increase precision
- `setPrecision(bits: number, mode = Round.halfEven): float` — Set precision
- `capPrecision(bits: number, mode = Round.halfEven): float` — Cap precision
- `shift(bits: number): float` — multiply by powers of 2

#### Rounding and Integer Conversion
- `toInt(mode = Round.trunc): bigint` — Convert to integer
- `frac(): float` — Fractional part
- `floor(): float` — Round down
- `ceil(): float` — Round up
- `trunc(): float` — Truncate toward zero
- `round(): float` — Round to nearest

#### Comparison
- `compare(other: float | number | bigint): number` — Compare values
- `lt0(): boolean` — Less than zero
- `le(other: float | number | bigint): boolean` — Less than or equal
- `eq(other: float | number | bigint): boolean` — Equal
- `ne(other: float | number | bigint): boolean` — Not equal
- `ge(other: float | number | bigint): boolean` — Greater than or equal
- `lt(other: float | number | bigint): boolean` — Less than
- `gt(other: float | number | bigint): boolean` — Greater than

#### Utility
- `neg(): float` — Negate
- `abs(): float` — Absolute value
- `toString(base?: number, max_digits?: number): string` — String representation

#### Static Properties
- `float.zero` — 0
- `float.one` — 1
- `float.two` — 2
- `float.Infinity` — Infinity

### functions

- `max(...values: float[]): float` - Maximum of all values
- `min(...values: float[]): float` - Minimum of all values
- `random(bits: number): float` - Random number with given number of bits (use Math.random internally)
- `pi(bits: number): float` - Pi to the given number of binary places

These correspond to the `Math` functions, but compute results to `bits` binary places.
- `sin(x: float | number | bigint, bits: number): float`
- `cos(x: float | number | bigint, bits: number): float`
- `tan(x: float | number | bigint, bits: number): float`
- `asin(x: float | number | bigint, bits: number): float`
- `acos(x: float | number | bigint, bits: number): float`
- `atan(x: float | number | bigint, bits: number): float`
- `log(x: float | number | bigint, nbits: number): float`
- `exp(x: float | number | bigint, bits: number): float`

## On Precision

`float.from(string)` returns a float with the minimum precision required to hold the number of digits provided. Fractional values that do not have exact binary representations may not have the accuracy you expect. Use `float.fromString` to supply the precision at time of parsing.

eg.
```typescript
const a = big.float.from('0.1');
console.log(a.toString()); // outputs '0.1' - huzzah
console.log(a.addPrecision(100).toString()); // outputs '0.0625' - wtf?
const b = big.float.fromString('0.1', 100);
console.log(b.toString());
console.log(b.addPrecision(100).toString()); // still outputs '0.1'
console.log(b.addPrecision(100).toString(10, Infinity)); // outputs 0.0999999999999999999999999999999704177160542120572970601788019 (only ~30 valid digits)

```

- `add`/`sub`/`mod`/`divmod` return floats with the higher precision of `this` and `other`
- `mul` returns a float with the precision of `this` plus the precision of `other`
- `div` returns a float with the precision of `this`
- `sqrt` returns a float with half the precision of `this`
- `pow` returns a float with n times the precision of `this`
- `root` returns a float with the precision of `this` divided by `base`
- `frac`/`floor`/`ceil`/`trunc`/`round` return floats with the same precision as `this`

To compute a result to a higher precision, increase the precision of the `this` (or `other`) before calling a method.

To reduce the precision of a result, use `capPrecision` or `setPrecision` on it.
 
e.g.
```typescript
const sqrt2 = big.float.from(2).setPrecision(2000).sqrt(); // Calculate sqrt2 to 1000 binary places
const ab = a.mul(b).capPrecision(100); // cap the precision to 100 binary places
```

## License
MIT
