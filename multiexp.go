package multiexp

import (
	"context"
	"fmt"
	"math/big"
)

const defaultWordChunkSize = 2

var (
	big1  = big.NewInt(1)
	masks = [_W]Word{}
)

func init() {
	for i := 0; i < _W; i++ {
		masks[i] = 1 << i
	}
}

// DoubleExp sets z1 = x**y1 mod |m|, z2 = x**y2 mod |m| ... (i.e. the sign of m is ignored), and returns z1, z2.
// If m == nil or m == 0, z = x**y unless y <= 0 then z = 1. If m != 0, y < 0,
// and x and m are not relatively prime, z is unchanged and nil is returned.
//
// DoubleExp is not a cryptographically constant-time operation.
func DoubleExp(x *big.Int, y2 [2]*big.Int, m *big.Int) [2]*big.Int {
	// make sure x > 1, m is not nil, and m > 0, otherwise, use default Exp function
	if x.Cmp(big1) <= 0 || m == nil || m.Sign() <= 0 {
		return defaultExp2(x, m, [2]*big.Int{y2[0], y2[1]})
	}
	// make sure y1 and y2 are positive
	if y2[0].Sign() <= 0 || y2[1].Sign() <= 0 {
		return defaultExp2(x, m, y2)
	}
	// make sure m is odd
	if m.Bit(0) != 1 {
		return defaultExp2(x, m, y2)
	}
	xWords, y1Words, y2Words, mWords := newNat(x), newNat(y2[0]), newNat(y2[1]), newNat(m)
	return doubleExpNNMontgomery(xWords, y1Words, y2Words, mWords)
}

// defaultExp2 uses the default Exp function of big int to handle the edge cases that cannot be handled by DoubleExp in
// this library or cannot benefit from this library in terms of performance
func defaultExp2(x, m *big.Int, y2 [2]*big.Int) [2]*big.Int {
	fmt.Println("something wrong here, get into defaultExp2")
	var ret [2]*big.Int
	for i := range y2 {
		ret[i] = new(big.Int).Exp(x, y2[i], m)
	}
	return ret
}

// defaultExp4 uses the default Exp function of big int to handle the edge cases that cannot be handled by FourfoldExp in
// this library or cannot benefit from this library in terms of performance
func defaultExp4(x, m *big.Int, y4 [4]*big.Int) [4]*big.Int {
	var ret [4]*big.Int
	for i := range y4 {
		ret[i] = new(big.Int).Exp(x, y4[i], m)
	}
	return ret
}

// doubleExpNNMontgomery calculates x**y1 mod m and x**y2 mod m
// Uses Montgomery representation.
func doubleExpNNMontgomery(x, y1, y2, m nat) [2]*big.Int {
	power0, power1, k0, numWords := montgomerySetup(x, m)
	y1Extra, y2Extra, commonBits := gcw(y1, y2)
	mmValues := multiMontgomery(m, power0, power1, k0, numWords, []nat{y1Extra, y2Extra, commonBits})
	// calculate z1 and z2, 1st, 2nd and 3rd elements of mmValues correspond to y1Extra, y2Extra and commonBits
	temp := nat(nil).make(numWords)
	temp = temp.montgomery(mmValues[0], mmValues[2], m, k0, numWords)
	mmValues[0], temp = temp, mmValues[0]
	temp = temp.montgomery(mmValues[1], mmValues[2], m, k0, numWords)
	mmValues[1], temp = temp, mmValues[1]
	mmValues = mmValues[:2] //mm3 is useless now
	// convert to regular number
	// one = 1, with equal length to that of m
	one := make(nat, numWords)
	one[0] = 1
	for i := range mmValues {
		temp = temp.montgomery(mmValues[i], one, m, k0, numWords)
		mmValues[i], temp = temp, mmValues[i]
	}

	var ret [2]*big.Int
	for i := range mmValues {
		// One last reduction, just in case.
		// See golang.org/issue/13907.
		if mmValues[i].cmp(m) >= 0 {
			// Common case is m has high bit set; in that case,
			// since zz is the same length as m, there can be just
			// one multiple of m to remove. Just subtract.
			// We think that the subtraction should be sufficient in general,
			// so do that unconditionally, but double-check,
			// in case our beliefs are wrong.
			// The div is not expected to be reached.
			mmValues[i] = mmValues[i].sub(mmValues[i], m)
			if mmValues[i].cmp(m) >= 0 {
				_, mmValues[i] = nat(nil).div(nil, mmValues[i], m)
			}
		}
		// final normalization
		mmValues[i].norm()
		ret[i] = new(big.Int).SetBits(mmValues[i].intBits())
	}

	return ret
}

func montgomerySetup(x, m nat) (power0, power1 nat, k0 Word, numWords int) {
	numWords = len(m)
	xx := x

	// We want the lengths of x and m to be equal.
	// It is OK if x >= m as long as len(x) == len(m).
	if len(x) > numWords {
		_, x = nat(nil).div(nil, x, m)
		// Note: now len(x) <= numWords, not guaranteed ==.
	}
	if len(x) < numWords {
		rr := make(nat, numWords)
		copy(rr, x)
		x = rr
	}

	// Ideally the pre-computations would be performed outside, and reused
	// k0 = -m**-1 mod 2**_W. Algorithm from: Dumas, J.G. "On Newton–Raphson
	// Iteration for Multiplicative Inverses Modulo Prime Powers".
	k0 = 2 - m[0]
	t := m[0] - 1
	for i := 1; i < _W; i <<= 1 {
		t *= t
		k0 *= t + 1
	}
	k0 = -k0

	// RR = 2**(2*_W*len(m)) mod m
	RR := nat(nil).setWord(1)
	zz1 := nat(nil).shl(RR, uint(2*numWords*_W))
	_, RR = nat(nil).div(RR, zz1, m)
	if len(RR) < numWords {
		zz1 = zz1.make(numWords)
		copy(zz1, RR)
		RR = zz1
	}

	// one = 1, with equal length to that of m
	one := make(nat, numWords)
	one[0] = 1

	// power0 = x**0
	power0 = power0.montgomery(one, RR, m, k0, numWords)
	// power1 = x**1
	power1 = power1.montgomery(xx, RR, m, k0, numWords)
	return
}

// multiMontgomery calculates the modular montgomery exponent with result not normalized
func multiMontgomery(m, power0, power1 nat, k0 Word, numWords int, yList []nat) []nat {
	// initialize each value to be 1 (Montgomery 1)
	zList := make([]nat, len(yList))
	for i := range zList {
		zList[i] = zList[i].make(numWords)
		copy(zList[i], power0)
	}

	squaredPower := nat(nil).make(numWords)
	copy(squaredPower, power1)
	//	fmt.Println("squaredPower = ", squaredPower.String())

	maxWordLen := 1
	for i := range yList {
		if len(yList[i]) > maxWordLen {
			maxWordLen = len(yList[i])
		}
	}

	temp := nat(nil).make(numWords)
	for i := 0; i < maxWordLen; i++ {
		for j := 0; j < _W; j++ {
			for k := range yList {
				if len(yList[k]) <= i {
					continue
				}
				if (yList[k][i] & masks[j]) != masks[j] {
					continue
				}
				temp = temp.montgomery(zList[k], squaredPower, m, k0, numWords)
				zList[k], temp = temp, zList[k]
			}
			// montgomery must have the returned value not same as the input values
			// we have to use this temp as the middle variable
			temp = temp.montgomery(squaredPower, squaredPower, m, k0, numWords)
			squaredPower, temp = temp, squaredPower
		}
	}

	return zList
}

// multiMontgomeryPrecomputed calculates the modular montgomery exponent with result not normalized
func multiMontgomeryPrecomputed(m, power0 nat, k0 Word,
	numWords int, yList []nat, preTable *PreTable) []nat {
	// initialize each value to be 1 (Montgomery 1)
	z := make([]nat, len(yList))
	for i := range z {
		z[i] = z[i].make(numWords)
		copy(z[i], power0)
	}

	var temp nat
	temp = temp.make(numWords)
	//	fmt.Println("squaredPower = ", squaredPower.String())

	maxLen := 1
	for i := range yList {
		if len(yList[i]) > maxLen {
			maxLen = len(yList[i])
		}
	}

	for i := 0; i < maxLen; i++ {
		for j := 0; j < _W; j++ {
			for k := range yList {
				if len(yList[k]) <= i {
					continue
				}
				if (yList[k][i] & masks[j]) != masks[j] {
					continue
				}
				temp = temp.montgomery(z[k], preTable.table[i][j], m, k0, numWords)
				z[k], temp = temp, z[k]
			}
		}
	}
	return z
}

// FourfoldExp sets z1 = x**y1 mod |m|, z2 = x**y2 mod |m| ... (i.e. the sign of m is ignored), and returns z1, z2...
// In construction, many panic conditions. Use at your own risk!
//
// FourfoldExp is not a cryptographically constant-time operation.
func FourfoldExp(x, m *big.Int, y4 [4]*big.Int) [4]*big.Int {
	// make sure x > 1, m is not nil, and m > 0, otherwise, use default Exp function
	if x.Cmp(big1) <= 0 || m == nil || m.Sign() <= 0 {
		return defaultExp4(x, m, y4)
	}
	// make sure all the y4 elements are positive
	for i := range y4 {
		if y4[i].Sign() <= 0 {
			return defaultExp4(x, m, y4)
		}
	}
	// make sure m is odd
	if m.Bit(0) != 1 {
		return defaultExp4(x, m, y4)
	}
	xWords, mWords := newNat(x), newNat(m)
	return fourfoldExpNNMontgomery(xWords, mWords, y4)
}

// fourfoldExpNNMontgomery calculates x**y1 mod m and x**y2 mod m x**y3 mod m and x**y4 mod m
// Uses Montgomery representation.
func fourfoldExpNNMontgomery(x, m nat, y [4]*big.Int) [4]*big.Int {
	power0, power1, k0, numWords := montgomerySetup(x, m)
	// Zero round, find common bits of the four values
	//fmt.Println("test here, len = ", len([]nat{y[0].abs, y[1].abs, y[2].abs, y[3].abs}))
	gcwList := fourfoldGCW([4]nat{newNat(y[0]), newNat(y[1]), newNat(y[2]), newNat(y[3])})
	// First round, find common bits of the three values
	var cm012, cm013, cm023, cm123 nat
	cm012 = threefoldGCW(*(*[3]nat)(gcwList[:3]))
	cm013 = threefoldGCW([3]nat{gcwList[0], gcwList[1], gcwList[3]})
	cm023 = threefoldGCW([3]nat{gcwList[0], gcwList[2], gcwList[3]})
	cm123 = threefoldGCW(*(*[3]nat)(gcwList[1:4]))

	var cm01, cm23, cm02, cm13, cm03, cm12 nat
	gcwList[0], gcwList[1], cm01 = gcw(gcwList[0], gcwList[1])
	gcwList[2], gcwList[3], cm23 = gcw(gcwList[2], gcwList[3])
	gcwList[0], gcwList[2], cm02 = gcw(gcwList[0], gcwList[2])
	gcwList[1], gcwList[3], cm13 = gcw(gcwList[1], gcwList[3])
	gcwList[0], gcwList[3], cm03 = gcw(gcwList[0], gcwList[3])
	gcwList[1], gcwList[2], cm12 = gcw(gcwList[1], gcwList[2])

	z := multiMontgomery(m, power0, power1, k0, numWords,
		//      0-4      	  5     6      7       8     9     10     11    12    13    14
		append(gcwList[:], cm012, cm013, cm023, cm123, cm01, cm23, cm02, cm13, cm03, cm12),
	)

	// calculate the actual values
	var converted [4]nat
	converted[0] = assembleAndConvert(z[0], []nat{z[4], z[5], z[6], z[7], z[9], z[11], z[13]}, m, k0, numWords)
	converted[1] = assembleAndConvert(z[1], []nat{z[4], z[5], z[6], z[8], z[9], z[12], z[14]}, m, k0, numWords)
	converted[2] = assembleAndConvert(z[2], []nat{z[4], z[5], z[7], z[8], z[10], z[11], z[14]}, m, k0, numWords)
	converted[3] = assembleAndConvert(z[3], []nat{z[4], z[6], z[7], z[8], z[10], z[12], z[13]}, m, k0, numWords)

	var ret [4]*big.Int
	// normalize and set value
	for i := range ret {
		converted[i].norm()
		ret[i] = new(big.Int).SetBits(converted[i].intBits())
	}
	return ret
}

// ExpParallel computes x ** y mod |m| utilizing multiple CPU cores
// numRoutine specifies the number of routine for computing the result
func ExpParallel(x, y, m *big.Int, preTable *PreTable, numRoutine, wordChunkSize int) *big.Int {
	if preTable == nil {
		panic("precompute table is nil")
	}
	if preTable.Base.Cmp(x) != 0 {
		panic("precompute table not match: invalid base")
	}
	if preTable.Modulus.Cmp(m) != 0 {
		panic("precompute table not match: invalid modulus")
	}
	// make sure x > 1, m is not nil, m > 0, m is odd, and y is positive,
	// otherwise, use default Exp function
	if x.Cmp(big1) <= 0 || y.Sign() <= 0 || m == nil || m.Sign() <= 0 || m.Bit(0) != 1 {
		return new(big.Int).Exp(x, y, m)
	}
	if numRoutine <= 0 {
		numRoutine = 1
	}
	if wordChunkSize <= 0 {
		wordChunkSize = defaultWordChunkSize
	}
	xWords, yWords, mWords := newNat(x), newNat(y), newNat(m)
	zWords := expNNMontgomeryPrecomputedParallel(xWords, yWords, mWords, preTable, numRoutine, wordChunkSize)
	return new(big.Int).SetBits(zWords.intBits())
}

func expNNMontgomeryPrecomputedParallel(x, y, m nat, table *PreTable, numRoutines, wordChunkSize int) nat {
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()
	power0, _, k0, numWords := montgomerySetup(x, m)
	numPivots := len(y) / wordChunkSize
	if len(y)%wordChunkSize != 0 {
		numPivots++
	}
	pivots := make(chan int, numPivots)
	for i := 0; i < len(y); i += wordChunkSize {
		pivots <- i
	}
	outputs := make(chan nat, numRoutines)
	for i := 0; i < numRoutines; i++ {
		go table.routineExpNNMontgomery(ctx, power0, y, m, k0, wordChunkSize, pivots, outputs)
	}
	ret := power0
	temp := nat(nil).make(numWords)
	for out := range outputs {
		temp = temp.montgomery(ret, out, m, k0, numWords)
		ret, temp = temp, ret
		numRoutines--
		if numRoutines == 0 {
			close(pivots)
			close(outputs)
			break
		}
	}
	one := make(nat, numWords)
	one[0] = 1
	temp = temp.montgomery(ret, one, m, k0, numWords)
	ret, temp = temp, ret
	// final reduction
	if ret.cmp(m) >= 0 {
		ret = ret.sub(ret, m)
		if ret.cmp(m) >= 0 {
			_, ret = nat(nil).div(nil, ret, m)
		}
	}
	// normalization
	return ret.norm()
}
