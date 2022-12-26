package multiexp

import (
	"math/big"
)

var masks = [_W]big.Word{}

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
func DoubleExp(x, y1, y2, m *big.Int) []*big.Int {
	xWords := x.Bits()
	if len(xWords) == 0 {
		return ones(2)
	}
	if m == nil {
		return ones(2)
	}
	if x.Sign() <= 0 || y1.Sign() <= 0 || y2.Sign() <= 0 || m.Sign() <= 0 {
		return defaultExp(x, m, []*big.Int{y1, y2})
	}
	if len(xWords) == 0 || (len(xWords) == 1 && xWords[0] <= 1) {
		return ones(2)
	}
	// x > 1

	mWords := m.Bits() // m.abs may be nil for m == 0
	if len(mWords) == 0 {
		return ones(2)
	}
	// m > 1

	y1Words, y2Words := y1.Bits(), y2.Bits()
	// x**0 == 1 or x**1 == x
	if len(y1Words) == 0 || (len(y1Words) == 1 && y1Words[0] <= 1) ||
		len(y2Words) == 0 || (len(y2Words) == 1 && y2Words[0] <= 1) {
		return defaultExp(x, m, []*big.Int{y1, y2})
	}
	// y > 1

	// only consider odd number m
	if mWords[0]&1 == 1 {
		return doubleExpNNMontgomery(xWords, y1Words, y2Words, mWords)
	}

	// if m is even
	results := make([]*big.Int, 2)
	for idx := range results {
		results[idx] = new(big.Int)
		nat(results[idx].Bits()).rem(xWords, mWords)
	}
	return results
}

// ones inputs a slice length and returns a slice of *Int, with all values "1"
func ones(length int) []*big.Int {
	ret := make([]*big.Int, length)
	for i := range ret {
		ret[i] = new(big.Int).SetInt64(1)
	}
	return ret
}

// defaultExp uses the default Exp function of big int to handle the edge cases that cannot benefit from this library
// in terms of performance
func defaultExp(x, m *big.Int, yList []*big.Int) []*big.Int {
	ret := make([]*big.Int, len(yList))
	for i := range yList {
		ret[i] = new(big.Int).Exp(x, yList[i], m)
	}
	return ret
}

// doubleExpNNMontgomery calculates x**y1 mod m and x**y2 mod m
// Uses Montgomery representation.
func doubleExpNNMontgomery(x, y1, y2, m nat) []*big.Int {
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

	ret := make([]*big.Int, 2)
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
		ret[i] = new(big.Int).SetBits(mmValues[i])
	}

	return ret
}

func montgomerySetup(x, m nat) (power0, power1 nat, k0 big.Word, numWords int) {
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
func multiMontgomery(m, power0, power1 nat, k0 big.Word, numWords int, yList []nat) []nat {
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

// multiMontgomeryWithPreComputeTable calculates the modular montgomery exponent with result not normalized
func multiMontgomeryWithPreComputeTable(m, power0, power1 nat, k0 big.Word, numWords int, yList []nat, preTable *PreTable) []nat {
	// initialize each value to be 1 (Montgomery 1)
	z := make([]nat, len(yList))
	for i := range z {
		z[i] = z[i].make(numWords)
		copy(z[i], power0)
	}

	var squaredPower, temp nat
	squaredPower = squaredPower.make(numWords)
	temp = temp.make(numWords)
	copy(squaredPower, power1)
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
func FourfoldExp(x, m *big.Int, y []*big.Int) []*big.Int {
	xWords := x.Bits()
	if len(xWords) == 0 {
		return ones(4)
	}
	if x.Sign() <= 0 || m.Sign() <= 0 {
		return defaultExp(x, m, y)
	}
	if len(y)%2 != 0 {
		return defaultExp(x, m, y)
	}
	for i := range y {
		if y[i].Sign() <= 0 {
			return defaultExp(x, m, y)
		}
	}
	if len(xWords) == 1 && xWords[0] == 1 {
		return ones(len(y))
	}

	// x > 1

	if m == nil {
		return ones(len(y))
	}
	mWords := m.Bits() // m.abs may be nil for m == 0
	if len(mWords) == 0 {
		return ones(len(y))
	}
	// m > 1
	// y > 0

	if mWords[0]&1 != 1 {
		return defaultExp(x, m, y)
	}
	return fourfoldExpNNMontgomery(xWords, mWords, y)
}

// fourfoldExpNNMontgomery calculates x**y1 mod m and x**y2 mod m x**y3 mod m and x**y4 mod m
// Uses Montgomery representation.
func fourfoldExpNNMontgomery(x, m nat, y []*big.Int) []*big.Int {
	power0, power1, k0, numWords := montgomerySetup(x, m)
	// Zero round, find common bits of the four values
	//fmt.Println("test here, len = ", len([]nat{y[0].abs, y[1].abs, y[2].abs, y[3].abs}))
	gcwList := fourfoldGCW([]nat{y[0].Bits(), y[1].Bits(), y[2].Bits(), y[3].Bits()})
	// First round, find common bits of the three values
	var cm012, cm013, cm023, cm123 nat
	cm012 = threefoldGCW(gcwList[:3])
	cm013 = threefoldGCW([]nat{gcwList[0], gcwList[1], gcwList[3]})
	cm023 = threefoldGCW([]nat{gcwList[0], gcwList[2], gcwList[3]})
	cm123 = threefoldGCW(gcwList[1:4])

	var cm01, cm23, cm02, cm13, cm03, cm12 nat
	gcwList[0], gcwList[1], cm01 = gcw(gcwList[0], gcwList[1])
	gcwList[2], gcwList[3], cm23 = gcw(gcwList[2], gcwList[3])
	gcwList[0], gcwList[2], cm02 = gcw(gcwList[0], gcwList[2])
	gcwList[1], gcwList[3], cm13 = gcw(gcwList[1], gcwList[3])
	gcwList[0], gcwList[3], cm03 = gcw(gcwList[0], gcwList[3])
	gcwList[1], gcwList[2], cm12 = gcw(gcwList[1], gcwList[2])
	//                                                                    0-4	  5     6      7       8     9     10     11    12    13    14
	z := multiMontgomery(m, power0, power1, k0, numWords, append(gcwList, cm012, cm013, cm023, cm123, cm01, cm23, cm02, cm13, cm03, cm12))

	// calculate the actual values
	assembleAndConvert(&z[0], []nat{z[4], z[5], z[6], z[7], z[9], z[11], z[13]}, m, k0, numWords)
	assembleAndConvert(&z[1], []nat{z[4], z[5], z[6], z[8], z[9], z[12], z[14]}, m, k0, numWords)
	assembleAndConvert(&z[2], []nat{z[4], z[5], z[7], z[8], z[10], z[11], z[14]}, m, k0, numWords)
	assembleAndConvert(&z[3], []nat{z[4], z[6], z[7], z[8], z[10], z[12], z[13]}, m, k0, numWords)

	z = z[:4] //the rest are useless now

	ret := make([]*big.Int, 4)
	// normalize and set value
	for i := range z {
		z[i].norm()
		ret[i] = new(big.Int).SetBits(z[i])
	}
	return ret
}
