package multiexp

import (
	. "math/big"
)

// defaultExp handles the conditions this lib do not support or cannot do better than the golang's exp
func defaultExp(x, m *Int, yList []*Int) []*Int {
	ret := make([]*Int, len(yList))
	for i := range yList {
		ret[i] = new(Int).Exp(x, yList[i], m)
	}
	return ret
}

// DoubleExp sets z1 = x**y1 mod |m|, z2 = x**y2 mod |m| ... (i.e. the sign of m is ignored), and returns z1, z2.
// If m == nil or m == 0, z = x**y unless y <= 0 then z = 1. If m != 0, y < 0,
// and x and m are not relatively prime, z is unchanged and nil is returned.
//
// DoubleExp is not a cryptographically constant-time operation.
func DoubleExp(x, y1, y2, m *Int) []*Int {
	xWords := x.Bits()
	if len(xWords) == 0 {
		return ones(2)
	}
	if m == nil {
		return ones(2)
	}
	if x.Sign() <= 0 || y1.Sign() <= 0 || y2.Sign() <= 0 || m.Sign() <= 0 {
		return defaultExp(x, m, []*Int{y1, y2})
	}
	if len(xWords) <= 1 {
		return ones(2)
	}
	// x > 1

	mWords := m.Bits() // m.abs may be nil for m == 0
	if len(mWords) == 0 {
		return ones(2)
	}
	// m > 1

	y1Words := y1.Bits()
	y2Words := y2.Bits()
	// x**0 == 1 or x**1 == x
	if len(y1Words) <= 1 || len(y2Words) <= 1 {
		return defaultExp(x, m, []*Int{y1, y2})
	}
	// y > 1

	// only consider odd number m
	if mWords[0] == 1 {
		return doubleExpNNMontgomery(xWords, y1Words, y2Words, mWords)
	}

	// if m is even
	results := make([]*Int, 2)
	for idx := range results {
		results[idx] = new(Int)
		nat(results[idx].Bits()).rem(xWords, mWords)
	}
	return results
}

// ones inputs a slice length and returns a slice of *Int, with all values "1"
func ones(length int) []*Int {
	ret := make([]*Int, length)
	for i := range ret {
		ret[i] = new(Int).SetInt64(1)
	}
	return ret
}

// doubleExpNNMontgomery calculates x**y1 mod m and x**y2 mod m
// Uses Montgomery representation.
func doubleExpNNMontgomery(x, y1, y2, m nat) []*Int {
	numWords := len(m)

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
	k0 := 2 - m[0]
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

	// powers[i] contains x^i
	//var powers [2]nat
	//powers[0] = powers[0].montgomery(one, RR, m, k0, numWords)
	//powers[1] = powers[1].montgomery(x, RR, m, k0, numWords)
	var power0, power1 nat
	power0 = power0.montgomery(one, RR, m, k0, numWords)
	power1 = power1.montgomery(x, RR, m, k0, numWords)

	y1New, y2New, y3 := gcb(y1, y2)
	z := multiMontgomery(m, power0, power1, k0, numWords, []nat{y1New, y2New, y3})
	// calculate z1 and z2
	temp := nat(nil).make(numWords)
	temp = temp.montgomery(z[0], z[2], m, k0, numWords)
	z[0], temp = temp, z[0]
	temp = temp.montgomery(z[1], z[2], m, k0, numWords)
	z[1], temp = temp, z[1]
	z = z[:2] //z3 is useless now
	// convert to regular number
	for i := range z {
		temp = temp.montgomery(z[i], one, m, k0, numWords)
		z[i], temp = temp, z[i]
	}
	for i := range z {
		// One last reduction, just in case.
		// See golang.org/issue/13907.
		if z[i].cmp(m) >= 0 {
			// Common case is m has high bit set; in that case,
			// since zz is the same length as m, there can be just
			// one multiple of m to remove. Just subtract.
			// We think that the subtract should be sufficient in general,
			// so do that unconditionally, but double-check,
			// in case our beliefs are wrong.
			// The div is not expected to be reached.
			z[i] = z[i].sub(z[i], m)
			if z[i].cmp(m) >= 0 {
				_, z[i] = nat(nil).div(nil, z[i], m)
			}
		}
	}

	ret := make([]*Int, 2)
	for i := range ret {
		ret[i] = new(Int)
	}

	// normalize and set value
	for i := range z {
		z[i].norm()
		ret[i].SetBits(z[i])
	}
	return ret
}

// multiMontgomery calculates the modular montgomery exponent with result not normalized
func multiMontgomery(m, power0, power1 nat, k0 Word, numWords int, yList []nat) []nat {
	// initialize each value to be 1 (Montgomery 1)
	z := make([]nat, len(yList))
	for i := range z {
		z[i] = z[i].make(numWords)
		copy(z[i], power0)
	}

	var squaredPower nat
	squaredPower = squaredPower.make(numWords)
	copy(squaredPower, power1)
	//	fmt.Println("squaredPower = ", squaredPower.String())

	maxLen := 1
	for i := range yList {
		if len(yList[i]) > maxLen {
			maxLen = len(yList[i])
		}
	}

	temp := nat(nil).make(numWords)
	for i := 0; i < maxLen; i++ {
		for j := 0; j < _W; j++ {
			mask := Word(1 << j)
			for k := range yList {
				if len(yList[k]) > i {
					if (yList[k][i] & mask) == mask {
						temp = temp.montgomery(z[k], squaredPower, m, k0, numWords)
						z[k], temp = temp, z[k]
					}
				}
			}
			// montgomery must have the returned value not same as the input values
			// we have to use this temp as the middle variable
			temp = temp.montgomery(squaredPower, squaredPower, m, k0, numWords)
			squaredPower, temp = temp, squaredPower
		}
	}
	return z
}

// multiMontgomeryWithPreComputeTable calculates the modular montgomery exponent with result not normalized
func multiMontgomeryWithPreComputeTable(m, power0, power1 nat, k0 Word, numWords int, y []nat, pretable *PreTable) []nat {
	// initialize each value to be 1 (Montgomery 1)
	z := make([]nat, len(y))
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
	for i := range y {
		if len(y[i]) > maxLen {
			maxLen = len(y[i])
		}
	}

	for i := 0; i < maxLen; i++ {
		for j := 0; j < _W; j++ {
			mask := Word(1 << j)
			for k := range y {
				if len(y[k]) > i {
					if (y[k][i] & mask) == mask {
						temp = temp.montgomery(z[k], pretable.table[i][j], m, k0, numWords)
						z[k], temp = temp, z[k]
					}
				}
			}
			// // montgomery must have the returned value not same as the input values
			// // we have to use this temp as the middle variable
			// temp = temp.montgomery(squaredPower, squaredPower, m, k0, numWords)
			// squaredPower, temp = temp, squaredPower
		}
	}
	return z
}

// FourfoldExp sets z1 = x**y1 mod |m|, z2 = x**y2 mod |m| ... (i.e. the sign of m is ignored), and returns z1, z2...
// In construction, many panic conditions. Use at your own risk!
//
// FourfoldExp is not a cryptographically constant-time operation.
func FourfoldExp(x, m *Int, y []*Int) []*Int {
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
func fourfoldExpNNMontgomery(x, m nat, y []*Int) []*Int {
	numWords := len(m)

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
	k0 := 2 - m[0]
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

	// powers[i] contains x^i
	var powers [2]nat
	powers[0] = powers[0].montgomery(one, RR, m, k0, numWords)
	powers[1] = powers[1].montgomery(x, RR, m, k0, numWords)

	// Zero round, find common bits of the four values
	//fmt.Println("test here, len = ", len([]nat{y[0].abs, y[1].abs, y[2].abs, y[3].abs}))
	yNewList := fourfoldGcb([]nat{y[0].Bits(), y[1].Bits(), y[2].Bits(), y[3].Bits()})
	// First round, find common bits of the three values
	var cm012, cm013, cm023, cm123 nat
	cm012 = threefoldGcb(yNewList[:3])
	cm013 = threefoldGcb([]nat{yNewList[0], yNewList[1], yNewList[3]})
	cm023 = threefoldGcb([]nat{yNewList[0], yNewList[2], yNewList[3]})
	cm123 = threefoldGcb(yNewList[1:4])

	var cm01, cm23, cm02, cm13, cm03, cm12 nat
	yNewList[0], yNewList[1], cm01 = gcb(yNewList[0], yNewList[1])
	yNewList[2], yNewList[3], cm23 = gcb(yNewList[2], yNewList[3])
	yNewList[0], yNewList[2], cm02 = gcb(yNewList[0], yNewList[2])
	yNewList[1], yNewList[3], cm13 = gcb(yNewList[1], yNewList[3])
	yNewList[0], yNewList[3], cm03 = gcb(yNewList[0], yNewList[3])
	yNewList[1], yNewList[2], cm12 = gcb(yNewList[1], yNewList[2])
	//                                                                    0-4	  5     6      7       8     9     10     11    12    13    14
	z := multiMontgomery(m, powers[0], powers[1], k0, numWords, append(yNewList, cm012, cm013, cm023, cm123, cm01, cm23, cm02, cm13, cm03, cm12))

	// calculate the actual values
	assembleAndConvert(&z[0], []nat{z[4], z[5], z[6], z[7], z[9], z[11], z[13]}, m, k0, numWords)
	assembleAndConvert(&z[1], []nat{z[4], z[5], z[6], z[8], z[9], z[12], z[14]}, m, k0, numWords)
	assembleAndConvert(&z[2], []nat{z[4], z[5], z[7], z[8], z[10], z[11], z[14]}, m, k0, numWords)
	assembleAndConvert(&z[3], []nat{z[4], z[6], z[7], z[8], z[10], z[12], z[13]}, m, k0, numWords)

	z = z[:4] //the rest are useless now
	ret := make([]*Int, 4)
	for i := range ret {
		ret[i] = new(Int)
	}

	// normalize and set value
	for i := range z {
		z[i].norm()
		ret[i].SetBits(z[i])
	}
	return ret
}
