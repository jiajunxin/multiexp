package multiexp

import . "math/big"

// GCB inputs two positive integer a and b, calculates the greatest common
func gcb(a, b nat) (nat, nat, nat) {
	var minBitLen int
	if len(a) > len(b) {
		minBitLen = len(b)
	} else {
		minBitLen = len(a)
	}
	var aNew, bNew, c nat

	aNew = aNew.make(len(a))
	bNew = bNew.make(len(b))
	c = c.make(minBitLen)
	for i := 0; i < minBitLen; i++ {
		c[i] = commonBits(a[i], b[i])
		aNew[i] = a[i] - c[i]
		bNew[i] = b[i] - c[i]
	}
	if len(a) > len(b) {
		for i := minBitLen; i < len(a); i++ {
			aNew[i] = a[i]
		}
	} else {
		for i := minBitLen; i < len(b); i++ {
			bNew[i] = b[i]
		}
	}
	return aNew, bNew, c
}

// fourfoldGcb inputs four positive integer a, b, c, d and calculates the greatest common
// the last element in output is the common bits integer
func fourfoldGcb(input []nat) []nat {
	if len(input) != 4 {
		panic("fourfoldGcb require the input size to be 4")
	}
	maxBitLen := 0
	minBitLen := len(input[0])
	for i := 0; i < 4; i++ {
		if maxBitLen < len(input[i]) {
			maxBitLen = len(input[i])
		}
		if minBitLen > len(input[i]) {
			minBitLen = len(input[i])
		}
	}

	var output [5]nat
	for i := 0; i < 4; i++ {
		output[i] = output[i].make(len(input[i]))
	}
	output[4] = output[4].make(minBitLen)
	for i := 0; i < minBitLen; i++ {
		output[4][i] = fourfoldCommonBits(input[0][i], input[1][i], input[2][i], input[3][i])
		output[0][i] = input[0][i] - output[4][i]
		output[1][i] = input[1][i] - output[4][i]
		output[2][i] = input[2][i] - output[4][i]
		output[3][i] = input[3][i] - output[4][i]
	}
	for i := 0; i < 4; i++ {
		if len(output[i]) > minBitLen {
			for j := minBitLen; j < len(output[i]); j++ {
				output[i][j] = input[i][j]
			}
		}
	}

	return output[:]
}

// threefoldGcb inputs three positive integer a, b, c and calculates the greatest common
// the last element in output is the common bits integer
func threefoldGcb(input []nat) nat {
	if len(input) != 3 {
		panic("threefoldGcb require the input size to be 3")
	}
	maxBitLen := 0
	minBitLen := len(input[0])
	for i := 0; i < 3; i++ {
		if maxBitLen < len(input[i]) {
			maxBitLen = len(input[i])
		}
		if minBitLen > len(input[i]) {
			minBitLen = len(input[i])
		}
	}

	var output nat
	output = output.make(minBitLen)
	for i := 0; i < minBitLen; i++ {
		output[i] = threefoldCommonBits(input[0][i], input[1][i], input[2][i])
		input[0][i] = input[0][i] - output[i]
		input[1][i] = input[1][i] - output[i]
		input[2][i] = input[2][i] - output[i]
	}
	return output[:]
}

func commonBits(a, b Word) Word {
	var ret uint
	ret = 0
	var mask uint
	for i := 0; i < _W; i++ {
		mask = uint(1 << i)
		if ((uint(a) & mask) == mask) && ((uint(b) & mask) == mask) {
			//fmt.Println("i == ", i, "mask = ", mask)
			ret = uint(ret) | mask
		}
	}
	return Word(ret)
}

func threefoldCommonBits(a, b, c Word) Word {
	var ret uint
	ret = 0
	var mask uint
	for i := 0; i < _W; i++ {
		mask = uint(1 << i)
		if ((uint(a) & mask) == mask) && ((uint(b) & mask) == mask) && ((uint(c) & mask) == mask) {
			ret = uint(ret) | mask
		}
	}
	return Word(ret)
}

func fourfoldCommonBits(a, b, c, d Word) Word {
	var ret uint
	ret = 0
	var mask uint
	for i := 0; i < _W; i++ {
		mask = uint(1 << i)
		if ((uint(a) & mask) == mask) && ((uint(b) & mask) == mask) && ((uint(c) & mask) == mask) && ((uint(d) & mask) == mask) {
			ret = uint(ret) | mask
		}
	}
	return Word(ret)
}
