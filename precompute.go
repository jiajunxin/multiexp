package multiexp

import (
	"context"
	"fmt"
	"sync"
	"time"

	"math/big"
)

type PreTable struct {
	Base          *big.Int
	Modulus       *big.Int
	TableSize     int
	table         [][]nat
	wordChunkSize int
}

func NewPrecomputeTable(base, modular *big.Int, tableSize int) *PreTable {
	if tableSize <= 0 {
		return nil
	}
	if base == nil || modular == nil {
		return nil
	}
	if base.Sign() <= 0 || modular.Sign() <= 0 {
		return nil
	}

	x := base.Bits()
	if len(x) == 0 {
		return nil
	}
	if len(x) == 1 && x[0] == 1 {
		return nil
	}
	// x > 1

	m := modular.Bits() // m.abs may be nil for m == 0
	_, power1, k0, numWords := montgomerySetup(x, m)
	if numWords == 0 {
		return nil
	}

	var temp, squaredPower nat
	temp = temp.make(numWords)
	squaredPower = squaredPower.make(numWords)
	copy(squaredPower, power1)
	preTable := make([][]nat, tableSize)
	for i := range preTable {
		preTable[i] = make([]nat, _W)
		for j := range preTable[i] {
			preTable[i][j] = preTable[i][j].make(numWords)
		}
	}

	for i := 0; i < tableSize; i++ {
		for j := 0; j < _W; j++ {
			// montgomery must have the returned value not same as the input values
			// we have to use this temp as the middle variable
			copy(preTable[i][j], squaredPower)
			temp = temp.montgomery(squaredPower, squaredPower, m, k0, numWords)
			squaredPower, temp = temp, squaredPower
		}
	}

	return &PreTable{
		Base:      base,
		Modulus:   modular,
		TableSize: tableSize,
		table:     preTable,
	}
}

func (p *PreTable) routineExpNNMontgomery(ctx context.Context, power0, m nat, k0 big.Word, numWords int,
	inputChan chan input, outputChan chan nat) {
	temp := nat(nil).make(numWords)
	for {
		select {
		case <-ctx.Done():
			return
		case in := <-inputChan:
			ret := nat(nil).make(numWords)
			copy(ret, power0)
			for i := 0; i < in.wordLen; i++ {
				for j := 0; j < _W; j++ {
					if (in.wordChunk[i] & masks[j]) != masks[j] {
						continue
					}
					temp = temp.montgomery(ret, p.table[i+in.pivot][j], m, k0, numWords)
					ret, temp = temp, ret
				}
			}
			outputChan <- ret
		default:
			continue
		}
	}
}

// FourfoldExpWithPreComputeTableParallel sets z1 = x**y1 mod |m|, z2 = x**y2 mod |m| ... (i.e. the sign of m is ignored), and returns z1, z2...
// In construction, many panic conditions. Use at your own risk!
// Use at most 4 threads for now.
// FourfoldExpWithPreComputeTableParallel is not a cryptographically constant-time operation.
func FourfoldExpWithPreComputeTableParallel(x, m *big.Int, yList []*big.Int, preTable *PreTable) []*big.Int {
	if x.Sign() < 0 {
		panic("invalid x: negative value")
	}
	if x.Cmp(big1) <= 0 {
		return defaultExp(x, m, yList)
	}
	if m == nil {
		panic("invalid m: nil value")
	}
	if m.Sign() <= 0 {
		panic("invalid m: non-positive value")
	}
	if len(yList) != 4 {
		panic("invalid yList element size")
	}
	for i := range yList {
		if yList[i].Sign() <= 0 {
			panic("invalid yList: non-positive value")
		}
	}
	if m.Bit(0) != 1 {
		panic("The input modular is not an odd number")
	}
	// check if the table is same as the input parameters
	if preTable.Base.Cmp(x) != 0 || preTable.Modulus.Cmp(m) != 0 {
		panic("The input table does not match the input")
	}
	xWords, mWords := x.Bits(), m.Bits()
	return fourfoldExpNNMontgomeryWithPreComputeTableParallel(xWords, mWords, yList, preTable)
}

// fourfoldExpNNMontgomery calculates x**y1 mod m and x**y2 mod m x**y3 mod m and x**y4 mod m
// Uses Montgomery representation.
func fourfoldExpNNMontgomeryWithPreComputeTableParallel(x, m nat, y []*big.Int, preTable *PreTable) []*big.Int {
	power0, _, k0, numWords := montgomerySetup(x, m)

	gcwList := fourfoldGCW([]nat{y[0].Bits(), y[1].Bits(), y[2].Bits(), y[3].Bits()})

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
	c0 := make(chan []nat)
	c1 := make(chan []nat)
	c2 := make(chan []nat)
	c3 := make(chan []nat)
	go multiMontgomeryWithPreComputeTableWithChan(m, power0, k0, numWords, gcwList[:4], preTable, c0)
	go multiMontgomeryWithPreComputeTableWithChan(m, power0, k0, numWords, []nat{gcwList[4], cm012, cm013, cm023}, preTable, c1)
	go multiMontgomeryWithPreComputeTableWithChan(m, power0, k0, numWords, []nat{cm123, cm01, cm23, cm02}, preTable, c2)
	go multiMontgomeryWithPreComputeTableWithChan(m, power0, k0, numWords, []nat{cm13, cm03, cm12}, preTable, c3)

	z1 := <-c0
	z2 := <-c1
	z3 := <-c2
	z4 := <-c3
	z := append(z1, z2...)
	z = append(z, z3...)
	z = append(z, z4...)
	//                                                                    0-4	  5     6      7       8     9     10     11    12    13    14
	//z := multiMontgomeryWithPreComputeTable(RR, m, powers[0], powers[1], k0, numWords, append(gcwList, cm012, cm013, cm023, cm123, cm01, cm23, cm02, cm13, cm03, cm12), preTable)
	// calculate the actual values

	var wg sync.WaitGroup
	wg.Add(4)
	go assembleAndConvertWithWG(&z[0], []nat{z[4], z[5], z[6], z[7], z[9], z[11], z[13]}, m, k0, numWords, &wg)
	go assembleAndConvertWithWG(&z[1], []nat{z[4], z[5], z[6], z[8], z[9], z[12], z[14]}, m, k0, numWords, &wg)
	go assembleAndConvertWithWG(&z[2], []nat{z[4], z[5], z[7], z[8], z[10], z[11], z[14]}, m, k0, numWords, &wg)
	go assembleAndConvertWithWG(&z[3], []nat{z[4], z[6], z[7], z[8], z[10], z[12], z[13]}, m, k0, numWords, &wg)
	wg.Wait()

	z = z[:4]

	ret := make([]*big.Int, 4)
	// normalize and set value
	for i := range z {
		z[i].norm()
		ret[i] = new(big.Int).SetBits(z[i])
	}
	return ret
}

func assembleAndConvert(prod *nat, set []nat, mm nat, k0 big.Word, numWords int) {
	var temp nat
	temp = temp.make(numWords)
	var m nat
	m = m.make(numWords)
	copy(m, mm)
	for i := range set {
		temp = temp.montgomery(*prod, set[i], m, k0, numWords)

		*prod, temp = temp, *prod
		//fmt.Println("prod", i, " = ", prod.String())
	}

	// one = 1, with equal length to that of m
	one := make(nat, numWords)
	one[0] = 1
	// convert to regular number
	temp = temp.montgomery(*prod, one, m, k0, numWords)
	*prod, temp = temp, *prod
	//fmt.Println("prod convert = ", prod.String())
	// One last reduction, just in case.
	// See golang.org/issue/13907.
	if prod.cmp(m) >= 0 {
		// Common case is m has high bit set; in that case,
		// since zz is the same length as m, there can be just
		// one multiple of m to remove. Just subtract.
		// We think that the subtraction should be sufficient in general,
		// so do that unconditionally, but double-check,
		// in case our beliefs are wrong.
		// The div is not expected to be reached.
		*prod = (*prod).sub(*prod, m)
		if prod.cmp(m) >= 0 {
			_, *prod = nat(nil).div(nil, *prod, m)
		}
	}
}

func assembleAndConvertWithWG(prod *nat, set []nat, mm nat, k0 big.Word, numWords int, wg *sync.WaitGroup) {
	var temp nat
	temp = temp.make(numWords)
	var m nat
	m = m.make(numWords)
	copy(m, mm)
	for i := range set {
		temp = temp.montgomery(*prod, set[i], m, k0, numWords)

		*prod, temp = temp, *prod
		//fmt.Println("prod", i, " = ", prod.String())
	}

	// one = 1, with equal length to that of m
	one := make(nat, numWords)
	one[0] = 1
	// convert to regular number
	temp = temp.montgomery(*prod, one, m, k0, numWords)
	*prod, temp = temp, *prod
	//fmt.Println("prod convert = ", prod.String())
	// One last reduction, just in case.
	// See golang.org/issue/13907.
	if prod.cmp(m) >= 0 {
		// Common case is m has high bit set; in that case,
		// since zz is the same length as m, there can be just
		// one multiple of m to remove. Just subtract.
		// We think that the subtraction should be sufficient in general,
		// so do that unconditionally, but double-check,
		// in case our beliefs are wrong.
		// The div is not expected to be reached.
		*prod = (*prod).sub(*prod, m)
		if prod.cmp(m) >= 0 {
			_, *prod = nat(nil).div(nil, *prod, m)
		}
	}
	wg.Done()
}

// multiMontgomeryWithPreComputeTableWithChan calculates the modular montgomery exponent with result not normalized
func multiMontgomeryWithPreComputeTableWithChan(m, power0 nat, k0 big.Word, numWords int,
	y []nat, preTable *PreTable, c chan []nat) {
	startingTime := time.Now().UTC()

	// initialize each value to be 1 (Montgomery 1)
	z := make([]nat, len(y))
	for i := range z {
		z[i] = z[i].make(numWords)
		copy(z[i], power0)
	}

	maxLen := 1
	for i := range y {
		if len(y[i]) > maxLen {
			maxLen = len(y[i])
		}
	}

	temp := nat(nil).make(numWords)
	for i := 0; i < maxLen; i++ {
		for j := 0; j < _W; j++ {
			//mask := Word(1 << j)
			for k := range y {
				if len(y[k]) <= i {
					continue
				}
				if (y[k][i] & masks[j]) != masks[j] {
					continue
				}
				temp = temp.montgomery(z[k], preTable.table[i][j], m, k0, numWords)
				z[k], temp = temp, z[k]
			}
			// // montgomery must have the returned value not same as the input values
			// // we have to use this temp as the middle variable
			// temp = temp.montgomery(squaredPower, squaredPower, m, k0, numWords)
			// squaredPower, temp = temp, squaredPower
		}
	}
	duration := time.Now().UTC().Sub(startingTime)
	fmt.Println("inside multiMontgomeryWithPreComputeTableWithChan, len(y) = ", len(y))
	fmt.Printf("Running multiMontgomeryWithPreComputeTableWithChan Takes [%.3f] Seconds \n", duration.Seconds())
	c <- z
}
