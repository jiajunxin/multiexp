package multiexp

import (
	"context"
	"fmt"

	"math/big"
	"math/bits"
)

// PreTable is the pre-computation table for multi-exponentiation
type PreTable struct {
	Base      *big.Int
	Modulus   *big.Int
	TableSize int
	table     [][_W]nat
}

func GetTableSize(table *PreTable) {
	fmt.Println("The table size = ", table.TableSize, "rows, ", bits.UintSize, " columns, each element size = ", bits.UintSize)
	fmt.Println("Totally ", table.TableSize*bits.UintSize*bits.UintSize/8, "bytes")
}

// NewPrecomputeTable creates a pre-computation table for multi-exponentiation
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

	x := newNat(base)
	if len(x) == 0 {
		return nil
	}
	if len(x) == 1 && x[0] == 1 {
		return nil
	}
	// x > 1

	m := newNat(modular) // m.abs may be nil for m == 0
	_, power1, k0, numWords := montgomerySetup(x, m)
	if numWords == 0 {
		return nil
	}

	var temp, squaredPower nat
	temp = temp.make(numWords)
	squaredPower = squaredPower.make(numWords)
	copy(squaredPower, power1)
	preTable := make([][_W]nat, tableSize)
	for i := range preTable {
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

func (p *PreTable) routineExpNNMontgomery(ctx context.Context, power0, y, m nat, k0 Word, wordChunkSize int,
	pivots <-chan int, outputs chan<- nat) {
	numWords := len(m)
	ret := nat(nil).make(numWords)
	copy(ret, power0)
	temp := nat(nil).make(numWords)
	for {
		select {
		case <-ctx.Done():
			return
		case l := <-pivots:
			r := l + wordChunkSize
			if r > len(y) {
				r = len(y)
			}
			for i := l; i < r; i++ {
				for j := 0; j < _W; j++ {
					if (y[i] & masks[j]) != masks[j] {
						continue
					}
					temp = temp.montgomery(ret, p.table[i][j], m, k0, numWords)
					ret, temp = temp, ret
				}
			}
		default:
			if len(pivots) == 0 {
				outputs <- ret
				return
			}
		}
	}
}

// FourfoldExpPrecomputedParallel sets z1 = x**y1 mod |m|, z2 = x**y2 mod |m| ... (i.e. the sign of m is ignored), and returns z1, z2...
// In construction, many panic conditions. Use at your own risk!
// Use at most 4 threads for now.
// FourfoldExpPrecomputedParallel is not a cryptographically constant-time operation.
func FourfoldExpPrecomputedParallel(x, m *big.Int, y4 [4]*big.Int, preTable *PreTable) [4]*big.Int {
	if x.Sign() < 0 {
		panic("invalid x: negative value")
	}
	if x.Cmp(big1) <= 0 {
		return defaultExp4(x, m, y4)
	}
	if m == nil {
		panic("invalid m: nil value")
	}
	if m.Sign() <= 0 {
		panic("invalid m: non-positive value")
	}
	for i := range y4 {
		if y4[i].Sign() <= 0 {
			panic("invalid y4: non-positive value")
		}
	}
	if m.Bit(0) != 1 {
		panic("The input modular is not an odd number")
	}
	// check if the table is same as the input parameters
	if preTable.Base.Cmp(x) != 0 || preTable.Modulus.Cmp(m) != 0 {
		panic("The input table does not match the input")
	}
	xWords, mWords := newNat(x), newNat(m)
	return fourfoldExpNNMontgomeryPrecomputedParallel(xWords, mWords, y4, preTable)
}

// FourfoldExpPrecomputed sets z1 = x**y1 mod |m|, z2 = x**y2 mod |m| ... (i.e. the sign of m is ignored), and returns z1, z2...
// In construction, many panic conditions. Use at your own risk!
// Use single thread
// FourfoldExpPrecomputed is not a cryptographically constant-time operation.
func FourfoldExpPrecomputed(x, m *big.Int, y4 [4]*big.Int, preTable *PreTable) [4]*big.Int {
	if x.Sign() < 0 {
		panic("invalid x: negative value")
	}
	if x.Cmp(big1) <= 0 {
		return defaultExp4(x, m, y4)
	}
	if m == nil {
		panic("invalid m: nil value")
	}
	if m.Sign() <= 0 {
		panic("invalid m: non-positive value")
	}
	for i := range y4 {
		if y4[i].Sign() <= 0 {
			panic("invalid y4: non-positive value")
		}
	}
	if m.Bit(0) != 1 {
		panic("The input modular is not an odd number")
	}
	// check if the table is same as the input parameters
	if preTable.Base.Cmp(x) != 0 || preTable.Modulus.Cmp(m) != 0 {
		panic("The input table does not match the input")
	}
	xWords, mWords := newNat(x), newNat(m)
	return fourfoldExpNNMontgomeryPrecomputed(xWords, mWords, y4, preTable)
}

// fourfoldExpNNMontgomery calculates x**y1 mod m and x**y2 mod m x**y3 mod m and x**y4 mod m
// Uses Montgomery representation.
func fourfoldExpNNMontgomeryPrecomputedParallel(x, m nat, y4 [4]*big.Int, preTable *PreTable) [4]*big.Int {
	power0, _, k0, numWords := montgomerySetup(x, m)

	gcwList := fourfoldGCW([4]nat{newNat(y4[0]), newNat(y4[1]), newNat(y4[2]), newNat(y4[3])})

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
	var c4 [4]chan []nat
	for i := range c4 {
		c4[i] = make(chan []nat)
	}
	go multiMontgomeryPrecomputedChan(m, power0, k0, numWords, gcwList[:4], preTable, c4[0])
	go multiMontgomeryPrecomputedChan(m, power0, k0, numWords, []nat{gcwList[4], cm012, cm013, cm023}, preTable, c4[1])
	go multiMontgomeryPrecomputedChan(m, power0, k0, numWords, []nat{cm123, cm01, cm23, cm02}, preTable, c4[2])
	go multiMontgomeryPrecomputedChan(m, power0, k0, numWords, []nat{cm13, cm03, cm12}, preTable, c4[3])

	var z []nat
	for i := range c4 {
		z = append(z, <-c4[i]...)
	}
	// z := multiMontgomeryPrecomputed(RR, m, powers[0], powers[1], k0, numWords, append(gcwList, cm012, cm013, cm023, cm123, cm01, cm23, cm02, cm13, cm03, cm12), preTable)
	// calculate the actual values

	var outputs [4]chan nat
	for i := range outputs {
		outputs[i] = make(chan nat)
	}
	go assembleAndConvertChan(z[0], []nat{z[4], z[5], z[6], z[7], z[9], z[11], z[13]}, m, k0, numWords, outputs[0])
	go assembleAndConvertChan(z[1], []nat{z[4], z[5], z[6], z[8], z[9], z[12], z[14]}, m, k0, numWords, outputs[1])
	go assembleAndConvertChan(z[2], []nat{z[4], z[5], z[7], z[8], z[10], z[11], z[14]}, m, k0, numWords, outputs[2])
	go assembleAndConvertChan(z[3], []nat{z[4], z[6], z[7], z[8], z[10], z[12], z[13]}, m, k0, numWords, outputs[3])

	var ret [4]*big.Int
	// normalize and set value
	for i := range ret {
		output := <-outputs[i]
		output.norm()
		ret[i] = new(big.Int).SetBits(output.intBits())
	}
	return ret
}

// fourfoldExpNNMontgomery calculates x**y1 mod m and x**y2 mod m x**y3 mod m and x**y4 mod m
// Uses Montgomery representation.
func fourfoldExpNNMontgomeryPrecomputed(x, m nat, y4 [4]*big.Int, preTable *PreTable) [4]*big.Int {
	power0, _, k0, numWords := montgomerySetup(x, m)

	gcwList := fourfoldGCW([4]nat{newNat(y4[0]), newNat(y4[1]), newNat(y4[2]), newNat(y4[3])})

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
	// var c4 [4]chan []nat
	// for i := range c4 {
	// 	c4[i] = make(chan []nat)
	// }
	// multiMontgomeryPrecomputedChan(m, power0, k0, numWords, gcwList[:4], preTable, c4[0])
	// multiMontgomeryPrecomputedChan(m, power0, k0, numWords, []nat{gcwList[4], cm012, cm013, cm023}, preTable, c4[1])
	// multiMontgomeryPrecomputedChan(m, power0, k0, numWords, []nat{cm123, cm01, cm23, cm02}, preTable, c4[2])
	// multiMontgomeryPrecomputedChan(m, power0, k0, numWords, []nat{cm13, cm03, cm12}, preTable, c4[3])

	// var z []nat
	// for i := range c4 {
	// 	z = append(z, <-c4[i]...)
	// }
	z := multiMontgomeryPrecomputed(m, power0, k0, numWords, append(gcwList[:], cm012, cm013, cm023, cm123, cm01, cm23, cm02, cm13, cm03, cm12), preTable)
	// calculate the actual values

	var outputs [4]nat

	outputs[0] = assembleAndConvert(z[0], []nat{z[4], z[5], z[6], z[7], z[9], z[11], z[13]}, m, k0, numWords)
	outputs[1] = assembleAndConvert(z[1], []nat{z[4], z[5], z[6], z[8], z[9], z[12], z[14]}, m, k0, numWords)
	outputs[2] = assembleAndConvert(z[2], []nat{z[4], z[5], z[7], z[8], z[10], z[11], z[14]}, m, k0, numWords)
	outputs[3] = assembleAndConvert(z[3], []nat{z[4], z[6], z[7], z[8], z[10], z[12], z[13]}, m, k0, numWords)

	var ret [4]*big.Int
	// normalize and set value
	for i := range ret {
		output := outputs[i]
		output.norm()
		ret[i] = new(big.Int).SetBits(output.intBits())
	}
	return ret
}

func assembleAndConvert(prod nat, set []nat, mm nat, k0 Word, numWords int) nat {
	temp := nat(nil).make(numWords)
	m := nat(nil).make(numWords)
	copy(m, mm)
	for i := range set {
		temp = temp.montgomery(prod, set[i], m, k0, numWords)
		prod, temp = temp, prod
	}

	// one = 1, with equal length to that of m
	one := make(nat, numWords)
	one[0] = 1
	// convert to regular number
	temp = temp.montgomery(prod, one, m, k0, numWords)
	prod, temp = temp, prod
	// one last reduction, just in case.
	if prod.cmp(m) >= 0 {
		prod = prod.sub(prod, m)
		if prod.cmp(m) >= 0 {
			_, prod = nat(nil).div(nil, prod, m)
		}
	}
	return prod
}

func assembleAndConvertChan(prod nat, set []nat, mm nat, k0 Word, numWords int, output chan<- nat) {
	output <- assembleAndConvert(prod, set, mm, k0, numWords)
}

// multiMontgomeryPrecomputedChan calculates the modular montgomery exponent with result not normalized
func multiMontgomeryPrecomputedChan(m, power0 nat, k0 Word, numWords int,
	y []nat, preTable *PreTable, c chan []nat) {
	//startingTime := time.Now().UTC()

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
		}
	}
	//duration := time.Now().UTC().Sub(startingTime)
	// fmt.Println("inside multiMontgomeryPrecomputedChan, len(y) = ", len(y))
	// fmt.Printf("Running multiMontgomeryPrecomputedChan Takes [%.3f] Seconds \n", duration.Seconds())
	c <- z
}
