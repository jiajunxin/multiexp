package multiexp

import (
	"crypto/rand"
	"math/big"
	"sync"
	"testing"
)

const numTestBits = 20000
const numTestGroupBits = 2048

var (
	benchRandLimit          *big.Int
	benchRandGroupLimit     *big.Int
	onceBenchRandLimit      sync.Once
	onceBenchRandGroupLimit sync.Once
	g, mod                  *big.Int
	xList                   []*big.Int
	onceBenchParameters     sync.Once
	table                   *PreTable
	onceBenchTable          sync.Once
)

func getBenchRandLimit() *big.Int {
	onceBenchRandLimit.Do(func() {
		benchRandLimit = new(big.Int).SetInt64(1)
		benchRandLimit.Lsh(benchRandLimit, numTestBits)
	})
	return benchRandLimit
}

func getBenchGroupLimit() *big.Int {
	onceBenchRandGroupLimit.Do(func() {
		benchRandGroupLimit = new(big.Int).SetInt64(1)
		benchRandGroupLimit.Lsh(benchRandGroupLimit, numTestGroupBits)
	})
	return benchRandGroupLimit
}

// this is used to test different random g, mod and exp
// We separate it because we also need the static case for the precomputations
func getDifferentBenchParameters(numX int) []*big.Int {
	var xListRan []*big.Int
	for i := 0; i < 4; i++ {
		x := new(big.Int)
		x, _ = rand.Int(rand.Reader, getBenchRandLimit())
		xListRan = append(xListRan, x)
	}
	if numX < 0 || numX > len(xList) {
		numX = len(xList)
	}
	return xListRan[:numX]
}

func getBenchParameters(numX int) (*big.Int, *big.Int, []*big.Int) {
	onceBenchParameters.Do(func() {
		g, mod = new(big.Int), new(big.Int)
		g, _ = rand.Int(rand.Reader, getBenchGroupLimit())
		mod = getValidModulus(rand.Reader, getBenchGroupLimit())
		for i := 0; i < 4; i++ {
			x := new(big.Int)
			x, _ = rand.Int(rand.Reader, getBenchRandLimit())
			xList = append(xList, x)
		}
	})
	if numX < 0 || numX > len(xList) {
		numX = len(xList)
	}
	return g, mod, xList[:numX]
}

func getBenchPrecomputeTable() *PreTable {
	onceBenchTable.Do(func() {
		g, n, _ := getBenchParameters(0)
		randLmtLen := (getBenchRandLimit().BitLen() / _W) + 1
		table = NewPrecomputeTable(g, n, randLmtLen)
	})
	return table
}

func BenchmarkOriginalDoubleExp(b *testing.B) {
	g, n, _ := getBenchParameters(1)
	b.ResetTimer()
	var result big.Int
	for i := 0; i < b.N; i++ {
		xListRan := getDifferentBenchParameters(2)
		b.StartTimer()
		result.Exp(g, xListRan[0], n)
		result.Exp(g, xListRan[1], n)
		b.StopTimer()
	}
}

func BenchmarkDoubleExp(b *testing.B) {
	g, n, _ := getBenchParameters(1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		xListRan := getDifferentBenchParameters(2)
		x2 := (*[2]*big.Int)(xListRan)
		b.StartTimer()
		DoubleExp(g, *x2, n)
		b.StopTimer()
	}
}

func BenchmarkOriginalFourfoldExp(b *testing.B) {
	g, n, _ := getBenchParameters(1)
	b.ResetTimer()
	var result big.Int
	for i := 0; i < b.N; i++ {
		xListRan := getDifferentBenchParameters(4)
		b.StartTimer()
		result.Exp(g, xListRan[0], n)
		result.Exp(g, xListRan[1], n)
		result.Exp(g, xListRan[2], n)
		result.Exp(g, xListRan[3], n)
		b.StopTimer()
	}
}

func BenchmarkFourfoldExp(b *testing.B) {
	g, n, _ := getBenchParameters(1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		xListRan := getDifferentBenchParameters(4)
		x4 := (*[4]*big.Int)(xListRan)
		b.StartTimer()
		FourfoldExp(g, n, *x4)
		b.StopTimer()
	}
}

func BenchmarkFourfoldExpWithTable(b *testing.B) {
	g, n, _ := getBenchParameters(1)

	maxLen := (numTestBits / _W) + 1
	table := NewPrecomputeTable(g, n, maxLen)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		xListRan := getDifferentBenchParameters(4)
		x4 := (*[4]*big.Int)(xListRan)
		b.StartTimer()
		FourfoldExpPrecomputed(g, n, *x4, table)
		b.StopTimer()
	}
}

func BenchmarkDefaultExp(b *testing.B) {
	g, n, xList := getBenchParameters(1)
	result := new(big.Int)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result.Exp(g, xList[0], n)
	}
}

func BenchmarkExpParallel1(b *testing.B) {
	g, n, xList := getBenchParameters(1)
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, xList[0], n, table, 1, 0)
	}
}

func BenchmarkExpParallel2(b *testing.B) {
	g, n, xList := getBenchParameters(1)
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, xList[0], n, table, 2, 0)
	}
}

func BenchmarkExpParallel4(b *testing.B) {
	g, n, xList := getBenchParameters(1)
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, xList[0], n, table, 4, 0)
	}
}

func BenchmarkExpParallel8(b *testing.B) {
	g, n, xList := getBenchParameters(1)
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, xList[0], n, table, 8, 0)
	}
}

func BenchmarkExpParallel16(b *testing.B) {
	g, n, xList := getBenchParameters(1)
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, xList[0], n, table, 16, 0)
	}
}
