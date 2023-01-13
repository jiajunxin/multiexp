package multiexp

import (
	"crypto/rand"
	"math/big"
	"sync"
	"testing"
)

const numTestBits = 200000
const numTestGroupBits = 2000

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
	g, n, xList := getBenchParameters(2)
	b.ResetTimer()
	var result big.Int
	for i := 0; i < b.N; i++ {
		result.Exp(g, xList[0], n)
		result.Exp(g, xList[1], n)
	}
}

func BenchmarkDoubleExp(b *testing.B) {
	g, n, xList := getBenchParameters(2)
	x2 := (*[2]*big.Int)(xList)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		DoubleExp(g, *x2, n)
	}
}

func BenchmarkOriginalFourfoldExp(b *testing.B) {
	g, n, xList := getBenchParameters(4)
	b.ResetTimer()
	var result big.Int
	for i := 0; i < b.N; i++ {
		result.Exp(g, xList[0], n)
		result.Exp(g, xList[1], n)
		result.Exp(g, xList[2], n)
		result.Exp(g, xList[3], n)
	}
}

func BenchmarkFourfoldExp(b *testing.B) {
	g, n, xList := getBenchParameters(4)
	x4 := (*[4]*big.Int)(xList)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FourfoldExp(g, n, *x4)
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
