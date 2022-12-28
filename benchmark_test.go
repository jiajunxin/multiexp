package multiexp

import (
	"crypto/rand"
	"math/big"
	"sync"
	"testing"
)

const numBits = 20000

var (
	benchRandLimit      *big.Int
	onceBenchRandLimit  sync.Once
	g, x, mod           *big.Int
	onceBenchParameters sync.Once
	table               *PreTable
	onceBenchTable      sync.Once
)

func getBenchRandLimit() *big.Int {
	onceBenchRandLimit.Do(func() {
		benchRandLimit = new(big.Int).SetInt64(1)
		benchRandLimit.Lsh(benchRandLimit, numBits)
	})
	return benchRandLimit
}

func getBenchParameters() (*big.Int, *big.Int, *big.Int) {
	onceBenchParameters.Do(func() {
		g, x, mod = new(big.Int), new(big.Int), new(big.Int)
		g, _ = rand.Int(rand.Reader, getBenchRandLimit())
		x, _ = rand.Int(rand.Reader, getBenchRandLimit())
		mod = getValidModulus(rand.Reader, getBenchRandLimit())
	})
	return g, x, mod
}

func getBenchPrecomputeTable() *PreTable {
	onceBenchTable.Do(func() {
		g, _, N := getBenchParameters()
		randLmtLen := (getBenchRandLimit().BitLen() / GetWidth()) + 1
		table = NewPrecomputeTable(g, N, randLmtLen)
	})
	return table
}

func BenchmarkOriginalDoubleExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getBenchRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*big.Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getBenchRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getBenchRandLimit())
	b.ResetTimer()
	var result big.Int
	for i := 0; i < b.N; i++ {
		result.Exp(g, x[i], N)
		result.Exp(g, x[i+1], N)
	}
}

func BenchmarkDoubleExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getBenchRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*big.Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getBenchRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getBenchRandLimit())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		DoubleExp(g, x[i], x[i+1], N)
	}
}

func BenchmarkOriginalFourfoldExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getBenchRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*big.Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getBenchRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getBenchRandLimit())
	b.ResetTimer()
	var result big.Int
	for i := 0; i < b.N; i++ {
		result.Exp(g, x[i], N)
		result.Exp(g, x[i+1], N)
		result.Exp(g, x[i+2], N)
		result.Exp(g, x[i+3], N)
	}
}

func BenchmarkFourfoldExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getBenchRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*big.Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getBenchRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getBenchRandLimit())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FourfoldExp(g, N, x[i:i+4])
	}
}

func BenchmarkDefaultExp(b *testing.B) {
	g, x, N := getBenchParameters()
	result := new(big.Int)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result.Exp(g, x, N)
	}
}

func BenchmarkExpParallel1(b *testing.B) {
	g, x, N := getBenchParameters()
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, x, N, table, 1)
	}
}

func BenchmarkExpParallel4(b *testing.B) {
	g, x, N := getBenchParameters()
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, x, N, table, 4)
	}
}

func BenchmarkExpParallel8(b *testing.B) {
	g, x, N := getBenchParameters()
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, x, N, table, 8)
	}
}

func BenchmarkExpParallel16(b *testing.B) {
	g, x, N := getBenchParameters()
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, x, N, table, 16)
	}
}

func BenchmarkExpParallel32(b *testing.B) {
	g, x, N := getBenchParameters()
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, x, N, table, 32)
	}
}

func BenchmarkExpParallel64(b *testing.B) {
	g, x, N := getBenchParameters()
	table := getBenchPrecomputeTable()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ExpParallel(g, x, N, table, 64)
	}
}
