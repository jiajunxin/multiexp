package multiexp

import (
	"crypto/rand"
	"sync"
	"testing"

	. "math/big"
)

const numBits = 20000

var (
	randLimit     *Int
	onceRandLimit sync.Once
)

func getRandLimit() *Int {
	onceRandLimit.Do(func() {
		randLimit = new(Int).SetInt64(1)
		randLimit.Lsh(randLimit, numBits)
	})
	return randLimit
}

func BenchmarkOriginalDoubleExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getRandLimit())
	b.ResetTimer()
	var result Int
	for i := 0; i < b.N; i++ {
		result.Exp(g, x[i], N)
		result.Exp(g, x[i+1], N)
	}
}

func BenchmarkDoubleExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getRandLimit())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		DoubleExp(g, x[i], x[i+1], N)
	}
}

func BenchmarkOriginalFourfoldExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getRandLimit())
	b.ResetTimer()
	var result Int
	for i := 0; i < b.N; i++ {
		result.Exp(g, x[i], N)
		result.Exp(g, x[i+1], N)
		result.Exp(g, x[i+2], N)
		result.Exp(g, x[i+3], N)
	}
}

func BenchmarkFourfoldExp(b *testing.B) {
	g, err := rand.Int(rand.Reader, getRandLimit())
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, getRandLimit())
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, getRandLimit())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FourfoldExp(g, N, x[i:i+4])
	}
}
