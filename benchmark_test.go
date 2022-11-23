package multiexp

import (
	"crypto/rand"
	. "math/big"
	"testing"
)

const numofBits = 20000

func BenchmarkOriginalDoubleExp(b *testing.B) {
	var max Int
	max.SetInt64(1)
	max.Lsh(&max, numofBits)
	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, &max)
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, &max)
	b.ResetTimer()
	var result Int
	for i := 0; i < b.N; i++ {
		result.Exp(g, x[i], N)
		result.Exp(g, x[i+1], N)
	}
}

func BenchmarkDoubleExp(b *testing.B) {
	var max Int
	max.SetInt64(1)
	max.Lsh(&max, numofBits)
	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, &max)
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, &max)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		DoubleExp(g, x[i], x[i+1], N)
	}
}

func BenchmarkOriginalFourfoldExp(b *testing.B) {
	var max Int
	max.SetInt64(1)
	max.Lsh(&max, numofBits)
	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, &max)
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, &max)
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
	var max Int
	max.SetInt64(1)
	max.Lsh(&max, numofBits)
	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		b.Errorf(err.Error())
	}
	x := make([]*Int, b.N+4)
	for i := range x {
		x[i], err = rand.Int(rand.Reader, &max)
		if err != nil {
			b.Errorf(err.Error())
		}
	}
	N := getValidModulus(rand.Reader, &max)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FourFoldExp(g, N, x[i:i+4])
	}
}
