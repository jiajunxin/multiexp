package multiexp

import (
	"crypto/rand"
	"fmt"
	"io"
	"math/big"
	"testing"

	. "math/big"
)

func getValidModulus(r io.Reader, max *Int) *Int {
	for true {
		N, err := rand.Int(r, max)
		if err != nil {
			panic(err)
		}
		if N.Bits()[0]&1 == 1 {
			return N
		}
	}
	return nil
}

func TestDoubleExp(t *testing.T) {
	var max Int
	max.SetInt64(1000000)

	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x1, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x2, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}

	N := getValidModulus(rand.Reader, &max)

	result := DoubleExp(g, x1, x2, N)
	var result2 Int
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
}

func TestFourfoldExp(t *testing.T) {
	var max big.Int
	max.SetInt64(1000000)

	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x1, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x2, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x3, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x4, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	N := getValidModulus(rand.Reader, &max)

	result := FourfoldExp(g, N, []*Int{x1, x2, x3, x4})
	var result2 Int
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x3, N)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x4, N)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	g.SetInt64(1000000)
	x1.SetInt64(2000000)
	x2.SetInt64(3000000)
	x3.SetInt64(4000000)
	x4.SetInt64(5000000)
	N.SetInt64(2000001)
	result = FourfoldExp(g, N, []*Int{x1, x2, x3, x4})
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x3, N)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x4, N)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
}

func TestFourfoldExpParallel(t *testing.T) {
	var max big.Int
	// We need max to be larger to make the precompute actually work.
	max.SetInt64(1000000000) //2^30
	max.Mul(&max, &max)      //2^60
	max.Mul(&max, &max)      //2^120

	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x1, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x2, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x3, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	x4, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	N := getValidModulus(rand.Reader, &max)
	maxLen := (max.BitLen() / GetWidth()) + 1
	fmt.Println("BitLen = ", max.BitLen())
	fmt.Println("maxLen = ", maxLen)
	table := PreCompute(g, N, maxLen)
	result := FourfoldExpWithPreComputeTableParallel(g, N, []*Int{x1, x2, x3, x4}, table)
	var result2 Int
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
	result2.Exp(g, x3, N)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
	result2.Exp(g, x4, N)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
	g.SetInt64(1000000)
	x1.SetInt64(2000000)
	x2.SetInt64(3000000)
	x3.SetInt64(4000000)
	x4.SetInt64(5000000)
	N.SetInt64(2000001)
	table = PreCompute(g, N, maxLen)
	result = FourfoldExpWithPreComputeTableParallel(g, N, []*Int{x1, x2, x3, x4}, table)
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
	result2.Exp(g, x3, N)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
	result2.Exp(g, x4, N)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for FourfoldExpParallel")
	}
}
