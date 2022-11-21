package multiexp

import (
	"crypto/rand"
	"fmt"
	"math/big"
	"testing"
)

func TestDoubleExp(t *testing.T) {
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
	N, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	fmt.Println("g = ", g.String())
	fmt.Println("x1 = ", x1.String())
	fmt.Println("x2 = ", x2.String())
	fmt.Println("N = ", N.String())

	result := DoubleExp(g, x1, x2, N)
	var result2 big.Int
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
}
