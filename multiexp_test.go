package multiexp

import (
	"crypto/rand"
	"fmt"
	"io"
	"math/big"
	"reflect"
	"testing"
)

func getValidModulus(r io.Reader, max *big.Int) *big.Int {
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

	n := getValidModulus(rand.Reader, &max)

	result := DoubleExp(g, [2]*big.Int{x1, x2}, n)
	var result2 big.Int
	result2.Exp(g, x1, n)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
	result2.Exp(g, x2, n)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
}

func TestDoubleExpwithProd(t *testing.T) {
	setSize := 999
	var max, prod1, prod2 big.Int
	max.SetInt64(1000000)
	prod1.SetInt64(1)
	prod2.SetInt64(1)

	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	for i := 0; i < setSize; i++ {
		x1, err := rand.Int(rand.Reader, &max)
		if err != nil {
			t.Errorf(err.Error())
		}
		prod1.Mul(&prod1, x1)
	}

	for i := 0; i < setSize; i++ {
		x2, err := rand.Int(rand.Reader, &max)
		if err != nil {
			t.Errorf(err.Error())
		}
		prod2.Mul(&prod2, x2)
	}

	n := getValidModulus(rand.Reader, &max)

	result := DoubleExp(g, [2]*big.Int{&prod1, &prod2}, n)
	var two, temp1, temp2 big.Int
	two.SetInt64(2)
	temp1.Mod(&prod1, &two)
	temp2.Mod(&prod2, &two)
	fmt.Println("temp1 = ", temp1.String())
	fmt.Println("temp2 = ", temp2.String())
	var result2 big.Int
	result2.Exp(g, &prod1, n)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
	result2.Exp(g, &prod2, n)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
}

func TestDoubleExpwithProd2(t *testing.T) {
	setSize := 999
	var max, prod1, prod2 big.Int
	max.SetInt64(1000000)
	prod1.SetInt64(1)
	prod2.SetInt64(1)

	g, err := rand.Int(rand.Reader, &max)
	if err != nil {
		t.Errorf(err.Error())
	}
	for i := 0; i < setSize; i++ {
		x1 := getPrime256()
		prod1.Mul(&prod1, x1)
		x2 := getPrime256()
		prod2.Mul(&prod2, x2)
	}

	n := getValidModulus(rand.Reader, &max)

	result := DoubleExp(g, [2]*big.Int{&prod1, &prod2}, n)
	var result2 big.Int
	result2.Exp(g, &prod1, n)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
	result2.Exp(g, &prod2, n)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for DoubleExp")
	}
}

func getPrime256() *big.Int {
	flag := false
	for !flag {
		ranNum, err := rand.Prime(rand.Reader, 256)
		if err != nil {
			panic(err)
		}
		flag = ranNum.ProbablyPrime(80)
		if !flag {
			continue
		}
		return ranNum
	}
	return nil
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
	n := getValidModulus(rand.Reader, &max)

	result := FourfoldExp(g, n, [4]*big.Int{x1, x2, x3, x4})
	var result2 big.Int
	result2.Exp(g, x1, n)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x2, n)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x3, n)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x4, n)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	g.SetInt64(1000000)
	x1.SetInt64(2000000)
	x2.SetInt64(3000000)
	x3.SetInt64(4000000)
	x4.SetInt64(5000000)
	n.SetInt64(2000001)
	result = FourfoldExp(g, n, [4]*big.Int{x1, x2, x3, x4})
	result2.Exp(g, x1, n)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x2, n)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x3, n)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
	result2.Exp(g, x4, n)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for FourfoldExp")
	}
}

func TestFourfoldExpwithTable(t *testing.T) {
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
	maxLen := (max.BitLen() / _W) + 1
	// fmt.Println("BitLen = ", max.BitLen())
	// fmt.Println("maxLen = ", maxLen)
	table := NewPrecomputeTable(g, N, maxLen)
	result := FourfoldExpPrecomputed(g, N, [4]*big.Int{x1, x2, x3, x4}, table)
	var result2 big.Int
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
	}
	result2.Exp(g, x3, N)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
	}
	result2.Exp(g, x4, N)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
	}
	g.SetInt64(1000000)
	x1.SetInt64(2000000)
	x2.SetInt64(3000000)
	x3.SetInt64(4000000)
	x4.SetInt64(5000000)
	N.SetInt64(2000001)
	table = NewPrecomputeTable(g, N, maxLen)
	result = FourfoldExpPrecomputed(g, N, [4]*big.Int{x1, x2, x3, x4}, table)
	result2.Exp(g, x1, N)
	if result2.Cmp(result[0]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
	}
	result2.Exp(g, x2, N)
	if result2.Cmp(result[1]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
	}
	result2.Exp(g, x3, N)
	if result2.Cmp(result[2]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
	}
	result2.Exp(g, x4, N)
	if result2.Cmp(result[3]) != 0 {
		t.Errorf("Wrong result for TestFourfoldExpwithTable")
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
	maxLen := (max.BitLen() / _W) + 1
	// fmt.Println("BitLen = ", max.BitLen())
	// fmt.Println("maxLen = ", maxLen)
	table := NewPrecomputeTable(g, N, maxLen)
	result := FourfoldExpPrecomputedParallel(g, N, [4]*big.Int{x1, x2, x3, x4}, table)
	var result2 big.Int
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
	table = NewPrecomputeTable(g, N, maxLen)
	result = FourfoldExpPrecomputedParallel(g, N, [4]*big.Int{x1, x2, x3, x4}, table)
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

func TestExpParallel(t *testing.T) {
	//randLimit := new(big.Int)
	//randLimit.SetInt64(1000000000)
	//g, err := rand.Int(rand.Reader, randLimit)
	//if err != nil {
	//	t.Errorf(err.Error())
	//}
	//x, err := rand.Int(rand.Reader, randLimit)
	//if err != nil {
	//	t.Errorf(err.Error())
	//}
	//n := getValidModulus(rand.Reader, randLimit)
	g, n, xList := getBenchParameters(1)
	table := getBenchPrecomputeTable()
	type args struct {
		x             *big.Int
		y             *big.Int
		m             *big.Int
		preTable      *PreTable
		numRoutine    int
		wordChunkSize int
	}
	tests := []struct {
		name string
		args args
		want *big.Int
	}{
		{
			name: "TestExpParallel",
			args: args{
				x:          g,
				y:          xList[0],
				m:          n,
				preTable:   table,
				numRoutine: 4,
			},
			want: new(big.Int).Exp(g, xList[0], n),
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := ExpParallel(tt.args.x, tt.args.y, tt.args.m, tt.args.preTable, tt.args.numRoutine, tt.args.wordChunkSize); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("ExpParallel() = %v, want %v", got, tt.want)
			}
		})
	}
}
