package multiexp

import (
	"fmt"
)

// GetWidth returns the width of uint in this system
func GetWidth() int {
	return _W
}

func StatforInt(input nat) {
	fmt.Println("Nat len of input = ", len(input))
	var counter uint64
	for i := range input {
		counter += Bit1Counter(input[i])
	}
	fmt.Println("Input has bit '1' = ", counter)
}

func Bit1Counter(input Word) uint64 {
	var ret uint64
	var mask uint
	for i := 0; i < 32; i++ {
		mask = uint(1 << i)
		if (uint(input) & mask) == mask {
			ret++
		}
	}
	return ret
}
