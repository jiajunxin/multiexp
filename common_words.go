package multiexp

// gcw inputs two positive integer a and b, calculates the most common words
// i.e. a = 11011111, b = 11100000, most common word(s) = 11000000
func gcw(a, b nat) (nat, nat, nat) {
	aExtra := nat(nil).make(len(a))
	bExtra := nat(nil).make(len(b))
	var minWordLen int
	if len(a) > len(b) {
		minWordLen = len(b)
		for i := minWordLen; i < len(a); i++ {
			aExtra[i] = a[i]
		}
	} else {
		minWordLen = len(a)
		for i := minWordLen; i < len(b); i++ {
			bExtra[i] = b[i]
		}
	}

	commonWords := nat(nil).make(minWordLen)
	for i := 0; i < minWordLen; i++ {
		commonWords[i] = a[i] & b[i]
		aExtra[i] = a[i] - commonWords[i]
		bExtra[i] = b[i] - commonWords[i]
	}

	return aExtra, bExtra, commonWords
}

// fourfoldGCW inputs four positive integer a, b, c, d and calculates the greatest common words
// the last element in output is the common word slice
func fourfoldGCW(input [4]nat) [5]nat {
	maxWordLen := 0
	minWordLen := len(input[0])
	for i := 0; i < 4; i++ {
		if maxWordLen < len(input[i]) {
			maxWordLen = len(input[i])
		}
		if minWordLen > len(input[i]) {
			minWordLen = len(input[i])
		}
	}

	var outputs [5]nat
	for i := 0; i < 4; i++ {
		outputs[i] = outputs[i].make(len(input[i]))
	}
	outputs[4] = outputs[4].make(minWordLen)
	for i := 0; i < minWordLen; i++ {
		outputs[4][i] = input[0][i] & input[1][i] & input[2][i] & input[3][i]
		outputs[0][i] = input[0][i] - outputs[4][i]
		outputs[1][i] = input[1][i] - outputs[4][i]
		outputs[2][i] = input[2][i] - outputs[4][i]
		outputs[3][i] = input[3][i] - outputs[4][i]
	}
	for i := 0; i < 4; i++ {
		if len(outputs[i]) > minWordLen {
			for j := minWordLen; j < len(outputs[i]); j++ {
				outputs[i][j] = input[i][j]
			}
		}
	}

	return outputs
}

// threefoldGcb inputs three positive integer a, b, c and calculates the greatest common words
// the last element in output is the common word slice
func threefoldGCW(input [3]nat) nat {
	maxWordLen := 0
	minWordLen := len(input[0])
	for i := 0; i < 3; i++ {
		if maxWordLen < len(input[i]) {
			maxWordLen = len(input[i])
		}
		if minWordLen > len(input[i]) {
			minWordLen = len(input[i])
		}
	}

	output := nat(nil).make(minWordLen)
	for i := 0; i < minWordLen; i++ {
		output[i] = input[0][i] & input[1][i] & input[2][i]
		input[0][i] = input[0][i] - output[i]
		input[1][i] = input[1][i] - output[i]
		input[2][i] = input[2][i] - output[i]
	}
	return output
}
