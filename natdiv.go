// Copyright 2009 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package multiexp

import (
	"math/bits"
)

// div returns q, r such that q = ⌊u/v⌋ and r = u%v = u - q·v.
// It uses z and z2 as the storage for q and r.
func (z nat) div(z2, u, v nat) (q, r nat) {
	if len(v) == 0 {
		panic("division by zero")
	}

	if u.cmp(v) < 0 {
		q = z[:0]
		r = z2.set(u)
		return
	}

	if len(v) == 1 {
		// Short division: long optimized for a single-word divisor.
		// In that case, the 2-by-1 guess is all we need at each step.
		var r2 Word
		q, r2 = z.divW(u, v[0])
		r = z2.setWord(r2)
		return
	}

	q, r = z.divLarge(z2, u, v)
	return
}

// divW returns q, r such that q = ⌊x/y⌋ and r = x%y = x - q·y.
// It uses z as the storage for q.
// Note that y is a single digit (Word), not a big number.
func (z nat) divW(x nat, y Word) (q nat, r Word) {
	m := len(x)
	switch {
	case y == 0:
		panic("division by zero")
	case y == 1:
		q = z.set(x) // result is x
		return
	case m == 0:
		q = z[:0] // result is 0
		return
	}
	// m > 0
	z = z.make(m)
	r = divWVW(z, 0, x, y)
	q = z.norm()
	return
}

// divWVW overwrites z with ⌊x/y⌋, returning the remainder r.
// The caller must ensure that len(z) = len(x).
func divWVW(z []Word, xn Word, x []Word, y Word) (r Word) {
	r = xn
	if len(x) == 1 {
		qq, rr := bits.Div(uint(r), uint(x[0]), uint(y))
		z[0] = Word(qq)
		return Word(rr)
	}
	rec := reciprocalWord(y)
	for i := len(z) - 1; i >= 0; i-- {
		z[i], r = divWW(r, x[i], y, rec)
	}
	return r
}

// div returns q, r such that q = ⌊uIn/vIn⌋ and r = uIn%vIn = uIn - q·vIn.
// It uses z and u as the storage for q and r.
// The caller must ensure that len(vIn) ≥ 2 (use divW otherwise)
// and that len(uIn) ≥ len(vIn) (the answer is 0, uIn otherwise).
func (z nat) divLarge(u, uIn, vIn nat) (q, r nat) {
	n := len(vIn)
	m := len(uIn) - n

	// Scale the inputs so vIn's top bit is 1 (see “Scaling Inputs” above).
	// vIn is treated as a read-only input (it may be in use by another
	// goroutine), so we must make a copy.
	// uIn is copied to u.
	shift := nlz(vIn[n-1])
	vp := getNat(n)
	v := *vp
	shlVU(v, vIn, shift)
	u = u.make(len(uIn) + 1)
	u[len(uIn)] = shlVU(u[0:len(uIn)], uIn, shift)

	// The caller should not pass aliased z and u, since those are
	// the two different outputs, but correct just in case.
	if alias(z, u) {
		z = nil
	}
	q = z.make(m + 1)

	// Use basic or recursive long division depending on size.
	if n < divRecursiveThreshold {
		q.divBasic(u, v)
	} else {
		q.divRecursive(u, v)
	}
	putNat(vp)

	q = q.norm()

	// Undo scaling of remainder.
	shrVU(u, u, shift)
	r = u.norm()

	return q, r
}

// divBasic implements long division as described above.
// It overwrites q with ⌊u/v⌋ and overwrites u with the remainder r.
// q must be large enough to hold ⌊u/v⌋.
func (q nat) divBasic(u, v nat) {
	n := len(v)
	m := len(u) - n

	qhatvp := getNat(n + 1)
	qhatv := *qhatvp

	// Set up for divWW below, precomputing reciprocal argument.
	vn1 := v[n-1]
	rec := reciprocalWord(vn1)

	// Compute each digit of quotient.
	for j := m; j >= 0; j-- {
		// Compute the 2-by-1 guess q̂.
		// The first iteration must invent a leading 0 for u.
		qhat := Word(_M)
		var ujn Word
		if j+n < len(u) {
			ujn = u[j+n]
		}

		// ujn ≤ vn1, or else q̂ would be more than one digit.
		// For ujn == vn1, we set q̂ to the max digit M above.
		// Otherwise, we compute the 2-by-1 guess.
		if ujn != vn1 {
			var rhat Word
			qhat, rhat = divWW(ujn, u[j+n-1], vn1, rec)

			// Refine q̂ to a 3-by-2 guess. See “Refining Guesses” above.
			vn2 := v[n-2]
			x1, x2 := mulWW(qhat, vn2)
			ujn2 := u[j+n-2]
			for greaterThan(x1, x2, rhat, ujn2) { // x1x2 > r̂ u[j+n-2]
				qhat--
				prevRhat := rhat
				rhat += vn1
				// If r̂  overflows, then
				// r̂ u[j+n-2]v[n-1] is now definitely > x1 x2.
				if rhat < prevRhat {
					break
				}
				// TODO(rsc): No need for a full mulWW.
				// x2 += vn2; if x2 overflows, x1++
				x1, x2 = mulWW(qhat, vn2)
			}
		}

		// Compute q̂·v.
		qhatv[n] = mulAddVWW(qhatv[0:n], v, qhat, 0)
		qhl := len(qhatv)
		if j+qhl > len(u) && qhatv[n] == 0 {
			qhl--
		}

		// Subtract q̂·v from the current section of u.
		// If it underflows, q̂·v > u, which we fix up
		// by decrementing q̂ and adding v back.
		c := subVV(u[j:j+qhl], u[j:], qhatv)
		if c != 0 {
			c := addVV(u[j:j+n], u[j:], v)
			// If n == qhl, the carry from subVV and the carry from addVV
			// cancel out and don't affect u[j+n].
			if n < qhl {
				u[j+n] += c
			}
			qhat--
		}

		// Save quotient digit.
		// Caller may know the top digit is zero and not leave room for it.
		if j == m && m == len(q) && qhat == 0 {
			continue
		}
		q[j] = qhat
	}

	putNat(qhatvp)
}

// greaterThan reports whether the two digit numbers x1 x2 > y1 y2.
// TODO(rsc): In contradiction to most of this file, x1 is the high
// digit and x2 is the low digit. This should be fixed.
func greaterThan(x1, x2, y1, y2 Word) bool {
	return x1 > y1 || x1 == y1 && x2 > y2
}

// divRecursiveThreshold is the number of divisor digits
// at which point divRecursive is faster than divBasic.
const divRecursiveThreshold = 100

// divRecursive implements recursive division as described above.
// It overwrites z with ⌊u/v⌋ and overwrites u with the remainder r.
// z must be large enough to hold ⌊u/v⌋.
// This function is just for allocating and freeing temporaries
// around divRecursiveStep, the real implementation.
func (z nat) divRecursive(u, v nat) {
	// Recursion depth is (much) less than 2 log₂(len(v)).
	// Allocate a slice of temporaries to be reused across recursion,
	// plus one extra temporary not live across the recursion.
	recDepth := 2 * bits.Len(uint(len(v)))
	tmp := getNat(3 * len(v))
	temps := make([]*nat, recDepth)

	z.clear()
	z.divRecursiveStep(u, v, 0, tmp, temps)

	// Free temporaries.
	for _, n := range temps {
		if n != nil {
			putNat(n)
		}
	}
	putNat(tmp)
}

// divRecursiveStep is the actual implementation of recursive division.
// It adds ⌊u/v⌋ to z and overwrites u with the remainder r.
// z must be large enough to hold ⌊u/v⌋.
// It uses temps[depth] (allocating if needed) as a temporary live across
// the recursive call. It also uses tmp, but not live across the recursion.
func (z nat) divRecursiveStep(u, v nat, depth int, tmp *nat, temps []*nat) {
	// u is a subsection of the original and may have leading zeros.
	// TODO(rsc): The v = v.norm() is useless and should be removed.
	// We know (and require) that v's top digit is ≥ B/2.
	u = u.norm()
	v = v.norm()
	if len(u) == 0 {
		z.clear()
		return
	}

	// Fall back to basic division if the problem is now small enough.
	n := len(v)
	if n < divRecursiveThreshold {
		z.divBasic(u, v)
		return
	}

	// Nothing to do if u is shorter than v (implies u < v).
	m := len(u) - n
	if m < 0 {
		return
	}

	// We consider B digits in a row as a single wide digit.
	// (See “Recursive Division” above.)
	//
	// TODO(rsc): rename B to Wide, to avoid confusion with _B,
	// which is something entirely different.
	// TODO(rsc): Look into whether using ⌈n/2⌉ is better than ⌊n/2⌋.
	B := n / 2

	// Allocate a nat for qhat below.
	if temps[depth] == nil {
		temps[depth] = getNat(n) // TODO(rsc): Can be just B+1.
	} else {
		*temps[depth] = temps[depth].make(B + 1)
	}

	// Compute each wide digit of the quotient.
	//
	// TODO(rsc): Change the loop to be
	//	for j := (m+B-1)/B*B; j > 0; j -= B {
	// which will make the final step a regular step, letting us
	// delete what amounts to an extra copy of the loop body below.
	j := m
	for j > B {
		// Divide u[j-B:j+n] (3 wide digits) by v (2 wide digits).
		// First make the 2-by-1-wide-digit guess using a recursive call.
		// Then extend the guess to the full 3-by-2 (see “Refining Guesses”).
		//
		// For the 2-by-1-wide-digit guess, instead of doing 2B-by-B-digit,
		// we use a (2B+1)-by-(B+1) digit, which handles the possibility that
		// the result has an extra leading 1 digit as well as guaranteeing
		// that the computed q̂ will be off by at most 1 instead of 2.

		// s is the number of digits to drop from the 3B- and 2B-digit chunks.
		// We drop B-1 to be left with 2B+1 and B+1.
		s := B - 1

		// uu is the up-to-3B-digit section of u we are working on.
		uu := u[j-B:]

		// Compute the 2-by-1 guess q̂, leaving r̂ in uu[s:B+n].
		qhat := *temps[depth]
		qhat.clear()
		qhat.divRecursiveStep(uu[s:B+n], v[s:], depth+1, tmp, temps)
		qhat = qhat.norm()

		// Extend to a 3-by-2 quotient and remainder.
		// Because divRecursiveStep overwrote the top part of uu with
		// the remainder r̂, the full uu already contains the equivalent
		// of r̂·B + uₙ₋₂ from the “Refining Guesses” discussion.
		// Subtracting q̂·vₙ₋₂ from it will compute the full-length remainder.
		// If that subtraction underflows, q̂·v > u, which we fix up
		// by decrementing q̂ and adding v back, same as in long division.

		// TODO(rsc): Instead of subtract and fix-up, this code is computing
		// q̂·vₙ₋₂ and decrementing q̂ until that product is ≤ u.
		// But we can do the subtraction directly, as in the comment above
		// and in long division, because we know that q̂ is wrong by at most one.
		qhatv := tmp.make(3 * n)
		qhatv.clear()
		qhatv = qhatv.mul(qhat, v[:s])
		for i := 0; i < 2; i++ {
			e := qhatv.cmp(uu.norm())
			if e <= 0 {
				break
			}
			subVW(qhat, qhat, 1)
			c := subVV(qhatv[:s], qhatv[:s], v[:s])
			if len(qhatv) > s {
				subVW(qhatv[s:], qhatv[s:], c)
			}
			addAt(uu[s:], v[s:], 0)
		}
		if qhatv.cmp(uu.norm()) > 0 {
			panic("impossible")
		}
		c := subVV(uu[:len(qhatv)], uu[:len(qhatv)], qhatv)
		if c > 0 {
			subVW(uu[len(qhatv):], uu[len(qhatv):], c)
		}
		addAt(z, qhat, j-B)
		j -= B
	}

	// TODO(rsc): Rewrite loop as described above and delete all this code.

	// Now u < (v<<B), compute lower bits in the same way.
	// Choose shift = B-1 again.
	s := B - 1
	qhat := *temps[depth]
	qhat.clear()
	qhat.divRecursiveStep(u[s:].norm(), v[s:], depth+1, tmp, temps)
	qhat = qhat.norm()
	qhatv := tmp.make(3 * n)
	qhatv.clear()
	qhatv = qhatv.mul(qhat, v[:s])
	// Set the correct remainder as before.
	for i := 0; i < 2; i++ {
		if e := qhatv.cmp(u.norm()); e > 0 {
			subVW(qhat, qhat, 1)
			c := subVV(qhatv[:s], qhatv[:s], v[:s])
			if len(qhatv) > s {
				subVW(qhatv[s:], qhatv[s:], c)
			}
			addAt(u[s:], v[s:], 0)
		}
	}
	if qhatv.cmp(u.norm()) > 0 {
		panic("impossible")
	}
	c := subVV(u[0:len(qhatv)], u[0:len(qhatv)], qhatv)
	if c > 0 {
		c = subVW(u[len(qhatv):], u[len(qhatv):], c)
	}
	if c > 0 {
		panic("impossible")
	}

	// Done!
	addAt(z, qhat.norm(), 0)
}
