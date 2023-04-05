# Fast big-integer multi-output exponentiation

We observe that to calculate $g^x \mod N$ and $g^y \mod N$, we can combine them using the same squaring of $g$.
Furthermore, we can extract the common ``1" bits of $x$ and $y$ as $z$ and have $x^{\prime}= x-z$, $y^{\prime}= y-z$ to save common multiplications during the process.

Based on the same idea, we can further combine the process of $g^{x_1} \mod N$, $g^{x_2} \mod N$, $g^{x_3} \mod N$, $g^{x_4} \mod N$ by extracting all combinations of their common bits.
We first extract common bits for $\{x_1, x_2, x_3, x_4\}$ and update the original values; then we extract common bits of $\{x_1, x_2, x_3\}$, $\{x_1, x_2, x_4\}$, $\{x_1, x_3, x_4\}$ and $\{x_2, x_3, x_4\}$ respectively and update the original values;
finally, we extract $\binom{4}{2}$ combinations for each pair of values and update.
In this way, we get four updated values $\{x_1^{\prime}, x_2^{\prime}, x_3^{\prime}, x_4^{\prime}\}$ that are sparse in terms of "1" bit, and their original common "1" bits that can be calculated at once. 

In this library, we optimize multi-output exponentiations till 4 numbers. For larger values, combining more exponents becomes harder
For example, if we want to find all the common bits of 4 integers, the sum of all the combinations is $\binom{4}{2}+\binom{4}{3}+\binom{4}{4}=O(2^4)$.
e.g., for the common bits of $k$ integers, the sum of all the combinations is $O(2^k)$, which makes implementation much more complicated.

The most advanced method to calculate the modular exponent is Montgomery modular multiplication with a binary expression of exponents.
We modify the Montgomery modular multiplication based on Golang's Big library.
Besides, we also implement precomputation tables for even better speed-up. 

## Benchmark
We generate a random base $g$ and modulus $N$ with 2048 bits. To support Montgomery modular multiplication, we restrict $N$ to be odd, which is sufficient for RSA groups and Montgomery modular multiplication.
We use the same $g$ and $N$ for each run and random numbers as exponents, whose range varies from 2K-2M bits, and we report the average across 10 runs. 
The following table showcases our results in terms of function completion time and the percentage performance gain between our ''optimized'' execution over the ''naive'' one.

|                            | Time(ms) for # of bits=2,000 | Time(ms) for # of bits=20,000 | Time(ms) for # of bits=200,000 | Time(ms) for # of bits=2,000,000 |
| -------------------------- | ---------------------------- | ----------------------------- | ------------------------------ | -------------------------------- |
| 2×NaiveExponent            | 11.1                         | 94.9                          | 911.9                          | 9443.7                           |
| DoubleExponent             | 7.6                          | 66.9                          | 646.4                          | 6640.3                           |
| 4×NaiveExponent            | 22.3                         | 187.5                         | 1827.6                         | 18830.5                          |
| FourfoldExponent           | 9.9                          | 77.1                          | 733.9                          | 7482.0                           |
| precomputeFourfoldExponent | 5.5                          | 36.9                          | 359.7                          | 3596.5                           |

Looking at the specific functions, first, we calculate $g^{x_1} \mod N$ and $g^{x_2} \mod N$ ''naively'' by calculating each exponent separately. 
We refer to this as the 2×NaiveExponent function and we compare it against our optimized version, DoubleExponent.
Similarly, for random values $x_1,x_2,x_3,x_4$, we calculate $g^{x_1} \mod N$, $g^{x_2} \mod N$, $g^{x_3} \mod N$ and $g^{x_4} \mod N$ ``naively'' and we refer to this function as 4×NaiveExponent. 
We compare it against our optimized FourfoldExponent function.
The result of 2,000 bits exponent is the average of 10,000 runs, the result of 20,000 bits exponent is the average of 1,000 runs,
the result of 200,000 bits exponent is the average of 100 runs and the result of 2,000,000 bits exponent is the average of 10 runs.

We observe that with our optimizations DoubleExponent is around 30% faster and FourfoldExponent around 60% faster than their original counterparts.
Additionally, we report the performance of FourfoldExponent when combined with a precomputation table that includes a precomputation of every single bit; we refer to this as the precomputeFourfoldExponent function.
Precomputation tables with more elaborate settings (e.g., including four combinations of every two precomputed bits) lead to faster calculation and larger table sizes.
We list the precomputation table size with respect to the maximum exponent bits supported.

| Max exponent bits | $2^{20}$ | $2^{22}$ | $2^{24}$ | $2^{26}$ | $2^{28}$ |
| ----------------- | -------- | -------- | -------- | -------- | -------- |
| Table Size (MB)   | 8        | 32       | 128      | 512      | 2048     |
