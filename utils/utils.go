package main

import (
	"fmt"
	"math"
	"strconv"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// Duplicate + zero padding
// length h vector n times
// Length of a must be less than or equal to h
func vecDuplicate(a []float64, n int, h int) []float64 {
	m := len(a)
	if m > h {
		panic(fmt.Errorf("length of a is greater than h"))
	}
	if m < h {
		for k := 0; k < h-m; k++ {
			a = append(a, 0)
		}
	}
	for i := 1; i < n; i++ {
		for j := 0; j < m; j++ {
			a = append(a, a[j])
		}
		if m < h {
			for k := 0; k < h-m; k++ {
				a = append(a, 0)
			}
		}
	}
	return a
}

// Rounding of vector
func roundVec(a []float64) []float64 {
	m := len(a)
	b := make([]float64, m)
	for i := 0; i < m; i++ {
		b[i] = math.Round(a[i])
	}
	return b
}

// Multiply a scalar to a vector
func scalarVecMult(s float64, M []float64) []float64 {
	C := make([]float64, len(M))
	for i := 0; i < len(M); i++ {
		C[i] = s * M[i]
	}
	return C
}

// mod q of float vector
// Each component of b belongs to [0, q).
func modVecFloat(a []float64, q uint64) []uint64 {
	qf := float64(q)
	b := make([]uint64, len(a))
	for i := 0; i < len(a); i++ {
		b[i] = uint64(a[i] - math.Floor(a[i]/qf)*qf)
	}

	return b
}

// Result: [A, B] of size n x (T+1)
func appendVecToMat(A [][]float64, B [][]float64) [][]float64 {
	// A : n x T
	// B : n x 1
	m := len(B)
	vec := make([]float64, m)
	for i := 0; i < m; i++ {
		vec[i] = B[i][0]
	}
	A = append(A, vec)
	return A
}

func matMult(A [][]float64, B [][]float64) [][]float64 {
	// A : m x n
	// B : n x l
	m := len(A)
	n := len(A[0])
	n1 := len(B)
	l := len(B[0])

	if n != n1 {
		panic(fmt.Errorf("matrix dimension do not match"))
	}

	C := make([][]float64, m)
	for i := 0; i < m; i++ {
		C[i] = make([]float64, l)
	}

	for i := 0; i < m; i++ {
		for j := 0; j < l; j++ {
			tmp := 0.0
			for k := 0; k < n; k++ {
				tmp = tmp + A[i][k]*B[k][j]
			}
			C[i][j] = tmp
		}
	}
	return C
}

func matAdd(A [][]float64, B [][]float64) [][]float64 {
	// A : m x n
	// B : m x n
	m := len(A)
	n := len(A[0])

	C := make([][]float64, m)
	for i := 0; i < m; i++ {
		C[i] = make([]float64, n)
	}

	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			C[i][j] = A[i][j] + B[i][j]
		}
	}
	return C
}

// A must have 1 column
func mat2vec(A [][]float64) []float64 {
	b := make([]float64, len(A))
	if len(A[0]) > 1 {
		panic(fmt.Errorf("input is not a column vector"))
	}

	for i := 0; i < len(A); i++ {
		b[i] = A[i][0]
	}

	return b
}

// Inner sum of vector([]uint64) elements (mod q)
func vecSumUint(a []uint64, T uint64, bredparams [2]uint64) uint64 {
	sum := uint64(0)
	for i := 0; i < len(a); i++ {
		sum += a[i]
	}
	sum = ring.BRedAdd(sum, T, bredparams)
	return sum
}

// Each component of b belongs to (-T/2, T/2).
func signFloat(a float64, T uint64) float64 {
	q := float64(T)
	b := a - math.Floor((a+q/2.0)/q)*q
	return b
}

func average(a []float64) float64 {
	avg := float64(0)
	for i := 0; i < len(a); i++ {
		avg += a[i]
	}
	avg = avg / float64(len(a))
	return avg
}

func subVec(v1 []float64, v2 []float64) []float64 {
	vReturn := make([]float64, len(v1))
	for i := range v1 {
		vReturn[i] = v1[i] - v2[i]
	}
	return vReturn
}

func vec2norm(v []float64) float64 {
	tmp := 0.0
	for i := range v {
		tmp = tmp + v[i]*v[i]
	}
	return math.Sqrt(tmp)
}

func mat2string(A [][]float64) [][]string {
	Astr := make([][]string, len(A))
	for i := 0; i < len(A); i++ {
		Astr[i] = make([]string, len(A[0]))
		for j := range A[0] {
			Astr[i][j] = strconv.FormatFloat(A[i][j], 'f', -1, 64)
		}
	}
	return Astr
}
