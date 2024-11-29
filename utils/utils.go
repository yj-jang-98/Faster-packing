package utils

import (
	"fmt"
	"math"
	"strconv"
)

// Rounding of vector
func RoundVec(a []float64) []float64 {
	m := len(a)
	b := make([]float64, m)
	for i := 0; i < m; i++ {
		b[i] = math.Round(a[i])
	}
	return b
}

// Convert a matrix to a string
func MatToString(A [][]float64) [][]string {
	Astr := make([][]string, len(A))
	for i := 0; i < len(A); i++ {
		Astr[i] = make([]string, len(A[0]))
		for j := range A[0] {
			Astr[i][j] = strconv.FormatFloat(A[i][j], 'f', -1, 64)
		}
	}
	return Astr
}

func VecToString(A []float64) [][]string {
	Astr := make([][]string, len(A))
	for i := 0; i < len(A); i++ {
		Astr[i] = make([]string, 1)
		Astr[i][0] = strconv.FormatFloat(A[i], 'f', -1, 64)
	}
	return Astr
}

// A must have 1 column
func MatToVec(A [][]float64) []float64 {
	b := make([]float64, len(A))
	if len(A[0]) > 1 {
		panic(fmt.Errorf("input is not a column vector"))
	}

	for i := 0; i < len(A); i++ {
		b[i] = A[i][0]
	}

	return b
}

// Result: [A, B] of size n x (T+1)
func AppendVecToMat(A [][]float64, B [][]float64) [][]float64 {
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

// Duplicate + zero padding
// length h vector n times
// Length of a must be less than or equal to h
func VecDuplicate(a []float64, n int, h int) []float64 {
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

// mod q of float vector
// Each component of b belongs to [0, q).
func ModVecFloat(a []float64, q uint64) []uint64 {
	qf := float64(q)
	b := make([]uint64, len(a))
	for i := 0; i < len(a); i++ {
		b[i] = uint64(a[i] - math.Floor(a[i]/qf)*qf)
	}

	return b
}

// Each component of b belongs to (-T/2, T/2).
func SignFloat(a float64, T uint64) float64 {
	q := float64(T)
	b := a - math.Floor((a+q/2.0)/q)*q
	return b
}

func Average(a []float64) float64 {
	avg := float64(0)
	for i := 0; i < len(a); i++ {
		avg += a[i]
	}
	avg = avg / float64(len(a))
	return avg
}
