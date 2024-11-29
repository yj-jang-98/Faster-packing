package utils

import (
	"fmt"
)

// Add a matrix to a matrix
func MatAdd(M1 [][]float64, M2 [][]float64) [][]float64 {
	// M1 : m x n
	// M2 : m x n

	m := len(M1)
	n := len(M1[0])

	MReturn := make([][]float64, m)
	for i := 0; i < m; i++ {
		MReturn[i] = make([]float64, n)
	}

	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			MReturn[i][j] = M1[i][j] + M2[i][j]
		}
	}
	return MReturn
}

// Multiply a scalar to a matrix
func ScalarMatMult(s float64, M [][]float64) [][]float64 {
	// M : m x n

	m := len(M)
	n := len(M[0])

	MReturn := make([][]float64, m)
	for i := 0; i < m; i++ {
		MReturn[i] = make([]float64, n)
		for j := range M[i] {
			MReturn[i][j] = s * M[i][j]
		}
	}
	return MReturn
}

// Multiply a matrix to a matrix
func MatMatMult(M1 [][]float64, M2 [][]float64) [][]float64 {
	// A : m x n
	// B : n x l

	m := len(M1)
	n := len(M1[0])
	n1 := len(M2)
	l := len(M2[0])

	if n != n1 {
		panic(fmt.Errorf("matrix dimension do not match"))
	}

	MReturn := make([][]float64, m)
	for i := 0; i < m; i++ {
		MReturn[i] = make([]float64, l)
	}

	for i := 0; i < m; i++ {
		for j := 0; j < l; j++ {
			tmp := 0.0
			for k := 0; k < n; k++ {
				tmp = tmp + M1[i][k]*M2[k][j]
			}
			MReturn[i][j] = tmp
		}
	}
	return MReturn
}

func MatVecMult(M [][]float64, v []float64) []float64 {
	// A : m x n
	// B : n x l

	m := len(M)
	n := len(M[0])
	n1 := len(v)

	if n != n1 {
		panic(fmt.Errorf("Matrix dimension don't match"))
	}

	vReturn := make([]float64, m)

	for i := 0; i < m; i++ {
		tmp := 0.0
		for k := 0; k < n; k++ {
			tmp = tmp + M[i][k]*v[k]
		}
		vReturn[i] = tmp
	}
	return vReturn
}
