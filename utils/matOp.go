package utils

import (
	"fmt"
	"math"
	"strconv"
)

// Multiply a scalar to a matrix
// Input
// - s : scalar float
// - M : m x n float matrix
// Output
// - MOut: m x n float matrix
func ScalarMatMult(s float64, M [][]float64) [][]float64 {
	row := len(M)
	col := len(M[0])
	MOut := make([][]float64, row)
	for r := 0; r < row; r++ {
		MOut[r] = make([]float64, col)
		for c := 0; c < col; c++ {
			MOut[r][c] = s * M[r][c]
		}
	}
	return MOut
}

// Multiply a matrix to a vector
// Input
// - M : m x n float matrix
// - v : n x 1 float vector
// Output
// - vOut: m x l float vector
func MatVecMult(M [][]float64, v []float64) []float64 {
	row := len(M)
	col := len(M[0])
	colv := len(v)

	if col != colv {
		panic(fmt.Errorf("Matrix dimension don't match"))
	}

	vOut := make([]float64, row)
	for r := 0; r < row; r++ {
		tmp := 0.0
		for k := 0; k < col; k++ {
			tmp = tmp + M[r][k]*v[k]
		}
		vOut[r] = tmp
	}
	return vOut
}

// Mod q of float matrix
// Each component of MOut belongs to [0, q).
// Input
// - M : m x n float matrix
// Output
// - MOut: m x n uint64 matrix
func ModMatFloat(M [][]float64, q uint64) [][]uint64 {
	qf := float64(q)
	row := len(M)
	col := len(M[0])
	MOut := make([][]uint64, row)
	for r := 0; r < row; r++ {
		MOut[r] = make([]uint64, col)
		for c := 0; c < col; c++ {
			MOut[r][c] = uint64(M[r][c] - math.Floor(M[r][c]/qf)*qf)
		}
	}
	return MOut
}

// Convert a matrix to a vector
// The matrix must have 1 column
// Input
// - M : m x 1 float matrix
// Output
// - vOut : m x 1 float vector
func MatToVec(M [][]float64) []float64 {
	row := len(M)
	col := len(M[0])
	if col > 1 {
		panic(fmt.Errorf("input is not a column vector"))
	}

	vOut := make([]float64, row)
	for r := 0; r < row; r++ {
		vOut[r] = M[r][0]
	}

	return vOut
}

// Convert a matrix to a string
// Input
// - M : m x n float matrix
// Output
// - strOut : m x n string matrix
func MatToString(M [][]float64) [][]string {
	row := len(M)
	col := len(M[0])
	strOut := make([][]string, row)
	for r := 0; r < row; r++ {
		strOut[r] = make([]string, col)
		for c := 0; c < col; c++ {
			strOut[r][c] = strconv.FormatFloat(M[r][c], 'f', -1, 64)
		}
	}
	return strOut
}
