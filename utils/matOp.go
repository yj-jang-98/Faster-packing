package utils

import (
	"fmt"
	"math"
	"strconv"
)

// Add a matrix to a matrix
// Input
// - M1 : m x n float matrix
// - M2 : m x n float matrix
// Output
// - MOut: m x n float matrix
func MatAdd(M1 [][]float64, M2 [][]float64) [][]float64 {
	row := len(M1)
	col := len(M1[0])
	MOut := make([][]float64, row)
	for r := 0; r < row; r++ {
		MOut[r] = make([]float64, col)
		for c := 0; c < col; c++ {
			MOut[r][c] = M1[r][c] + M2[r][c]
		}
	}
	return MOut
}

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

// Multiply a matrix to a matrix
// Input
// - M1 : m x n float matrix
// - M2 : n x l float matrix
// Output
// - MOut: m x l float matrix
func MatMatMult(M1 [][]float64, M2 [][]float64) [][]float64 {
	row1 := len(M1)
	col1 := len(M1[0])
	row2 := len(M2)
	col2 := len(M2[0])

	if col1 != row2 {
		panic(fmt.Errorf("matrix dimension do not match"))
	}

	MOut := make([][]float64, row1)
	for r := 0; r < row1; r++ {
		MOut[r] = make([]float64, col2)
		for c := 0; c < col2; c++ {
			tmp := 0.0
			for k := 0; k < col1; k++ {
				tmp = tmp + M1[r][k]*M2[k][c]
			}
			MOut[r][c] = tmp
		}
	}
	return MOut
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

// Append a vector to a matrix (column-wise)
// Input
// - M      : n x T float matrix
// - v      : n x 1 float matrix
// Output
// - MOUt = [M, v] : n x (T+1) float matrix
func AppendVecToMat(M [][]float64, v [][]float64) [][]float64 {
	vec := MatToVec(v)
	M = append(M, vec)
	return M
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
