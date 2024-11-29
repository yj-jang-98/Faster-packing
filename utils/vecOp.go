package utils

import (
	"fmt"
	"math"
	"strconv"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// Add a vector to a vector
// Input
// - v1 : m x 1 float vector
// - v2 : m x 1 float vector
// Output
// - vOut: m x 1 float vector
func VecAdd(v1 []float64, v2 []float64) []float64 {
	row := len(v1)
	vOut := make([]float64, row)
	for r := 0; r < row; r++ {
		vOut[r] = v1[r] + v2[r]
	}
	return vOut
}

// Subtract a vector to a vector
// Input
// - v1 : m x 1 float vector
// - v2 : m x 1 float vector
// Output
// - vOut: m x 1 float vector
func VecSub(v1 []float64, v2 []float64) []float64 {
	row := len(v1)
	vOut := make([]float64, row)
	for r := 0; r < row; r++ {
		vOut[r] = v1[r] - v2[r]
	}
	return vOut
}

// Multiply a scalar to a vector
// Input
// - s : scalar float
// - v : m x 1 float vector
// Output
// - vOut: m x 1 float vector
func ScalarVecMult(s float64, v []float64) []float64 {
	row := len(v)

	vOut := make([]float64, row)
	for r := 0; r < row; r++ {
		vOut[r] = s * v[r]
	}
	return vOut
}

// Inner sum of vector elements (mod q)
// Input
// - v : m x 1 uint64 vector
// Output
// - sum: scalar uint64
func VecSumUint(v []uint64, T uint64, bredparams [2]uint64) uint64 {
	row := len(v)
	sum := uint64(0)
	for r := 0; r < row; r++ {
		sum += v[r]
	}
	sum = ring.BRedAdd(sum, T, bredparams)
	return sum
}

// Mod q of float vector
// Each component of vOut belongs to [0, q).
// Input
// - v : m x 1 float vector
// Output
// - vOut: m x 1 uint64 vector
func ModVecFloat(v []float64, q uint64) []uint64 {
	qf := float64(q)
	row := len(v)
	vOut := make([]uint64, row)
	for r := 0; r < row; r++ {
		vOut[r] = uint64(v[r] - math.Floor(v[r]/qf)*qf)
	}
	return vOut
}

// ?????????
// Duplicate + zero padding
// length h vector n times
// Length of a must be less than or equal to h
// Input
// - v : n x 1 float vector
// Output
// - v: n x 1 float vector ???
func VecDuplicate(v []float64, n int, h int) []float64 {
	row := len(v)
	if row > h {
		panic(fmt.Errorf("length of a is greater than h"))
	}
	if row < h {
		for k := 0; k < h-row; k++ {
			v = append(v, 0)
		}
	}
	for i := 1; i < n; i++ {
		for j := 0; j < row; j++ {
			v = append(v, v[j])
		}
		if row < h {
			for k := 0; k < h-row; k++ {
				v = append(v, 0)
			}
		}
	}
	return v
}

// 2-norm of a vector
// Input
// - v : n x 1 float vector
// Output
// - s: scalar float
func Vec2Norm(v []float64) float64 {
	tmp := 0.0
	for i := range v {
		tmp = tmp + v[i]*v[i]
	}
	return math.Sqrt(tmp)
}

// Rounding of a vector
// Input
// - v : n x 1 float vector
// Output
// - v : n x 1 float vector
func RoundVec(v []float64) []float64 {
	row := len(v)
	vOut := make([]float64, row)
	for r := 0; r < row; r++ {
		vOut[r] = math.Round(v[r])
	}
	return vOut
}

// Convert a vector to a string
// Input
// - v : n x 1 float vector
// Output
// - strOut : n x 1 string matrix
func VecToString(v []float64) [][]string {
	row := len(v)
	strOut := make([][]string, row)
	for r := 0; r < row; r++ {
		strOut[r] = make([]string, 1)
		strOut[r][0] = strconv.FormatFloat(v[r], 'f', -1, 64)
	}
	return strOut
}
