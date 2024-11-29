package utils

import (
	"math"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// Subtract a vector from a vector
func VecSub(v1 []float64, v2 []float64) []float64 {
	// v1 : m x 1
	// v2 : m x 1

	m := len(v1)

	vReturn := make([]float64, m)
	for i := 0; i < m; i++ {
		vReturn[i] = v1[i] - v2[i]
	}
	return vReturn
}

func VecAdd(v1 []float64, v2 []float64) []float64 {
	// v1 : m x 1
	// v2 : m x 1

	m := len(v1)

	vReturn := make([]float64, m)

	for i := 0; i < m; i++ {
		vReturn[i] = v1[i] + v2[i]
	}
	return vReturn
}

// Multiply a scalar to a vector
func ScalarVecMult(s float64, v []float64) []float64 {
	// v : m x 1

	m := len(v)

	vReturn := make([]float64, m)
	for i := 0; i < m; i++ {
		vReturn[i] = s * v[i]
	}
	return vReturn
}

// Inner sum of vector([]uint64) elements (mod q)
func VecSumUint(a []uint64, T uint64, bredparams [2]uint64) uint64 {
	sum := uint64(0)
	for i := 0; i < len(a); i++ {
		sum += a[i]
	}
	sum = ring.BRedAdd(sum, T, bredparams)
	return sum
}

// 2-norm of a vector
func Vec2Norm(v []float64) float64 {
	tmp := 0.0
	for i := range v {
		tmp = tmp + v[i]*v[i]
	}
	return math.Sqrt(tmp)
}
