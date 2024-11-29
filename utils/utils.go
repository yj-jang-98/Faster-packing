package utils

import (
	"math"
)

// Maps a float to interval (-T/2, T/2).
// Input
// - s : float
// Output
// - sOut: float
func SignFloat(s float64, T uint64) float64 {
	Tf := float64(T)
	sOut := s - math.Floor((s+Tf/2.0)/Tf)*Tf
	return sOut
}

// Takes average of a vector v
// Input
// - v : m x 1
// Output
// - avg: scalar
func Average(v []float64) float64 {
	avg := float64(0)
	for i := 0; i < len(v); i++ {
		avg += v[i]
	}
	avg = avg / float64(len(v))
	return avg
}
