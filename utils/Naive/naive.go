package naive

import (
	"math"

	"github.com/CDSL-EncryptedControl/2024SICE/utils"
	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// Computes external product
// Input
// - ctRGSW: m x n RGSW matrix
// - ctRLWE: n x l RLWE vector
// Output
// - ctOut: m x l RLWE vector
func Mult(ctRLWE []*rlwe.Ciphertext, ctRGSW [][]*rgsw.Ciphertext, evaluator *rgsw.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) []*rlwe.Ciphertext {
	row := len(ctRGSW)
	col := len(ctRGSW[0])
	ctOut := make([]*rlwe.Ciphertext, row)
	tmpCt := rlwe.NewCiphertext(params, ctRLWE[0].Degree(), ctRLWE[0].Level())
	for r := 0; r < row; r++ {
		ctOut[r] = rlwe.NewCiphertext(params, ctRLWE[0].Degree(), ctRLWE[0].Level())
		for c := 0; c < col; c++ {
			evaluator.ExternalProduct(ctRLWE[c], ctRGSW[r][c], tmpCt)
			ringQ.Add(ctOut[r].Value[0], tmpCt.Value[0], ctOut[r].Value[0])
			ringQ.Add(ctOut[r].Value[1], tmpCt.Value[1], ctOut[r].Value[1])
		}
	}
	return ctOut
}

// Encrypts float vector into an RLWE ciphertext
// Input
// - v: n x 1 float vector
// Output
// - ctOut: n x l RLWE vector
func EncRlwe(v []float64, scale float64, encryptor rlwe.Encryptor, params rlwe.Parameters) []*rlwe.Ciphertext {
	var err error

	row := len(v)
	scaleV := utils.ScalarVecMult(scale, v)
	modV := utils.ModVecFloat(scaleV, params.Q()[0])

	ctOut := make([]*rlwe.Ciphertext, row)
	for r := 0; r < row; r++ {
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		for i := 0; i < params.N(); i++ {
			pt.Value.Coeffs[0][i] = modV[r]
		}
		ctOut[r], err = encryptor.EncryptNew(pt)
		if err != nil {
			panic(err)
		}
	}

	return ctOut
}

// Encrypts float matrix into an RGSW matrix
// Input
// - M: m x n float matrix
// Output
// - ctOut: m x n RGSW matrix
func EncRgsw(M [][]float64, encryptor *rgsw.Encryptor, levelQ int, levelP int, params rlwe.Parameters) [][]*rgsw.Ciphertext {
	row := len(M)
	col := len(M[0])
	modM := utils.ModMatFloat(M, params.Q()[0])

	ctOut := make([][]*rgsw.Ciphertext, row)
	for r := 0; r < row; r++ {
		ctOut[r] = make([]*rgsw.Ciphertext, col)
		for c := 0; c < col; c++ {
			pt := rlwe.NewPlaintext(params, params.MaxLevel())
			for i := 0; i < params.N(); i++ {
				pt.Value.Coeffs[0][i] = modM[r][c]
			}
			ctOut[r][c] = rgsw.NewCiphertext(params, levelQ, levelP, 0)
			encryptor.Encrypt(pt, ctOut[r][c])
		}
	}
	return ctOut
}

// 1) Decrypt RLWE vector -> integer vector
// 2) Map constant terms of integer vector [0,q/2) -> [-q/2, q/2)
// 3) Scale down
// Input
// - ctRLWE: n x 1 RLWE vector
// Output
// - valOut: n x 1 float vector
func Dec(ctRLWE []*rlwe.Ciphertext, decryptor rlwe.Decryptor, scale float64, params rlwe.Parameters) []float64 {
	row := len(ctRLWE)
	q := float64(params.Q()[0])
	valOut := make([]float64, row)
	for r := 0; r < row; r++ {
		pt := decryptor.DecryptNew(ctRLWE[r])
		if pt.IsNTT {
			params.RingQ().INTT(pt.Value, pt.Value)
		}
		// Constant terms
		val := float64(pt.Value.Coeffs[0][0])
		// Mapping to [-q/2, q/2)
		val = val - math.Floor((val+q/2.0)/q)*q
		// Scale down
		valOut[r] = val * scale
	}
	return valOut
}

// Add RLWE vectors
// Input
// - ctRLWE1: n x 1 RLWE vector
// - ctRLWE2: n x 1 RLWE vector
// Output
// - ctOut: n x 1 RLWE vector
func Add(ctRLWE1 []*rlwe.Ciphertext, ctRLWE2 []*rlwe.Ciphertext, params rlwe.Parameters) []*rlwe.Ciphertext {
	row := len(ctRLWE1)
	ctOut := make([]*rlwe.Ciphertext, row)
	for r := 0; r < row; r++ {
		ctOut[r] = rlwe.NewCiphertext(params, ctRLWE2[0].Degree(), ctRLWE2[0].Level())
		params.RingQ().Add(ctRLWE1[r].Value[0], ctRLWE2[r].Value[0], ctOut[r].Value[0])
		params.RingQ().Add(ctRLWE1[r].Value[1], ctRLWE2[r].Value[1], ctOut[r].Value[1])
	}

	return ctOut
}
