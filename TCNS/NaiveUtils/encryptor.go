package naive

import (
	"math"

	"github.com/CDSL-EncryptedControl/2024SICE/utils"
	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func externalProduct(ctB []*rlwe.Ciphertext, ctA [][]*rgsw.Ciphertext, evaluator *rgsw.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) []*rlwe.Ciphertext {
	// Computes the external product between ctA and ctB
	// ctA: m x n RGSW ciphertexts matrix
	// ctB: n x l RLWE ciphertexts vector
	// ctC: m x l RLWE ciphertexts vector

	row := len(ctA)    // m
	col := len(ctA[0]) // n
	ctC := make([]*rlwe.Ciphertext, row)
	tmpCt := rlwe.NewCiphertext(params, ctB[0].Degree(), ctB[0].Level())
	for r := 0; r < row; r++ {
		ctC[r] = rlwe.NewCiphertext(params, ctB[0].Degree(), ctB[0].Level())
		for c := 0; c < col; c++ {
			evaluator.ExternalProduct(ctB[c], ctA[r][c], tmpCt)
			ringQ.Add(ctC[r].Value[0], tmpCt.Value[0], ctC[r].Value[0])
			ringQ.Add(ctC[r].Value[1], tmpCt.Value[1], ctC[r].Value[1])
		}
	}
	return ctC
}

func encryptRlwe(A []float64, scale float64, encryptor rlwe.Encryptor, params rlwe.Parameters) []*rlwe.Ciphertext {
	// Encrypts an n-dimensional float vector A into an n-dimensional RLWE ciphertexts vector ctA after scaling
	var err error

	row := len(A)
	ctA := make([]*rlwe.Ciphertext, row)

	// Scale up. Scale should be chosen so that A_ is a vector consisting of integers in [-q/2, q/2)
	A_ := utils.ScalarVecMult(scale, A)
	modA := utils.ModVecFloat(A_, params.Q()[0])

	for r := 0; r < row; r++ {
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		for i := 0; i < params.N(); i++ {
			pt.Value.Coeffs[0][i] = modA[r]
		}
		ctA[r], err = encryptor.EncryptNew(pt)
		if err != nil {
			panic(err)
		}
	}

	return ctA
}

func encryptRgsw(A [][]float64, encryptor *rgsw.Encryptor, levelQ int, levelP int, params rlwe.Parameters) [][]*rgsw.Ciphertext {
	// Encrypts an m-by-n-dimensional float matrix A into an m-by-n-dimensional RGSW ciphertexts matrix ctA

	row := len(A)
	col := len(A[0])
	ctA := make([][]*rgsw.Ciphertext, row)
	modA := utils.ModMatFloat(A, params.Q()[0])
	for r := 0; r < row; r++ {
		ctA[r] = make([]*rgsw.Ciphertext, col)
		for c := 0; c < col; c++ {
			pt := rlwe.NewPlaintext(params, params.MaxLevel())
			for i := 0; i < params.N(); i++ {
				pt.Value.Coeffs[0][i] = modA[r][c]
			}
			ctA[r][c] = rgsw.NewCiphertext(params, levelQ, levelP, 0)
			encryptor.Encrypt(pt, ctA[r][c])
		}
	}
	return ctA
}

func decryptRlwe(ctA []*rlwe.Ciphertext, decryptor rlwe.Decryptor, scale float64, params rlwe.Parameters) []float64 {
	// 1) Decrypts an n-dimensional RLWE vector ctA and obtain an n-dimensional integer vector pt
	// 2) Maps the constant terms of pt from the set [0,q/2) back to [-q/2, q/2)
	// 3) Scale down and return decA

	row := len(ctA)
	q := float64(params.Q()[0])
	decA := make([]float64, row)
	for r := 0; r < row; r++ {
		pt := decryptor.DecryptNew(ctA[r])
		if pt.IsNTT {
			params.RingQ().INTT(pt.Value, pt.Value)
		}
		// Constant terms
		val := float64(pt.Value.Coeffs[0][0])
		// Mapping to [-q/2, q/2)
		val = val - math.Floor((val+q/2.0)/q)*q
		// Scale down
		decA[r] = val * scale
	}
	return decA
}

func ctAdd(ctA []*rlwe.Ciphertext, ctB []*rlwe.Ciphertext, params rlwe.Parameters) []*rlwe.Ciphertext {
	// Adds two m-dimensional RLWE ciphertexts vector ctA and ctB
	// A : m x 1
	// B : m x 1

	row := len(ctA)
	ctC := make([]*rlwe.Ciphertext, row)
	for r := 0; r < row; r++ {
		ctC[r] = rlwe.NewCiphertext(params, ctB[0].Degree(), ctB[0].Level())
		params.RingQ().Add(ctA[r].Value[0], ctB[r].Value[0], ctC[r].Value[0])
		params.RingQ().Add(ctA[r].Value[1], ctB[r].Value[1], ctC[r].Value[1])
	}

	return ctC
}
