package RGSW

import (
	"fmt"

	"github.com/CDSL-EncryptedControl/CDSL/utils"
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
func MultNaive(ctRLWE []*rlwe.Ciphertext, ctRGSW [][]*rgsw.Ciphertext, evaluator *rgsw.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) []*rlwe.Ciphertext {
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

// Computes external product with PACKING
// A x b = c
// ctRGSW[i] : RGSW encryption of i-th column of A
// ctRLWE[i] : RLWE encryption of i-th component of b
// ctOut     : RLWE encryption of Pack(c)
//
// Input
// - ctRGSW: n x 1 RGSW vector
// - ctRLWE: n x l RLWE vector
// Output
// - ctOut: RLWE ciphertext
func MultPack(ctRLWE []*rlwe.Ciphertext, ctRGSW []*rgsw.Ciphertext, evaluatorRGSW *rgsw.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) *rlwe.Ciphertext {
	row := len(ctRGSW)
	ctOut := rlwe.NewCiphertext(params, ctRLWE[0].Degree(), ctRLWE[0].Level())
	tmpCt := rlwe.NewCiphertext(params, ctRLWE[0].Degree(), ctRLWE[0].Level())
	for r := 0; r < row; r++ {
		evaluatorRGSW.ExternalProduct(ctRLWE[r], ctRGSW[r], tmpCt)
		ringQ.Add(ctOut.Value[0], tmpCt.Value[0], ctOut.Value[0])
		ringQ.Add(ctOut.Value[1], tmpCt.Value[1], ctOut.Value[1])
	}

	return ctOut
}

// Encrypts float matrix into an RGSW matrix
// Input
// - M: m x n float matrix
// Output
// - ctOut: m x n RGSW matrix
func Enc(M [][]float64, encryptor *rgsw.Encryptor, levelQ int, levelP int, params rlwe.Parameters) [][]*rgsw.Ciphertext {
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

func EncPackFx(F [][]float64, n int, encryptorRGSW *rgsw.Encryptor, levelQ int, levelP int, ringQ *ring.Ring, params rlwe.Parameters) *rgsw.Ciphertext {
	Fbar := make([][]float64, n)
	for r := 0; r < n; r++ {
		Fbar[r] = make([]float64, 1)
		if r < n-1 {
			Fbar[r][0] = F[r][0]
		} else {
			Fbar[r][0] = F[r][0] + 1
		}
	}
	fmt.Println(Fbar)
	modFbar := utils.ModMatFloat(Fbar, params.Q()[0])

	ctOut := rgsw.NewCiphertext(params, levelQ, levelP, 0)
	pt := rlwe.NewPlaintext(params, params.MaxLevel())
	for r := 0; r < n; r++ {
		// Store in the packing slots
		pt.Value.Coeffs[0][params.N()*r/n] = modFbar[r][0]
	}
	ringQ.NTT(pt.Value, pt.Value)
	encryptorRGSW.Encrypt(pt, ctOut)
	return ctOut
}

// Encrypts float matrix into a RGSW vector by PACKING each column
// Input
// - M: m x n float matrix
// Output
// - ctOut: n x 1 RGSW vector
// * ctOut[i] = RGSW encryption of i-th column of M
func EncPack(M [][]float64, tau int, encryptorRGSW *rgsw.Encryptor, levelQ int, levelP int, ringQ *ring.Ring, params rlwe.Parameters) []*rgsw.Ciphertext {
	row := len(M)
	col := len(M[0])
	modM := utils.ModMatFloat(M, params.Q()[0])

	ctOut := make([]*rgsw.Ciphertext, col)
	for c := 0; c < col; c++ {
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		for r := 0; r < row; r++ {
			// Store in the packing slots
			pt.Value.Coeffs[0][params.N()*r/tau] = modM[r][c]
		}
		ringQ.NTT(pt.Value, pt.Value)
		ctOut[c] = rgsw.NewCiphertext(params, levelQ, levelP, 0)
		encryptorRGSW.Encrypt(pt, ctOut[c])
	}
	return ctOut
}

// ================== New functions for CDC

func EncF(F [][]float64, encryptorRGSW *rgsw.Encryptor, levelQ int, levelP int, ringQ *ring.Ring, params rlwe.Parameters) *rgsw.Ciphertext {
	n := len(F)
	Fbar := make([][]float64, n)
	for r := 0; r < n; r++ {
		Fbar[r] = make([]float64, 1)
		if r < n-1 {
			Fbar[r][0] = F[r][0]
		} else {
			Fbar[r][0] = F[r][0] + 1
		}
	}
	modFbar := utils.ModMatFloat(Fbar, params.Q()[0])

	ctOut := rgsw.NewCiphertext(params, levelQ, levelP, 0)
	pt := rlwe.NewPlaintext(params, params.MaxLevel())
	for r := 0; r < n; r++ {
		// Store in the packing slots
		pt.Value.Coeffs[0][params.N()*r/n] = modFbar[r][0]
	}
	ringQ.NTT(pt.Value, pt.Value)
	encryptorRGSW.Encrypt(pt, ctOut)
	return ctOut
}

func EncH(H [][]float64, tau int, encryptorRGSW *rgsw.Encryptor, levelQ int, levelP int, ringQ *ring.Ring, params rlwe.Parameters) *rgsw.Ciphertext {

	m := len(H)
	n := len(H[0])
	modHbar := utils.ModMatFloat(H, params.Q()[0])
	ctOut := rgsw.NewCiphertext(params, levelQ, levelP, 0)
	pt := rlwe.NewPlaintext(params, params.MaxLevel())
	for c := 0; c < n; c++ {
		for r := 0; r < m; r++ {
			// Store in the packing slots
			if c == 0 {
				// First row
				pt.Value.Coeffs[0][r*params.N()/(n*tau)] = modHbar[r][c]
			} else {
				// Rest of the rows
				pt.Value.Coeffs[0][params.N()-c*params.N()/n+r*params.N()/(n*tau)] = params.Q()[0] - modHbar[r][c]
			}
		}
	}
	ringQ.NTT(pt.Value, pt.Value)
	encryptorRGSW.Encrypt(pt, ctOut)
	return ctOut
}

func EncG(G [][]float64, encryptorRGSW *rgsw.Encryptor, levelQ int, levelP int, ringQ *ring.Ring, params rlwe.Parameters) *rgsw.Ciphertext {
	n := len(G)
	p := len(G[0])

	modGbar := utils.ModMatFloat(G, params.Q()[0])

	ctOut := rgsw.NewCiphertext(params, levelQ, levelP, 0)
	pt := rlwe.NewPlaintext(params, params.MaxLevel())
	for r := 0; r < n; r++ {
		for c := 0; c < p; c++ {
			// Store in the packing slots
			if r == 0 {
				if c == 0 {
					// g_{0,0} 상수항
					pt.Value.Coeffs[0][0] = modGbar[r][c]
				} else {
					// -g_{0,1} ... -g_{0,p-1}
					pt.Value.Coeffs[0][params.N()-c] = params.Q()[0] - modGbar[r][c]
				}
			} else {
				// g_{i,j}
				pt.Value.Coeffs[0][r*params.N()/n-c] = modGbar[r][c]
			}
		}
	}
	ringQ.NTT(pt.Value, pt.Value)
	encryptorRGSW.Encrypt(pt, ctOut)
	return ctOut
}

func MultFx(ctTrace *rlwe.Ciphertext, ctRLWE *rlwe.Ciphertext, ctRGSW *rgsw.Ciphertext, monomial ring.Poly, evaluatorRGSW *rgsw.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) *rlwe.Ciphertext {
	// ctRGSW와 external product
	// ctRLWE에 X^{N/n} monomial 곱셈
	// 두 결과 덧셈

	ctOut := rlwe.NewCiphertext(params, ctRLWE.Degree(), ctRLWE.Level())
	ctTmp := rlwe.NewCiphertext(params, ctRLWE.Degree(), ctRLWE.Level())

	// // Monomial multiplication
	ringQ.MulCoeffsMontgomery(ctRLWE.Value[0], monomial, ctTmp.Value[0])
	ringQ.MulCoeffsMontgomery(ctRLWE.Value[1], monomial, ctTmp.Value[1])

	// External product
	evaluatorRGSW.ExternalProduct(ctTrace, ctRGSW, ctOut)

	ringQ.Add(ctOut.Value[0], ctTmp.Value[0], ctOut.Value[0])
	ringQ.Add(ctOut.Value[1], ctTmp.Value[1], ctOut.Value[1])

	return ctOut
}

func Mult(ctRLWE *rlwe.Ciphertext, ctRGSW *rgsw.Ciphertext, evaluatorRGSW *rgsw.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) *rlwe.Ciphertext {
	ctOut := rlwe.NewCiphertext(params, ctRLWE.Degree(), ctRLWE.Level())

	// External product
	evaluatorRGSW.ExternalProduct(ctRLWE, ctRGSW, ctOut)
	return ctOut
}
