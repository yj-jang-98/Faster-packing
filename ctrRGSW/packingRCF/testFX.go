package main

import (
	"fmt"
	"math"

	"time"

	utils "github.com/CDSL-EncryptedControl/CDSL/utils"
	RLWE "github.com/CDSL-EncryptedControl/CDSL/utils/core/RLWE"
	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

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

func Trace(ctRLWE *rlwe.Ciphertext, nStart int, nFinal int, evaluatorRLWE *rlwe.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) *rlwe.Ciphertext {
	// scale
	// Q - (Q+1)/tau
	scalar := params.Q()[0] - uint64((params.Q()[0]+1)/uint64(nStart))
	ctUnpack := rlwe.NewCiphertext(params, ctRLWE.Degree(), ctRLWE.Level())
	ringQ.MulScalar(ctRLWE.Value[0], scalar, ctUnpack.Value[0])
	ringQ.MulScalar(ctRLWE.Value[1], scalar, ctUnpack.Value[1])
	tmpCt := rlwe.NewCiphertext(params, ctRLWE.Degree(), ctRLWE.Level())
	for i := nStart; i > nFinal; i /= 2 {
		// Automorphism
		evaluatorRLWE.Automorphism(ctUnpack, uint64(i+1), tmpCt)

		ringQ.Add(ctUnpack.Value[0], tmpCt.Value[0], ctUnpack.Value[0])
		ringQ.Add(ctUnpack.Value[1], tmpCt.Value[1], ctUnpack.Value[1])
	}

	return ctUnpack
}

func main() {
	// *****************************************************************
	// ************************* User's choice *************************
	// *****************************************************************
	// ============== Encryption parameters ==============
	// Refer to ``Homomorphic encryption standard''
	params, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		// log2 of polynomial degree
		LogN: 12,
		// Size of ciphertext modulus (Q)
		LogQ: []int{56},
		// Size of special modulus (P)
		LogP:    []int{51},
		NTTFlag: true,
	})
	fmt.Println("Degree of polynomials:", params.N())
	fmt.Println("Ciphertext modulus:", params.QBigInt())
	fmt.Println("Special modulus:", params.PBigInt())
	// Default secret key distribution
	// Each coefficient in the polynomial is uniformly sampled in [-1, 0, 1]
	fmt.Println("Secret key distribution (Ternary):", params.Xs())
	// Default error distribution
	// Each coefficient in the polynomial is sampled according to a
	// discrete Gaussian distribution with standard deviation 3.2 and bound 19.2
	fmt.Println("Error distribution (Discrete Gaussian):", params.Xe())

	// ============== Pre-designed controller ==============
	// F must be an integer matrix
	F := [][]float64{
		{1, 1, 0, 0},
		{2, 0, 1, 0},
		{-1, 0, 0, 1},
		{2, 0, 0, 0},
	}

	// Controller initial state
	x_ini := []float64{
		1,
		2,
		-1,
		1,
	}
	// dimensions
	n := len(F)

	// ============== Quantization parameters ==============
	s := 1 / 1.0
	L := 1 / 1000000.0
	r := 1 / 10000.0
	fmt.Printf("Scaling parameters 1/L: %v, 1/s: %v, 1/r: %v \n", 1/L, 1/s, 1/r)
	// *****************************************************************
	// *****************************************************************

	// ============== Encryption settings ==============
	// Set parameters
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQ := params.RingQ()

	// Generate Galois elements for unpack
	galEls := make([]uint64, int(math.Log2(float64(n))))
	for i := 0; i < int(math.Log2(float64(n))); i++ {
		galEls[i] = uint64(n/int(math.Pow(2, float64(i))) + 1)
	}

	monomial := ringQ.NewPoly()
	monomial.Coeffs[0][params.N()-params.N()/n] = params.Q()[0] - 1
	ringQ.MForm(monomial, monomial)
	ringQ.NTT(monomial, monomial)

	// Generate keys
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evkRGSW := rlwe.NewMemEvaluationKeySet(rlk)
	evkRLWE := rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...)

	// Define encryptor and evaluator
	encryptorRLWE := rlwe.NewEncryptor(params, sk)
	decryptorRLWE := rlwe.NewDecryptor(params, sk)
	encryptorRGSW := rgsw.NewEncryptor(params, sk)
	evaluatorRGSW := rgsw.NewEvaluator(params, evkRGSW)
	evaluatorRLWE := rlwe.NewEvaluator(params, evkRLWE)

	// ==============  Encryption of controller ==============

	// Encryption
	ctF := EncPackFx(F, n, encryptorRGSW, levelQ, levelP, ringQ, params)

	// ============== Simulation ==============
	// Number of simulation steps
	iter := 3
	fmt.Printf("Number of iterations: %v\n", iter)

	// *****************
	// 1) unencrypted (original) controller
	// *****************

	xcUnenc := [][]float64{}
	xcUnenc = append(xcUnenc, x_ini)

	// Controller state
	x := x_ini
	fmt.Printf("Initial \n")
	fmt.Println(x)
	for i := 0; i < iter; i++ {
		x = utils.MatVecMult(F, x)
		xcUnenc = append(xcUnenc, x)
		fmt.Printf("iteration %v\n", i)
		fmt.Println(x)
	}

	// *****************
	// 2) encrypted controller
	// *****************

	// State and output storage
	xcEnc := [][]float64{}

	// Dimension: 1-by-(# of elements)
	xBar := utils.ScalVecMult(1/(r*s), x_ini)
	xCtPack := RLWE.EncPack(xBar, n, 1/L, *encryptorRLWE, ringQ, params)
	xc := RLWE.DecUnpack(xCtPack, n, n, *decryptorRLWE, r*s*L, ringQ, params)
	xcEnc = append(xcEnc, xc)
	fmt.Println("==========================")
	fmt.Printf("Initial \n")
	fmt.Println(xc)
	// For time check
	period := make([][]float64, iter)
	startPeriod := make([]time.Time, iter)

	for i := 0; i < iter; i++ {
		// FX
		startPeriod[i] = time.Now()
		xTrace := Trace(xCtPack, n, 1, evaluatorRLWE, ringQ, params)
		xCtPack = MultFx(xTrace, xCtPack, ctF, monomial, evaluatorRGSW, ringQ, params)
		period[i] = []float64{float64(time.Since(startPeriod[i]).Microseconds()) / 1000}
		xc = RLWE.DecUnpack(xCtPack, n, n, *decryptorRLWE, r*s*L, ringQ, params)
		fmt.Printf("iteration %v\n", i)
		// fmt.Println(xTraceDec)
		fmt.Println(xc)
		// Save data
		xcEnc = append(xcEnc, xc)
	}

	avgPeriod := utils.Average(utils.MatToVec(period))
	fmt.Println("Average elapsed time:", avgPeriod, "ms")

}
