package main

import (
	"fmt"
	"time"

	utils "github.com/CDSL-EncryptedControl/CDSL/utils"
	RGSW "github.com/CDSL-EncryptedControl/CDSL/utils/core/RGSW"
	RLWE "github.com/CDSL-EncryptedControl/CDSL/utils/core/RLWE"
	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func main() {
	// *****************************************************************
	// ************************* User's choice *************************
	// *****************************************************************
	// ============== Encryption parameters ==============
	// Refer to ``Homomorphic encryption standard''
	params, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		// log2 of polynomial degree
		LogN: 13,
		// Size of ciphertext modulus (Q)
		LogQ: []int{56},
		// Size of plaintext modulus (P)
		LogP:    []int{51},
		NTTFlag: true,
	})
	fmt.Println("Degree of polynomials:", params.N())
	fmt.Println("Ciphertext modulus:", params.QBigInt())
	fmt.Println("Special modulus:", params.PBigInt())
	// Default secret key distribution
	// Each coefficient in the polynomial is uniformly sampled in [-1, 0, 1]
	fmt.Println("Secret key distribution:", params.Xs())
	// Default error distribution
	// Each coefficient in the polynomial is sampled according to a
	// discrete Gaussian distribution with standard deviation 3.2 and bound 19.2
	fmt.Println("Error distribution:", params.Xe())

	// ============== Plant model ==============
	A := [][]float64{
		{0.2992, -0.1606, -0.8090, -0.5803},
		{0.2175, -0.3428, 0.3799, -0.3815},
		{-2.2822, -1.3596, -0.9895, -0.2337},
		{-0.6459, -0.0521, -0.9160, -0.1337},
	}
	B := [][]float64{
		{0.4652, -0.6905},
		{-0.0900, 0.3133},
		{0.1042, 0.1570},
		{0.0788, 0.0457},
	}
	C := [][]float64{
		{-1.1658, 0.1679, -0.2650, 0.1867},
		{-0.2356, -0.2303, -0.1040, -0.1006},
	}
	// Plant initial state
	xp_ini := []float64{
		1,
		1,
		1,
		1,
	}

	// ============== Pre-designed controller ==============
	// F must be an integer matrix
	F := [][]float64{
		{1, 1, 0, 0},
		{2, 0, 1, 0},
		{-1, 0, 0, 1},
		{2, 0, 0, 0},
	}
	G := [][]float64{
		{0, 1},
		{1, 3},
		{0, 5},
		{2, 1},
	}
	H := [][]float64{
		{1, 2, 3, 4},
		{0, 0, 1, 2},
	}
	// Controller initial state
	x_ini := []float64{
		1,
		2,
		-1,
		1,
	}

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

	// Generate keys
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	// Define encryptor and evaluator
	encryptorRLWE := rlwe.NewEncryptor(params, sk)
	decryptorRLWE := rlwe.NewDecryptor(params, sk)
	encryptorRGSW := rgsw.NewEncryptor(params, sk)
	evaluator := rgsw.NewEvaluator(params, evk)

	// ==============  Encryption of controller ==============
	// Quantization
	GBar := utils.ScalMatMult(1/s, G)
	HBar := utils.ScalMatMult(1/s, H)

	// Encryption
	ctF := RGSW.Enc(F, encryptorRGSW, levelQ, levelP, params)
	ctG := RGSW.Enc(GBar, encryptorRGSW, levelQ, levelP, params)
	ctH := RGSW.Enc(HBar, encryptorRGSW, levelQ, levelP, params)

	// ============== Simulation ==============
	// Number of simulation steps
	iter := 200
	fmt.Printf("Number of iterations: %v", iter)

	// *****************
	// 1) Plant + unencrypted (original) controller
	// *****************

	// Data storage
	yUnenc := [][]float64{}
	uUnenc := [][]float64{}
	xcUnenc := [][]float64{}
	xpUnenc := [][]float64{}

	xpUnenc = append(xpUnenc, xp_ini)
	xcUnenc = append(xcUnenc, x_ini)

	// Plant state
	xp := xp_ini
	// Controller state
	x := x_ini

	for i := 0; i < iter; i++ {
		y := utils.MatVecMult(C, xp)
		u := utils.MatVecMult(H, x)
		xp = utils.VecAdd(utils.MatVecMult(A, xp), utils.MatVecMult(B, u))
		x = utils.VecAdd(utils.MatVecMult(F, x), utils.MatVecMult(G, y))

		yUnenc = append(yUnenc, y)
		uUnenc = append(uUnenc, u)
		xcUnenc = append(xcUnenc, x)
		xpUnenc = append(xpUnenc, xp)
	}

	// *****************
	// 2) Plant + encrypted controller
	// *****************

	// State and output storage
	yEnc := [][]float64{}
	uEnc := [][]float64{}
	xpEnc := [][]float64{}
	xpEnc = append(xpEnc, xp_ini)

	// Plant state
	xp = xp_ini

	// Controller state encryption
	xBar := utils.ScalVecMult(1/(r*s), x_ini)
	xCt := RLWE.Enc(xBar, 1/L, *encryptorRLWE, ringQ, params)

	// For time check
	period := make([][]float64, iter)
	startPeriod := make([]time.Time, iter)

	for i := 0; i < iter; i++ {
		// **** Sensor ****
		// Plant output
		y := utils.MatVecMult(C, xp)

		startPeriod[i] = time.Now()

		// Quantize and encrypt plant output
		yBar := utils.RoundVec(utils.ScalVecMult(1/r, y))
		yCt := RLWE.Enc(yBar, 1/L, *encryptorRLWE, ringQ, params)

		// **** Encrypted Controller ****
		// Compute output
		uCt := RGSW.Mult(xCt, ctH, evaluator, ringQ, params)

		// **** Actuator ****
		// Decrypt output
		u := RLWE.Dec(uCt, *decryptorRLWE, r*s*s*L, ringQ, params)

		// **** Encrypted Controller ****
		// State update
		FxCt := RGSW.Mult(xCt, ctF, evaluator, ringQ, params)
		GyCt := RGSW.Mult(yCt, ctG, evaluator, ringQ, params)
		xCt = RLWE.AddVec(FxCt, GyCt, params)

		period[i] = []float64{float64(time.Since(startPeriod[i]).Microseconds()) / 1000}

		// **** Plant ****
		// State update
		xp = utils.VecAdd(utils.MatVecMult(A, xp), utils.MatVecMult(B, u))

		// Save data
		yEnc = append(yEnc, y)
		uEnc = append(uEnc, u)
		xpEnc = append(xpEnc, xp)
	}

	avgPeriod := utils.Average(utils.MatToVec(period))
	fmt.Println("Average elapsed time for a control period:", avgPeriod, "ms")

	// Compare plant input between 1) and 2)
	uDiff := make([][]float64, iter)
	for i := range uDiff {
		uDiff[i] = []float64{utils.Vec2Norm(utils.VecSub(uUnenc[i], uEnc[i]))}
	}

	// =========== Export data ===========

	// Plant state equipped with encrypted controller
	utils.DataExport(xpEnc, "./state.csv")

	// Plant intput from encrypted controller
	utils.DataExport(uEnc, "./uEnc.csv")

	// Plant output with encrypted controller
	utils.DataExport(yEnc, "./yEnc.csv")

	// Performance of encrypted controller
	utils.DataExport(uDiff, "./uDiff.csv")

	// Elapsed time
	utils.DataExport(period, "./period.csv")

}
