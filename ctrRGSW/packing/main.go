package main

import (
	"fmt"
	"math"

	"time"

	utils "github.com/CDSL-EncryptedControl/CDSL/utils"
	RGSW "github.com/CDSL-EncryptedControl/CDSL/utils/core/RGSW"
	RLWE "github.com/CDSL-EncryptedControl/CDSL/utils/core/RLWE"
	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
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

	A := [][]float64{
		{1, 0.0497869485895651, 0.00240013399977883, 3.99192178318363e-05},
		{0, 0.991480316476010, 0.0965381935679802, 0.00240013399977883},
		{0, -0.000612279081576232, 1.04210214304071, 0.0506998325801739},
		{0, -0.0246270901959133, 1.69541872243910, 1.04210214304071},
	}
	B := [][]float64{
		{0.00213051410434883},
		{0.0851968352399007},
		{0.00612279081576232},
		{0.246270901959133},
	}
	C := [][]float64{
		{1, 0, 0, 0},
	}
	// Plant initial state
	xp_ini := []float64{
		0,
		0,
		1,
		-1,
	}

	F := [][]float64{
		{1, 1, 0, 0, 0, 0, 0, 0},
		{13, 0, 1, 0, 0, 0, 0, 0},
		{4, 0, 0, 1, 0, 0, 0, 0},
		{-10, 0, 0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0},
		{0, 0, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 0},
	}
	// G := [][]float64{
	// 	{-640.464761022051},
	// 	{1715.36017949169},
	// 	{-1489.05705877506},
	// 	{389.491892102595},
	// 	{27.2282414240346},
	// 	{-0.604658260972883},
	// 	{-2.36399144578786},
	// 	{0.478391053561677},
	// }
	G := [][]float64{
		{-640.4648},
		{1715.3602},
		{-1489.0571},
		{389.4919},
		{27.2282},
		{-0.6047},
		{-2.3640},
		{0.4784},
	}

	H := [][]float64{
		{10, 0, 0, 0, 0, 0, 0, 0},
	}
	// Controller initial state
	x_ini := []float64{
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
	}

	// dimensions
	n := len(F)
	m := len(H)
	p := len(G[0])

	// ============== Quantization parameters ==============
	s := 1 / 10000.0
	L := 1 / 100000.0
	r := 1 / 100000.0
	fmt.Printf("Scaling parameters 1/L: %v, 1/s: %v, 1/r: %v \n", 1/L, 1/s, 1/r)
	// *****************************************************************
	// *****************************************************************

	// ============== Encryption settings ==============
	// Set parameters
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQ := params.RingQ()

	// Compute tau
	// least power of two greater than n, p_, and m
	maxDim := math.Max(math.Max(float64(n), float64(m)), float64(p))
	tau := int(math.Pow(2, math.Ceil(math.Log2(maxDim))))

	// Generate monomials for unpack
	logn := int(math.Log2(float64(tau)))
	monomials := make([]ring.Poly, logn)
	for i := 0; i < logn; i++ {
		monomials[i] = ringQ.NewPoly()
		idx := params.N() - params.N()/(1<<(i+1))
		monomials[i].Coeffs[0][idx] = 1
		ringQ.MForm(monomials[i], monomials[i])
		ringQ.NTT(monomials[i], monomials[i])
	}

	// tau = 8
	// galEls = 3개짜리
	// galEls = [8+1, 4+1, 2+1]

	// Generate Galois elements for unpack
	galEls := make([]uint64, int(math.Log2(float64(tau))))
	for i := 0; i < int(math.Log2(float64(tau))); i++ {
		galEls[i] = uint64(tau/int(math.Pow(2, float64(i))) + 1)
	}

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
	// Quantization
	GBar := utils.ScalMatMult(1/s, G)
	// HBar := utils.ScalMatMult(1/s, H)

	// Encryption
	// Dimension: 1-by-(# of columns)
	ctF := RGSW.EncPack(F, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctG := RGSW.EncPack(GBar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	// ctH := RGSW.EncPack(HBar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctH := RGSW.EncPack(H, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	// ============== Simulation ==============
	// Number of simulation steps
	iter := 100
	fmt.Printf("Number of iterations: %v\n", iter)

	// *****************
	// 1) Plant + unencrypted (original) controller
	// *****************

	// State and output storage
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

	// Dimension: 1-by-(# of elements)
	xBar := utils.ScalVecMult(1/(r*s), x_ini)
	xCtPack := RLWE.EncPack(xBar, tau, 1/L, *encryptorRLWE, ringQ, params)

	// For time check
	period := make([][]float64, iter)
	startPeriod := make([]time.Time, iter)

	for i := 0; i < iter; i++ {
		// **** Sensor ****
		// Plant output
		y := utils.MatVecMult(C, xp)

		startPeriod[i] = time.Now()

		// Quantize and encrypt
		yBar := utils.RoundVec(utils.ScalVecMult(1/r, y))
		yCtPack := RLWE.EncPack(yBar, tau, 1/L, *encryptorRLWE, ringQ, params)

		// **** Encrypted Controller ****
		// Unpack state
		xCt := RLWE.UnpackCt(xCtPack, n, tau, evaluatorRLWE, ringQ, monomials, params)

		// Unpack input
		yCt := RLWE.UnpackCt(yCtPack, p, tau, evaluatorRLWE, ringQ, monomials, params)

		// Compute output
		uCtPack := RGSW.MultPack(xCt, ctH, evaluatorRGSW, ringQ, params)

		// **** Actuator ****
		// Decrypt and Unapck
		u := RLWE.DecUnpack(uCtPack, m, tau, *decryptorRLWE, r*s*L, ringQ, params)

		// **** Encrypted Controller ****
		// State update
		FxCt := RGSW.MultPack(xCt, ctF, evaluatorRGSW, ringQ, params)
		GyCt := RGSW.MultPack(yCt, ctG, evaluatorRGSW, ringQ, params)
		xCtPack = RLWE.Add(FxCt, GyCt, params)

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

	// Export data ===============================================================

	// =========== Export data ===========

	// Plant state equipped with encrypted controller
	utils.DataExport(xpEnc, "./state.csv")

	// Plant intput from encrypted controller
	utils.DataExport(uEnc, "./uEnc.csv")

	// Plant intput from encrypted controller
	utils.DataExport(uUnenc, "./uUnenc.csv")

	// Plant output with encrypted controller
	utils.DataExport(yEnc, "./yEnc.csv")

	// Performance of encrypted controller
	utils.DataExport(uDiff, "./uDiff.csv")

	// Elapsed time
	utils.DataExport(period, "./period.csv")
}
