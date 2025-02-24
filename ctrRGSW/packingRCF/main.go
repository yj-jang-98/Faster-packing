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

	/// ============== Plant model ==============
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
	tau := int(math.Pow(2, math.Ceil(math.Log2(float64(m)))))

	monomial := ringQ.NewPoly()
	monomial.Coeffs[0][params.N()-params.N()/n] = params.Q()[0] - 1
	ringQ.MForm(monomial, monomial)
	ringQ.NTT(monomial, monomial)

	// Generate Galois elements for unpack
	// n+1 n/2 + 1 ... 2 + 1
	galEls := make([]uint64, int(math.Log2(float64(n*tau))))
	for i := 0; i < int(math.Log2(float64(n*tau))); i++ {
		galEls[i] = uint64((n*tau)/int(math.Pow(2, float64(i))) + 1)
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
	ctF := RGSW.EncF(F, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctG := RGSW.EncG(GBar, encryptorRGSW, levelQ, levelP, ringQ, params)
	// ctH := RGSW.EncH(HBar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctH := RGSW.EncH(H, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
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
	xCtPack := RLWE.EncPack(xBar, n, 1/L, *encryptorRLWE, ringQ, params)

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
		yCtPack := RLWE.Ency(yBar, 1/L, *encryptorRLWE, ringQ, params)

		// **** Encrypted Controller ****

		// Compute output
		xTraceOutput := RLWE.Trace(xCtPack, n*tau, n, evaluatorRLWE, ringQ, params)
		uCtPack := RGSW.Mult(xTraceOutput, ctH, evaluatorRGSW, ringQ, params)

		// **** Actuator ****
		// Decrypt and Unapck
		u := RLWE.DecUnpack(uCtPack, m, n*tau, *decryptorRLWE, r*s*L, ringQ, params)

		// **** Encrypted Controller ****
		// State update
		xTraceState := RLWE.Trace(xCtPack, n, 1, evaluatorRLWE, ringQ, params)
		FxCt := RGSW.MultFx(xTraceState, xCtPack, ctF, monomial, evaluatorRGSW, ringQ, params)
		GyCt := RGSW.Mult(yCtPack, ctG, evaluatorRGSW, ringQ, params)
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

	// Plant output with encrypted controller
	utils.DataExport(yEnc, "./yEnc.csv")

	// Performance of encrypted controller
	utils.DataExport(uDiff, "./uDiff.csv")

	// Elapsed time
	utils.DataExport(period, "./period.csv")
}
