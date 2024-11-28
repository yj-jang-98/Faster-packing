package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
)

func main() {

	// ************************* User's choice *************************
	// Encryption parameters
	// Refer to ``Homomorphic encryption standard''

	// log2 of polynomial degree
	logN := 12
	// Choose the size of plaintext modulus (2^ptSize)
	ptSize := uint64(26)
	// Choose the size of ciphertext modulus (2^ctSize)
	ctSize := int(74)

	// Plant model
	A := [][]float64{
		{1, 0.0020, 0.0663, 0.0047, 0.0076},
		{0, 1.0077, 2.0328, -0.5496, -0.0591},
		{0, 0.0478, 0.9850, -0.0205, -0.0092},
		{0, 0, 0, 0.3679, 0},
		{0, 0, 0, 0, 0.3679},
	}
	B := [][]float64{
		{0.0029, 0.0045},
		{-0.3178, -0.0323},
		{-0.0086, -0.0051},
		{0.6321, 0},
		{0, 0.6321},
	}
	C := [][]float64{
		{0, 1, 0, 0, 0},
		{0, -0.2680, 47.7600, -4.5600, 4.4500},
		{1, 0, 0, 0, 0},
		{0, 0, 0, 1, 0},
		{0, 0, 0, 0, 1},
	}
	// Plant initial state
	xp0 := [][]float64{
		{1},
		{-1},
		{0},
		{0.7},
		{1},
	}

	// Pre-designed controller
	F := [][]float64{
		{0.4064, 0.0014, 0.0004, 0.0053, 0.0009},
		{-0.1216, 0.2237, -1.0968, -0.0641, -0.1033},
		{0.0046, 0.0175, -0.0122, 0.0411, -0.0521},
		{0.3629, 0.3552, 2.3131, -0.0367, -0.0411},
		{-1.1877, -0.1964, -1.4755, 0.1370, 0.2368},
	}
	G := [][]float64{
		{0.0011, 0.0014, 0.5868, 0.0056, 0.0007},
		{0.6296, 0.0429, -0.0003, -0.1811, -0.1278},
		{0.0326, 0.0205, 0, 0.0337, -0.0480},
		{-0.0049, -0.0003, 0.0002, 0.1732, 0.0005},
		{-0.0037, 0.0003, 0, 0.0005, 0.1733},
	}
	H := [][]float64{
		{0.5743, 0.5544, 3.6332, -0.3636, -0.0668},
		{-1.8788, -0.3166, -2.3100, 0.2151, 0.0691},
	}

	// input-output representation (obtained by MATLAB)
	// u(k) = Hy[5]*y(k-1) + ... + Hy[1]*y(k-5) + Hu[5]*u(k-1) + ... + Hu[1]*u(k-5)
	// already vectorized
	Hy := [][]float64{
		{0.0003, 0.0009, -0.0124, -0.0001, 0.0001,
			-0.0006, -0.0006, 0.0241, 0.0002, -0.0002},
		{-0.0175, -0.0072, 0.0245, -0.0048, 0.0124,
			0.0131, 0.0049, -0.0480, 0.0022, -0.0094},
		{0.0639, 0.0144, 0.1916, -0.0100, -0.0333,
			-0.0388, -0.0088, -0.2579, 0.0038, 0.0198},
		{-0.0038, -0.0264, 0.0647, -0.0383, 0.0176,
			-0.0006, 0.0143, -0.4155, 0.0164, -0.0077},
		{0.4702, 0.0992, 0.3368, -0.0378, -0.2566,
			-0.2780, -0.0636, -1.1023, 0.0063, 0.1621},
	}
	Hu := [][]float64{
		{-0.0037, 0.0138, 0, 0, 0,
			-0.0021, -0.0177, 0, 0, 0},
		{-0.0118, -0.0463, 0, 0, 0,
			0.0265, 0.0612, 0, 0, 0},
		{-0.0117, -0.0265, 0, 0, 0,
			0.0105, 0.0425, 0, 0, 0},
		{0.0063, -0.0198, 0, 0, 0,
			-0.0085, 0.0285, 0, 0, 0},
		{0.0038, -0.0094, 0, 0, 0,
			-0.0071, 0.0147, 0, 0, 0},
	}
	// Controller initial state
	x0 := [][]float64{
		{-0.001},
		{0.013},
		{0.2},
		{-0.02},
		{0},
	}
	// initial [y(-n); ...; y(-1)] (obtained by MATLAB)
	yy0 := [][]float64{
		{0.1583, 0.0337, 0.4869, -0.0153, -0.0869},
		{0.0243, -0.0310, 0.3190, -0.0329, -0.0066},
		{0.6712, 0.1841, 1.0240, 0.0609, -0.4492},
		{-0.6974, -0.1527, -0.9954, 0.2786, 0.1280},
		{0.0887, 0.6957, 0.1828, 1.1950, -1.3772},
	}
	// initial [u(-n); ...; u(-1)] (obtained by MATLAB)
	uu0 := [][]float64{
		{0, 0},
		{0.2646, -0.5970},
		{0.1544, -0.5698},
		{0.9405, -1.6955},
		{-0.5578, 0.7293},
	}

	// Quantization parameters
	L := 0.00050
	s := 0.00010
	fmt.Println("Scaling parameters 1/L:", 1/L, "1/s:", 1/s)

	// ****************************************************************************************************

	// =========== Encryption settings ===========
	// Search a proper prime to set plaintext modulus
	primeGen := ring.NewNTTFriendlyPrimesGenerator(ptSize, uint64(math.Pow(2, float64(logN)+1)))
	ptModulus, _ := primeGen.NextAlternatingPrime()
	fmt.Println("Plaintext modulus:", ptModulus)

	// Create a chain of ciphertext modulus
	logQ := []int{int(math.Floor(float64(ctSize) * 0.5)), int(math.Ceil(float64(ctSize) * 0.5))}

	// Parameters satisfying 128-bit security
	// BGV scheme is used
	params, _ := bgv.NewParametersFromLiteral(bgv.ParametersLiteral{
		LogN:             logN,
		LogQ:             logQ,
		PlaintextModulus: ptModulus,
	})
	fmt.Println("Ciphertext modulus:", params.QBigInt())
	fmt.Println("Degree of polynomials:", params.N())

	// Generate secret key
	kgen := bgv.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()

	encryptor := bgv.NewEncryptor(params, sk)
	decryptor := bgv.NewDecryptor(params, sk)
	encoder := bgv.NewEncoder(params)
	eval := bgv.NewEvaluator(params, nil)

	bredparams := ring.GenBRedConstant(params.PlaintextModulus())

	// =========== Encryption of controller ===========
	// dimensions
	nx := len(A)
	ny := len(C)
	nu := len(B[0])
	h := ny

	// duplicate
	yy0vec := make([][]float64, nx)
	uu0vec := make([][]float64, nx)
	for i := 0; i < nx; i++ {
		yy0vec[i] = utils.vecDuplicate(yy0[i], nu, h)
		uu0vec[i] = utils.vecDuplicate(uu0[i], nu, h)
	}

	// Plaintext of past inputs and outputs
	ptY := make([]*rlwe.Plaintext, nx)
	ptU := make([]*rlwe.Plaintext, nx)
	// Plaintext of control parameters
	ptHy := make([]*rlwe.Plaintext, nx)
	ptHu := make([]*rlwe.Plaintext, nx)
	// Ciphertext of past inputs and outputs
	ctY := make([]*rlwe.Ciphertext, nx)
	ctU := make([]*rlwe.Ciphertext, nx)
	// Ciphertext of control parameters
	ctHy := make([]*rlwe.Ciphertext, nx)
	ctHu := make([]*rlwe.Ciphertext, nx)

	// Quantization - packing - encryption
	for i := 0; i < nx; i++ {
		ptY[i] = bgv.NewPlaintext(params, params.MaxLevel())
		encoder.Encode(modVecFloat(roundVec(scalarVecMult(1/L, yy0vec[i])), params.PlaintextModulus()), ptY[i])
		ctY[i], _ = encryptor.EncryptNew(ptY[i])

		ptU[i] = bgv.NewPlaintext(params, params.MaxLevel())
		encoder.Encode(modVecFloat(roundVec(scalarVecMult(1/L, uu0vec[i])), params.PlaintextModulus()), ptU[i])
		ctU[i], _ = encryptor.EncryptNew(ptU[i])

		ptHy[i] = bgv.NewPlaintext(params, params.MaxLevel())
		encoder.Encode(modVecFloat(roundVec(scalarVecMult(1/s, Hy[i])), params.PlaintextModulus()), ptHy[i])
		ctHy[i], _ = encryptor.EncryptNew(ptHy[i])

		ptHu[i] = bgv.NewPlaintext(params, params.MaxLevel())
		encoder.Encode(modVecFloat(roundVec(scalarVecMult(1/s, Hu[i])), params.PlaintextModulus()), ptHu[i])
		ctHu[i], _ = encryptor.EncryptNew(ptHu[i])
	}

	// =========== Simulation ===========
	// Number of simulation steps
	nsim := 100

	// 1) Plant + unencrypted (original) controller
	// Plant state
	xp := xp0
	// Controller state
	xc := x0
	// To save data
	yUnenc := [][]float64{}
	uUnenc := [][]float64{}
	xpUnenc := [][]float64{}
	xcUnenc := [][]float64{}
	xpUnenc = appendVecToMat(xpUnenc, xp0)
	xcUnenc = appendVecToMat(xcUnenc, x0)

	for i := 0; i < nsim; i++ {
		y := matMult(C, xp)
		u := matMult(H, xc)
		xp = matAdd(matMult(A, xp), matMult(B, u))
		xc = matAdd(matMult(F, xc), matMult(G, y))

		yUnenc = appendVecToMat(yUnenc, y)
		uUnenc = appendVecToMat(uUnenc, u)
		xpUnenc = appendVecToMat(xpUnenc, xp)
		xcUnenc = appendVecToMat(xcUnenc, xc)
	}

	// 2) Plant + encrypted controller
	// Plant state
	XP := xp0

	// Plant input
	U := make([][]float64, nu)
	// Unpacked and re-scaled u at actuator
	Uact := make([]uint64, params.N())
	// u after inner sum
	Usum := make([]uint64, nu)

	// To save data
	yEnc := [][]float64{}
	uEnc := [][]float64{}
	xpEnc := [][]float64{}
	xpEnc = appendVecToMat(xpEnc, xp0)

	// For time check
	period := make([][]float64, nsim)
	startPeriod := make([]time.Time, nsim)

	for i := 0; i < nsim; i++ {
		// Plant output
		Y := matMult(C, XP) // [][]float64

		startPeriod[i] = time.Now()

		// Sensor
		// Quantize and duplicate
		Ysens := modVecFloat(roundVec(scalarVecMult(1/L, vecDuplicate(mat2vec(Y), nu, h))), params.PlaintextModulus())
		Ypacked := bgv.NewPlaintext(params, params.MaxLevel())
		encoder.Encode(Ysens, Ypacked)
		Ycin, _ := encryptor.EncryptNew(Ypacked)

		// Encrypted controller
		Uout, _ := eval.MulNew(ctHy[0], ctY[0])
		eval.MulThenAdd(ctHu[0], ctU[0], Uout)
		for j := 1; j < nx; j++ {
			eval.MulThenAdd(ctHy[j], ctY[j], Uout)
			eval.MulThenAdd(ctHu[j], ctU[j], Uout)
		}

		// Actuator
		encoder.Decode(decryptor.DecryptNew(Uout), Uact)
		// Generate plant input
		for k := 0; k < nu; k++ {
			Usum[k] = vecSumUint(Uact[k*h:(k+1)*h], params.PlaintextModulus(), bredparams)
			U[k] = []float64{L * s * signFloat(float64(Usum[k]), params.PlaintextModulus())}
		}
		// Re-encryption
		Upacked := bgv.NewPlaintext(params, params.MaxLevel())
		encoder.Encode(modVecFloat(roundVec(scalarVecMult(1/L, vecDuplicate(mat2vec(U), nu, h))), params.PlaintextModulus()), Upacked)
		Ucin, _ := encryptor.EncryptNew(Upacked)

		// Controller state update
		ctY = append(ctY[1:], Ycin)
		ctU = append(ctU[1:], Ucin)

		period[i] = []float64{float64(time.Since(startPeriod[i]).Microseconds()) / 1000}

		// Plant state update
		XP = matAdd(matMult(A, XP), matMult(B, U))

		// Save data
		yEnc = appendVecToMat(yEnc, Y)
		uEnc = appendVecToMat(uEnc, U)
		xpEnc = appendVecToMat(xpEnc, XP)
	}

	avgPeriod := average(mat2vec(period))
	fmt.Println("Average elapsed time for a control period:", avgPeriod, "ms")

	// Compare plant input between 1) and 2)
	uDiff := make([][]float64, nsim)
	for i := range uDiff {
		uDiff[i] = []float64{vec2norm(subVec(uUnenc[i], uEnc[i]))}
	}

	// =========== Export data ===========

	// Plant state equipped with encrypted controller
	file1, err := os.Create("./state.csv")
	if err != nil {
		panic(err)
	}
	wr1 := csv.NewWriter(bufio.NewWriter(file1))
	wr1.WriteAll(mat2string(xpEnc))

	// Plant intput from encrypted controller
	file2, err := os.Create("./uEnc.csv")
	if err != nil {
		panic(err)
	}
	wr2 := csv.NewWriter(bufio.NewWriter(file2))
	wr2.WriteAll(mat2string(uEnc))

	// Plant output with encrypted controller
	file3, err := os.Create("./yEnc.csv")
	if err != nil {
		panic(err)
	}
	wr3 := csv.NewWriter(bufio.NewWriter(file3))
	wr3.WriteAll(mat2string(yEnc))

	// Performance of encrypted controller
	file4, err := os.Create("./uDiff.csv")
	if err != nil {
		panic(err)
	}
	wr4 := csv.NewWriter(bufio.NewWriter(file4))
	wr4.WriteAll(mat2string(uDiff))

	// Elapsed time
	file5, err := os.Create("./period.csv")
	if err != nil {
		panic(err)
	}
	wr5 := csv.NewWriter(bufio.NewWriter(file5))
	wr5.WriteAll(mat2string(period))
}
