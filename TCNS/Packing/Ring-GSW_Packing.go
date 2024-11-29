package main

import (
	"fmt"
	"math"
	"os"
	"time"

	"bufio"
	"encoding/csv"

	"github.com/CDSL-EncryptedControl/2024SICE/utils"
	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"

	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
	"gonum.org/v1/plot/vg/vgimg"
)

func modZq(a [][]float64, params rlwe.Parameters) [][]float64 {
	// Components of the matrix 'a' belongs to [-q/2, q/2)
	// Takes the modulo operation and maps all components to [0,q)

	q := float64(params.Q()[0])
	b := make([][]float64, len(a))
	for i := 0; i < len(a); i++ {
		b[i] = make([]float64, len(a[0]))
		for j := 0; j < len(a[0]); j++ {
			b[i][j] = a[i][j] - math.Floor(a[i][j]/q)*q
		}
	}
	return b
}

func modZqVec(a []float64, params rlwe.Parameters) []float64 {
	// Components of the vector 'a' belongs to [-q/2, q/2)
	// Takes the modulo operation and maps all components to [0,q)

	q := float64(params.Q()[0])
	b := make([]float64, len(a))
	for i := 0; i < len(a); i++ {
		b[i] = a[i] - math.Floor(a[i]/q)*q
	}
	return b
}

func externalProduct(ctB []*rlwe.Ciphertext, ctA []*rgsw.Ciphertext, evaluatorRGSW *rgsw.Evaluator, ringQ *ring.Ring, params rlwe.Parameters) *rlwe.Ciphertext {
	// Computes the external product between ctA and ctB
	// ctA: n-dimensional RLWE ciphertexts vector
	//// Each element is an RLWE encryption of each elements of vector A
	// ctB: n-dimensional RGSW ciphertexts vector
	//// Each element is an RGSW encryption of each columns of matrix B
	// ctC: 1-dimensional RLWE (packed) ciphertext

	row := len(ctA)
	packedCtC := rlwe.NewCiphertext(params, ctB[0].Degree(), ctB[0].Level())
	tmpCt := rlwe.NewCiphertext(params, ctB[0].Degree(), ctB[0].Level())
	for i := 0; i < row; i++ {
		evaluatorRGSW.ExternalProduct(ctB[i], ctA[i], tmpCt)
		ringQ.Add(packedCtC.Value[0], tmpCt.Value[0], packedCtC.Value[0])
		ringQ.Add(packedCtC.Value[1], tmpCt.Value[1], packedCtC.Value[1])
	}

	return packedCtC
}

func unpackPackedCt(packedCtC *rlwe.Ciphertext, n int, tau int, evaluatorRLWE *rlwe.Evaluator, ringQ *ring.Ring, monomials []ring.Poly, params rlwe.Parameters) []*rlwe.Ciphertext {
	// Unpacks a packed ciphertext and returns an n-dimensional RLWE ciphertexts vector

	scalar := params.Q()[0] - uint64((params.Q()[0]+1)/uint64(tau))
	ctUnpack := make([]*rlwe.Ciphertext, tau)
	ctOut := make([]*rlwe.Ciphertext, n)
	for i := 0; i < tau; i++ {
		ctUnpack[i] = rlwe.NewCiphertext(params, packedCtC.Degree(), packedCtC.Level())
	}
	tmpCt := rlwe.NewCiphertext(params, packedCtC.Degree(), packedCtC.Level())

	// Scaling
	ringQ.MulScalar(packedCtC.Value[0], scalar, packedCtC.Value[0])
	ringQ.MulScalar(packedCtC.Value[1], scalar, packedCtC.Value[1])

	ctUnpack[0] = packedCtC

	for i := tau; i > 1; i /= 2 {
		for j := 0; j < tau; j += i {
			// Automorphism
			evaluatorRLWE.Automorphism(ctUnpack[j], uint64(i+1), tmpCt)

			ringQ.Sub(tmpCt.Value[0], ctUnpack[j].Value[0], ctUnpack[i/2+j].Value[0])
			ringQ.Sub(tmpCt.Value[1], ctUnpack[j].Value[1], ctUnpack[i/2+j].Value[1])

			ringQ.Add(ctUnpack[j].Value[0], tmpCt.Value[0], ctUnpack[j].Value[0])
			ringQ.Add(ctUnpack[j].Value[1], tmpCt.Value[1], ctUnpack[j].Value[1])

			idx := int(math.Log2(float64(i))) - 1
			ringQ.MulCoeffsMontgomery(ctUnpack[i/2+j].Value[0], monomials[idx], ctUnpack[i/2+j].Value[0])
			ringQ.MulCoeffsMontgomery(ctUnpack[i/2+j].Value[1], monomials[idx], ctUnpack[i/2+j].Value[1])
		}
	}

	// Bit reverse
	j := 0
	for i := 1; i < tau; i += 1 {
		bit := tau >> 1
		for j >= bit {
			j -= bit
			bit >>= 1
		}
		j += bit
		if i < j {
			ctUnpack[i], ctUnpack[j] = ctUnpack[j], ctUnpack[i]
		}
	}

	// Takes the first n ciphertexts
	for j := 0; j < n; j += 1 {
		ctOut[j] = ctUnpack[j].CopyNew()
	}

	return ctOut
}

func encryptRlwe(A []float64, scale float64, encryptorRLWE rlwe.Encryptor, ringQ *ring.Ring, params rlwe.Parameters) []*rlwe.Ciphertext {
	// Encrypts an n-dimensional float vector A into an n-dimensional RLWE ciphertexts vector ctA after scaling
	var err error

	row := len(A)
	ctA := make([]*rlwe.Ciphertext, row)
	A_ := utils.ScalarVecMult(scale, A)
	modA := modZqVec(A_, params)
	for r := 0; r < row; r++ {
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		pt.Value.Coeffs[0][0] = uint64(modA[r])
		ringQ.NTT(pt.Value, pt.Value)
		ctA[r], err = encryptorRLWE.EncryptNew(pt)
		if err != nil {
			panic(err)
		}
	}

	return ctA
}

func encryptRgsw(A [][]float64, tau int, encryptorRGSW *rgsw.Encryptor, levelQ int, levelP int, ringQ *ring.Ring, params rlwe.Parameters) []*rgsw.Ciphertext {
	// Encrypts an m-by-n-dimensional float matrix A into an n-dimensional RGSW ciphertexts vector ctA

	row := len(A)
	col := len(A[0])
	ctA := make([]*rgsw.Ciphertext, col)
	modA := modZq(A, params)
	for c := 0; c < col; c++ {
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		for j := 0; j < row; j++ {
			// Store in the packing slots
			pt.Value.Coeffs[0][params.N()*j/tau] = uint64(modA[j][c])
		}
		ringQ.NTT(pt.Value, pt.Value)
		ctA[c] = rgsw.NewCiphertext(params, levelQ, levelP, 0)
		encryptorRGSW.Encrypt(pt, ctA[c])
	}
	return ctA
}

func decryptNewRlwe(ctA []*rlwe.Ciphertext, decryptorRLWE rlwe.Decryptor, scale float64, ringQ *ring.Ring, params rlwe.Parameters) []float64 {
	// 1) Decrypts an n-dimensional RLWE vector ctA and obtain an n-dimensional integer vector pt
	// 2) Maps the constant terms of pt from the set [0,q/2) back to [-q/2, q/2)
	// 3) Scale down and return decA

	row := len(ctA)
	q := float64(params.Q()[0])
	offset := uint64(q / (scale * 2.0))
	decA := make([]float64, row)
	for r := 0; r < row; r++ {
		ringQ.AddScalar(ctA[r].Value[0], offset, ctA[r].Value[0])
		pt := decryptorRLWE.DecryptNew(ctA[r])
		if pt.IsNTT {
			params.RingQ().INTT(pt.Value, pt.Value)
		}
		ringQ.SubScalar(ctA[r].Value[0], offset, ctA[r].Value[0])
		// Constant terms
		val := float64(pt.Value.Coeffs[0][0])
		// Mapping to [-q/2, q/2)
		val = val - math.Floor((val+q/2.0)/q)*q
		// Scale down
		decA[r] = val * scale
	}
	return decA
}

func ctAdd(ctA *rlwe.Ciphertext, ctB *rlwe.Ciphertext, params rlwe.Parameters) *rlwe.Ciphertext {
	// A : m x n
	// B : m x n
	ctC := rlwe.NewCiphertext(params, ctB.Degree(), ctB.Level())

	params.RingQ().Add(ctA.Value[0], ctB.Value[0], ctC.Value[0])
	params.RingQ().Add(ctA.Value[1], ctB.Value[1], ctC.Value[1])

	return ctC
}

func main() {
	params, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:    11,
		LogQ:    []int{56},
		LogP:    []int{47},
		NTTFlag: true,
	})
	fmt.Println("Degree N:", params.N())
	fmt.Println("Ciphertext modulus Q:", params.QBigInt(), "some prime close to 2^54")

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)

	// ======== Compute tau!! ========
	// least power of two greater than n, p_, and m
	tau := 4

	// Generate DFS index
	dfsId := make([]int, tau)
	for i := 0; i < tau; i++ {
		dfsId[i] = i
	}

	tmp := make([]int, tau)
	for i := 1; i < tau; i *= 2 {
		id := 0
		currBlock := tau / i
		nextBlock := currBlock / 2
		for j := 0; j < i; j++ {
			for k := 0; k < nextBlock; k++ {
				tmp[id] = dfsId[j*currBlock+2*k]
				tmp[nextBlock+id] = dfsId[j*currBlock+2*k+1]
				id++
			}
			id += nextBlock
		}

		for j := 0; j < tau; j++ {
			dfsId[j] = tmp[j]
		}
	}

	galEls := make([]uint64, int(math.Log2(float64(tau))))
	for i := 0; i < int(math.Log2(float64(tau))); i++ {
		galEls[i] = uint64(tau/int(math.Pow(2, float64(i))) + 1)
	}

	evkRGSW := rlwe.NewMemEvaluationKeySet(rlk)
	evkRLWE := rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...)

	encryptorRLWE := rlwe.NewEncryptor(params, sk)
	decryptorRLWE := rlwe.NewDecryptor(params, sk)
	encryptorRGSW := rgsw.NewEncryptor(params, sk)
	evaluatorRGSW := rgsw.NewEvaluator(params, evkRGSW)
	evaluatorRLWE := rlwe.NewEvaluator(params, evkRLWE)

	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQ := params.RingQ()

	// ======== Set Scale factors ========
	s := 1 / 10000.0
	L := 1 / 1000.0
	r := 1 / 1000.0

	// ======== Number of iterations ========
	iter := 1000

	// ======== Plant matrices ========
	A := [][]float64{
		{0.9984, 0, 0.0042, 0},
		{0, 0.9989, 0, -0.0033},
		{0, 0, 0.9958, 0},
		{0, 0, 0, 0.9967},
	}
	B := [][]float64{
		{0.0083, 0},
		{0, 0.0063},
		{0, 0.0048},
		{0.0031, 0},
	}
	C := [][]float64{
		{0.5, 0, 0, 0},
		{0, 0.5, 0, 0},
	}

	// ======== Controller matrices ========
	// F: n x n
	// G: n x p
	// H: m x n
	// R: n x m

	F := [][]float64{ // Must be an integer matrix
		{-1, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 2, 0},
		{0, 0, 0, 1},
	}

	G := [][]float64{
		{0.7160, -0.3828},
		{-0.8131, -1.4790},
		{0.6646, 1.1860},
		{0.0181, -0.0060},
	}

	R := [][]float64{
		{-1.7396, 0.3476},
		{0.2588, 1.3226},
		{0.5115, 2.4668},
		{0.0122, 0.0030},
	}

	H := [][]float64{
		{-0.8829, 0.0445, -0.0533, -0.0855},
		{0.1791, 0.2180, -0.2738, 0.0180},
	}

	// ======== Scale up G, R, and H to integers ========
	Gbar := utils.ScalarMatMult(1/s, G)
	Rbar := utils.ScalarMatMult(1/s, R)
	Hbar := utils.ScalarMatMult(1/s, H)

	// ======== Plant and Controller initial state ========
	xPlantInit := []float64{
		1,
		1,
		1,
		1,
	}
	xContInit := []float64{
		0.5,
		0.02,
		-1,
		0.9,
	}

	// ======== F,G,H,R RGSW encrpytion ========
	n := len(F)
	p_ := len(G[0])
	m := len(H)
	fmt.Printf("n \n %d\n", n)
	fmt.Printf("p_ \n %d \n", p_)
	fmt.Printf("m \n %d \n", m)

	// Dimension: 1-by-(# of columns)
	ctF := encryptRgsw(F, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctG := encryptRgsw(Gbar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctH := encryptRgsw(Hbar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctR := encryptRgsw(Rbar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)

	// ======== Run closed-loop without encryption ========
	fmt.Println("Nominal Loop Start")

	// State and output storage
	YOUT := [][]float64{}
	UOUT := [][]float64{}
	XCONT := [][]float64{}
	XPLANT := [][]float64{}

	// State initialization
	xPlantUnenc := xPlantInit
	xContUnenc := xContInit

	startUnenc := time.Now()
	for i := 0; i < iter; i++ {
		// Plant output
		yOut := utils.MatVecMult(C, xPlantUnenc)

		// Controller output
		uOut := utils.MatVecMult(H, xContUnenc)

		// Plant state update
		xPlantUnenc = utils.VecAdd(utils.MatVecMult(A, xPlantUnenc), utils.MatVecMult(B, uOut))

		// Controller state update
		xContUnenc = utils.VecAdd(utils.MatVecMult(F, xContUnenc), utils.MatVecMult(G, yOut))
		xContUnenc = utils.VecAdd(xContUnenc, utils.MatVecMult(R, uOut))

		// Append data
		YOUT = append(YOUT, yOut)
		UOUT = append(UOUT, uOut)
		XCONT = append(XCONT, xContUnenc)
		XPLANT = append(XPLANT, xPlantUnenc)

	}
	elapsedUnenc := time.Now().Sub(startUnenc)

	// ======== Run closed-loop with encryption ========
	fmt.Println("Encrypted Loop Start")

	// State and output storage
	YOUTENC := [][]float64{}
	UOUTENC := [][]float64{}
	XCONTENC := [][]float64{}
	XPLANTENC := [][]float64{}
	TENC := []float64{}

	// State initialization
	// Dimension: 1-by-(# of elements)
	xPlantEnc := xPlantInit
	xContEnc := utils.ScalarVecMult(1/(r*s), xContInit)
	ctxCont := encryptRlwe(xContEnc, 1/L, *encryptorRLWE, ringQ, params)
	ctTmp := rlwe.NewCiphertext(params, ctxCont[0].Degree(), ctxCont[0].Level())

	logn := int(math.Log2(float64(n)))
	monomials := make([]ring.Poly, logn)
	for i := 0; i < logn; i++ {
		monomials[i] = ringQ.NewPoly()
		idx := params.N() - params.N()/(1<<(i+1))
		monomials[i].Coeffs[0][idx] = 1
		ringQ.MForm(monomials[i], monomials[i])
		ringQ.NTT(monomials[i], monomials[i])
	}

	startEnc := time.Now()
	for i := 0; i < iter; i++ {
		startEncIter := time.Now()
		// Plant output
		yOut := utils.MatVecMult(C, xPlantEnc)

		// Quantize and encrypt plant output
		yOutRound := utils.RoundVec(utils.ScalarVecMult(1/r, yOut))
		ctyOut := encryptRlwe(yOutRound, 1/L, *encryptorRLWE, ringQ, params)

		// Controller output
		packedCtuOut := externalProduct(ctxCont, ctH, evaluatorRGSW, ringQ, params)

		// Unpack controller output
		unpackedCtuOut := unpackPackedCt(packedCtuOut, m, tau, evaluatorRLWE, ringQ, monomials, params)

		// Decrypt controller output and construct plant input
		uOut := decryptNewRlwe(unpackedCtuOut, *decryptorRLWE, r*s*s*L, ringQ, params)

		// Re-encrypt controller output
		ctuReEnc := encryptRlwe(uOut, 1/(r*L), *encryptorRLWE, ringQ, params)

		// Plant state update
		xPlantEnc = utils.VecAdd(utils.MatVecMult(A, xPlantEnc), utils.MatVecMult(B, uOut))

		// Controller state update
		ctFx := externalProduct(ctxCont, ctF, evaluatorRGSW, ringQ, params)
		ctGy := externalProduct(ctyOut, ctG, evaluatorRGSW, ringQ, params)
		ctRu := externalProduct(ctuReEnc, ctR, evaluatorRGSW, ringQ, params)
		ctTmp = ctAdd(ctFx, ctGy, params)
		ctTmp = ctAdd(ctTmp, ctRu, params)
		ctxCont = unpackPackedCt(ctTmp, n, tau, evaluatorRLWE, ringQ, monomials, params)

		elapsedEncIter := time.Now().Sub(startEncIter)
		// Decrypt plant output just for validation
		valyOut := decryptNewRlwe(ctyOut, *decryptorRLWE, r*L, ringQ, params)

		// Decrypt controller state just for validation
		valxCont := decryptNewRlwe(ctxCont, *decryptorRLWE, r*s*L, ringQ, params)

		// Append data
		YOUTENC = append(YOUTENC, valyOut)
		UOUTENC = append(UOUTENC, uOut)
		XCONTENC = append(XCONTENC, valxCont)
		XPLANTENC = append(XPLANTENC, xPlantEnc)
		TENC = append(TENC, float64(elapsedEncIter.Microseconds()))
	}
	elapsedEnc := time.Now().Sub(startEnc)

	// ======== For debugging ========
	// for i := 0; i < iter; i++ {
	// 	fmt.Println("=======================")
	// 	fmt.Printf("Plant output at iter %d: \n %v \n", i, YOUT[i])
	// 	fmt.Printf("Controller output at iter %d: \n %v \n", i, UOUT[i])
	// 	fmt.Printf("Plant state at iter %d: \n %v \n", i, XPLANT[i])
	// 	fmt.Printf("Controller state at iter %d: \n %v \n", i, XCONT[i])
	// 	// fmt.Println("=======================")
	// 	// fmt.Printf("Decrypted Plant output at iter %d: \n %v \n", i, YOUTENC[i])
	// 	// fmt.Printf("Decrypted controller output at iter %d: \n %v \n", i, UOUTENC[i])
	// 	// fmt.Printf("Decrypted Plant state at iter %d: \n %v \n", i, XPLANTENC[i])
	// 	// fmt.Printf("Decrypted controller state at iter %d: \n %v \n", i, XCONTENC[i])
	// }

	// ======== Simulation result ========
	fmt.Println("Iterations: ", iter)
	fmt.Println("Unenc total time: ", elapsedUnenc)
	fmt.Println("Enc total time: ", elapsedEnc)
	fmt.Printf("Unenc average time for one iteration: %v us \n", float64(elapsedUnenc.Microseconds())/float64(iter))
	fmt.Printf("Enc average time for one iteration: %v us \n", float64(elapsedEnc.Microseconds())/float64(iter))

	// **************************** Plot section *************************************
	// Plot 2-norm of the difference between plant state, controller state, plant output, and controller output
	rows := 2
	cols := 2
	plotsDiff := make([][]*plot.Plot, rows)
	for i := 0; i < rows; i++ {
		plotsDiff[i] = make([]*plot.Plot, cols)
		for j := 0; j < cols; j++ {
			p := plot.New()

			pts := make(plotter.XYs, iter)

			for k := range pts {
				pts[k].X = float64(k)
			}

			// First row
			if i == 0 {
				if j == 0 {
					for k := range pts {
						pts[k].Y = utils.Vec2Norm(utils.VecSub(XPLANT[k], XPLANTENC[k]))

					}
					p.Title.Text = "Plant State Difference"
				} else if j == 1 {
					for k := range pts {
						pts[k].Y = utils.Vec2Norm(utils.VecSub(YOUT[k], YOUTENC[k]))
					}
					p.Title.Text = "Plant Output Difference"
				}
			} else {

				// Second row
				if j == 0 {
					for k := range pts {
						pts[k].Y = utils.Vec2Norm(utils.VecSub(XCONT[k], XCONTENC[k]))

					}
					p.Title.Text = "Controller State Difference"
				} else if j == 1 {
					for k := range pts {
						pts[k].Y = utils.Vec2Norm(utils.VecSub(UOUT[k], UOUTENC[k]))
					}
					p.Title.Text = "Controller Output Difference"
				}
			}
			p.X.Label.Text = "iteration"
			p.Y.Label.Text = "2-norm value"
			// p.Y.Min = -1e-04
			// p.Y.Max = 1e-04
			p.Add(plotter.NewGrid())
			lLine, lPoints, _ := plotter.NewLinePoints(pts)
			p.Add(lLine, lPoints)
			plotsDiff[i][j] = p

		}
	}

	imgDiff := vgimg.New(vg.Points(1000), vg.Points(500))
	dcDiff := draw.New(imgDiff)

	t := draw.Tiles{
		Rows:      rows,
		Cols:      cols,
		PadX:      vg.Millimeter,
		PadY:      vg.Millimeter,
		PadTop:    vg.Points(50),
		PadBottom: vg.Points(50),
		PadLeft:   vg.Points(50),
		PadRight:  vg.Points(50),
	}

	canvasesDiff := plot.Align(plotsDiff, t, dcDiff)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			if plotsDiff[i][j] != nil {
				plotsDiff[i][j].Draw(canvasesDiff[i][j])
			}
		}
	}

	w, err := os.Create("Difference.png")
	if err != nil {
		panic(err)
	}
	defer w.Close()
	png := vgimg.PngCanvas{Canvas: imgDiff}
	if _, err := png.WriteTo(w); err != nil {
		panic(err)
	}

	// Export data ===============================================================

	fileUOUT, err := os.Create("./uUnenc.csv")
	if err != nil {
		panic(err)
	}
	wrUOUT := csv.NewWriter(bufio.NewWriter(fileUOUT))
	UOUTstr := utils.MatToString(UOUT)
	wrUOUT.WriteAll(UOUTstr)

	fileUOUTENC, err := os.Create("./uEnc.csv")
	if err != nil {
		panic(err)
	}
	wrUOUTENC := csv.NewWriter(bufio.NewWriter(fileUOUTENC))
	UOUTENCstr := utils.MatToString(UOUTENC)
	wrUOUTENC.WriteAll(UOUTENCstr)

	fileTENC, err := os.Create("./TEncPackLarge.csv")
	if err != nil {
		panic(err)
	}
	wrTENC := csv.NewWriter(bufio.NewWriter(fileTENC))
	TENCstr := utils.VecToString(TENC)
	wrTENC.WriteAll(TENCstr)
}
