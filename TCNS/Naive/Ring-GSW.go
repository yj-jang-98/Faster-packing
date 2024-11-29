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

	modA := modZqVec(A_, params)
	for r := 0; r < row; r++ {
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		for i := 0; i < params.N(); i++ {
			pt.Value.Coeffs[0][i] = uint64(modA[r])
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
	modA := modZq(A, params)
	for r := 0; r < row; r++ {
		ctA[r] = make([]*rgsw.Ciphertext, col)
		for c := 0; c < col; c++ {
			pt := rlwe.NewPlaintext(params, params.MaxLevel())
			for i := 0; i < params.N(); i++ {
				pt.Value.Coeffs[0][i] = uint64(modA[r][c])
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

func main() {
	params, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:    12,
		LogQ:    []int{56},
		LogP:    []int{42},
		NTTFlag: true,
	})
	fmt.Println("Degree N:", params.N())
	fmt.Println("Ciphertext modulus Q:", params.QBigInt(), "some prime close to 2^54")

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	encryptorRLWE := rlwe.NewEncryptor(params, sk)
	decryptorRLWE := rlwe.NewDecryptor(params, sk)
	encryptorRGSW := rgsw.NewEncryptor(params, sk)
	evaluator := rgsw.NewEvaluator(params, evk)

	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	// decompRNS := params.BaseRNSDecompositionVectorSize(levelQ, levelP)
	decompPw2 := params.BaseTwoDecompositionVectorSize(levelQ, levelP, 0)
	// ringQP := params.RingQP()
	ringQ := params.RingQ()

	fmt.Println("Ciphertext modulus:", params.QBigInt())
	fmt.Println("Degree of polynomials:", params.N())
	fmt.Println("base 2 decomposition:", decompPw2)

	// ======== Set Scale factors ========
	s := 1 / 10000.0
	L := 1 / 1000.0
	r := 1 / 1000.0

	// ======== Number of iterations ========
	iter := 500

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

	ctF := encryptRgsw(F, encryptorRGSW, levelQ, levelP, params)
	ctG := encryptRgsw(Gbar, encryptorRGSW, levelQ, levelP, params)
	ctH := encryptRgsw(Hbar, encryptorRGSW, levelQ, levelP, params)
	ctR := encryptRgsw(Rbar, encryptorRGSW, levelQ, levelP, params)

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
	xPlantEnc := xPlantInit
	xContEnc := utils.ScalarVecMult(1/(r*s), xContInit)
	ctxCont := encryptRlwe(xContEnc, 1/L, *encryptorRLWE, params)

	startEnc := time.Now()
	for i := 0; i < iter; i++ {
		startEncIter := time.Now()
		// Plant output
		yOut := utils.MatVecMult(C, xPlantEnc)

		// Quantize and encrypt plant output
		yOutRound := utils.RoundVec(utils.ScalarVecMult(1/r, yOut))
		ctyOut := encryptRlwe(yOutRound, 1/L, *encryptorRLWE, params)

		// Controller output
		ctuOut := externalProduct(ctxCont, ctH, evaluator, ringQ, params)

		// Decrypt controller output and construct plant input
		uOut := decryptRlwe(ctuOut, *decryptorRLWE, r*s*s*L, params)

		// Re-encrypt controller output
		ctuReEnc := encryptRlwe(uOut, 1/(r*L), *encryptorRLWE, params)

		// Plant state update
		xPlantEnc = utils.VecAdd(utils.MatVecMult(A, xPlantEnc), utils.MatVecMult(B, uOut))

		// Controller state update
		ctFx := externalProduct(ctxCont, ctF, evaluator, ringQ, params)
		ctGy := externalProduct(ctyOut, ctG, evaluator, ringQ, params)
		ctRu := externalProduct(ctuReEnc, ctR, evaluator, ringQ, params)
		ctxCont = ctAdd(ctFx, ctGy, params)
		ctxCont = ctAdd(ctxCont, ctRu, params)

		elapsedEncIter := time.Now().Sub(startEncIter)

		// Decrypt controller state just for validation
		valxCont := decryptRlwe(ctxCont, *decryptorRLWE, r*s*L, params)

		// Decrypt plant output just for validation
		valyOut := decryptRlwe(ctyOut, *decryptorRLWE, r*L, params)

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
	// 	fmt.Println("=======================")
	// 	fmt.Printf("Decrypted Plant output at iter %d: \n %v \n", i, YOUTENC[i])
	// 	fmt.Printf("Decrypted controller output at iter %d: \n %v \n", i, UOUTENC[i])
	// 	fmt.Printf("Decrypted Plant state at iter %d: \n %v \n", i, XPLANTENC[i])
	// 	fmt.Printf("Decrypted controller state at iter %d: \n %v \n", i, XCONTENC[i])
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

	fileTENC, err := os.Create("./TEncNaiveLarge.csv")
	if err != nil {
		panic(err)
	}
	wrTENC := csv.NewWriter(bufio.NewWriter(fileTENC))
	TENCstr := utils.VecToString(TENC)
	wrTENC.WriteAll(TENCstr)
}
